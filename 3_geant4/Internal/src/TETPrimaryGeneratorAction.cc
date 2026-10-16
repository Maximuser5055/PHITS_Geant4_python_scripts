//
// TETPrimaryGeneratorAction.cc
// Primary particles are allocated among selected source organ IDs in
// proportion to their phantom-geometry masses. Each source ID has its own
// bounding box and rejection sampler.
//

#include "TETPrimaryGeneratorAction.hh"

#include "G4Exception.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

TETPrimaryGeneratorAction::TETPrimaryGeneratorAction(
    TETModelImport* tetData,
    const std::vector<G4int>& internalSources)
    : G4VUserPrimaryGeneratorAction(),
      fParticleGun(nullptr),
      fLastPreparedRunID(-1),
      fPreparedEventCount(-1)
{
    if (internalSources.empty()) {
        G4Exception("TETPrimaryGeneratorAction::TETPrimaryGeneratorAction", "",
                    FatalErrorInArgument,
                    "No source organ IDs were provided.");
        return;
    }	

    // Create one source-region record for each requested organ ID. Keep the
    // caller's order for deterministic tie-breaking in largest-remainder
    // allocation.
    for (const G4int organID : internalSources) {
        if (fSourceIndexByID.find(organID) != fSourceIndexByID.end()) {
            continue; // Ignore duplicate IDs in the command line.
        }

        SourceRegion region;
        region.organID = organID;
        region.mass = 0.0;
        region.bBoxMin = G4ThreeVector(0., 0., 0.);
        region.bBoxDim = G4ThreeVector(0., 0., 0.);
        region.allocatedParticles = 0;
        region.cumulativeParticleEnd = 0;
        region.fractionalRemainder = 0.0;
        region.xMin = region.yMin = region.zMin = DBL_MAX;
        region.xMax = region.yMax = region.zMax = -DBL_MAX;

        fSourceIndexByID[organID] = fSourceRegions.size();
        fSourceRegions.push_back(region);
    }

    // Geometry mass map is the source of the mass fractions. The absolute
    // mass units cancel in m_i / sum(m_i).
    const auto& massMap = tetData->GetMassMap();
    for (auto& region : fSourceRegions) {
        const auto massIt = massMap.find(region.organID);
        if (massIt == massMap.end() || massIt->second <= 0.0) {
            std::ostringstream message;
            message << "No positive phantom-geometry mass was found for source organ ID "
                    << region.organID << ".";
            G4Exception("TETPrimaryGeneratorAction::TETPrimaryGeneratorAction", "",
                        FatalErrorInArgument, message.str().c_str());
            return;
        }
        region.mass = massIt->second;
    }

    // Collect tetrahedra and compute an independent bounding box for each ID.
    for (G4int i = 0; i < tetData->GetNumTetrahedron(); ++i) {
        const G4int organID = tetData->GetMaterialIndex(i);
        const auto indexIt = fSourceIndexByID.find(organID);
        if (indexIt == fSourceIndexByID.end()) continue;

        SourceRegion& region = fSourceRegions[indexIt->second];
        G4Tet* tetSolid = tetData->GetTetrahedron(i);
        region.tetrahedra.push_back(tetSolid);

        for (const auto& vertex : tetSolid->GetVertices()) {
            region.xMin = std::min(region.xMin, vertex.getX());
            region.xMax = std::max(region.xMax, vertex.getX());
            region.yMin = std::min(region.yMin, vertex.getY());
            region.yMax = std::max(region.yMax, vertex.getY());
            region.zMin = std::min(region.zMin, vertex.getZ());
            region.zMax = std::max(region.zMax, vertex.getZ());
        }
    }

    for (auto& region : fSourceRegions) {
        if (region.tetrahedra.empty()) {
            std::ostringstream message;
            message << "There are no tetrahedra for requested source organ ID "
                    << region.organID << ".";
            G4Exception("TETPrimaryGeneratorAction::TETPrimaryGeneratorAction", "",
                        FatalErrorInArgument, message.str().c_str());
            return;
        }
        region.bBoxMin = G4ThreeVector(region.xMin, region.yMin, region.zMin);
        region.bBoxDim = G4ThreeVector(region.xMax - region.xMin,
                                       region.yMax - region.yMin,
                                       region.zMax - region.zMin);
    }

    // Initialise the gamma particle gun. Energy is still set by the existing
    // macro/UI configuration if that is how the application is run.
    fParticleGun = new G4ParticleGun(1);
    G4ParticleDefinition* particle =
        G4ParticleTable::GetParticleTable()->FindParticle("gamma");
    fParticleGun->SetParticleDefinition(particle);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1., 0., 0.));
    fParticleGun->SetParticleEnergy(1. * MeV);
}

TETPrimaryGeneratorAction::~TETPrimaryGeneratorAction()
{
    delete fParticleGun;
}

void TETPrimaryGeneratorAction::PrepareAllocationForCurrentRun()
{
    G4RunManager* runManager = G4RunManager::GetRunManager();
    const G4Run* currentRun = runManager ? runManager->GetCurrentRun() : nullptr;
    if (!currentRun) {
        G4Exception("TETPrimaryGeneratorAction::PrepareAllocationForCurrentRun", "",
                    FatalException,
                    "Cannot allocate source particles because there is no current Geant4 run.");
        return;
    }

    const G4int runID = currentRun->GetRunID();
    const G4int nEvents = currentRun->GetNumberOfEventToBeProcessed();
    if (runID == fLastPreparedRunID && nEvents == fPreparedEventCount) return;

    if (nEvents < 0) {
        G4Exception("TETPrimaryGeneratorAction::PrepareAllocationForCurrentRun", "",
                    FatalException,
                    "Geant4 returned a negative number of events for the current run.");
        return;
    }

    G4double totalMass = 0.0;
    for (const auto& region : fSourceRegions) totalMass += region.mass;
    if (!(totalMass > 0.0)) {
        G4Exception("TETPrimaryGeneratorAction::PrepareAllocationForCurrentRun", "",
                    FatalException, "The total mass of the selected source regions is not positive.");
        return;
    }

    // Largest-remainder (Hamilton) apportionment: floor each ideal allocation,
    // then distribute the remaining events to the largest fractional parts.
    G4int allocated = 0;
    for (auto& region : fSourceRegions) {
        const G4double ideal = static_cast<G4double>(nEvents) * region.mass / totalMass;
        region.allocatedParticles = static_cast<G4int>(std::floor(ideal));
        region.fractionalRemainder = ideal - region.allocatedParticles;
        allocated += region.allocatedParticles;
    }

    G4int remaining = nEvents - allocated;
    std::vector<std::size_t> remainderOrder(fSourceRegions.size());
    std::iota(remainderOrder.begin(), remainderOrder.end(), 0);
    std::stable_sort(remainderOrder.begin(), remainderOrder.end(),
        [this](std::size_t a, std::size_t b) {
            return fSourceRegions[a].fractionalRemainder >
                   fSourceRegions[b].fractionalRemainder;
        });

    for (G4int i = 0; i < remaining; ++i) {
        ++fSourceRegions[remainderOrder[static_cast<std::size_t>(i) % remainderOrder.size()]]
              .allocatedParticles;
    }

    // Convert per-region counts to non-overlapping event-ID intervals. Event
    // IDs are assigned globally by Geant4, so this mapping remains exact and
    // race-free with multithreaded event processing.
    G4int cumulativeEnd = 0;
    for (auto& region : fSourceRegions) {
        cumulativeEnd += region.allocatedParticles;
        region.cumulativeParticleEnd = cumulativeEnd;
    }

    fLastPreparedRunID = runID;
    fPreparedEventCount = nEvents;

}

void TETPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    PrepareAllocationForCurrentRun();

    const G4int eventID = anEvent->GetEventID();
    const SourceRegion* selectedRegion = nullptr;
    for (const auto& region : fSourceRegions) {
        if (eventID < region.cumulativeParticleEnd) {
            selectedRegion = &region;
            break;
        }
    }

    if (!selectedRegion) {
        G4Exception("TETPrimaryGeneratorAction::GeneratePrimaries", "",
                    FatalException,
                    "Could not map this event ID to a source organ. Check run event counts.");
        return;
    }

    // Rejection-sample only within this organ ID's own bounding box and
    // tetrahedra. This avoids sampling the empty space between distant IDs.
    G4bool insideChk = false;
    G4ThreeVector pos;
    do {
        pos = G4ThreeVector(
            selectedRegion->bBoxMin.getX() + selectedRegion->bBoxDim.getX() * G4UniformRand(),
            selectedRegion->bBoxMin.getY() + selectedRegion->bBoxDim.getY() * G4UniformRand(),
            selectedRegion->bBoxMin.getZ() + selectedRegion->bBoxDim.getZ() * G4UniformRand());

        for (const auto* tet : selectedRegion->tetrahedra) {
            if (tet->Inside(pos) == kOutside) continue;
            insideChk = true;
            break;
        }
    } while (!insideChk);

    fParticleGun->SetParticlePosition(pos);
    fParticleGun->SetParticleMomentumDirection(G4RandomDirection());
    fParticleGun->GeneratePrimaryVertex(anEvent);
}
