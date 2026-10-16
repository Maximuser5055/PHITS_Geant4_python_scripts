//
// TETPrimaryGeneratorAction.cc
// Primary particles are allocated among selected source organ IDs in
// proportion to their phantom-geometry masses. Within each source ID,
// tetrahedra are selected proportional to volume and positions are generated
// uniformly using barycentric coordinates.
//

#include "TETPrimaryGeneratorAction.hh"

#include "G4Exception.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4Threading.hh"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

TETPrimaryGeneratorAction::TETPrimaryGeneratorAction(
    TETModelImport* tetData,
    const std::vector<G4int>& internalSources,
    const G4String& outputFilePath)
    : G4VUserPrimaryGeneratorAction(),
      fParticleGun(nullptr),
      fLastPreparedRunID(-1),
      fPreparedEventCount(-1),
      fOutputFilePath(outputFilePath)
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
        region.totalVolume = 0.0;
        region.allocatedParticles = 0;
        region.cumulativeParticleEnd = 0;
        region.fractionalRemainder = 0.0;

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

    // Collect tetrahedra for each requested source ID. Cache each tetrahedron's
    // vertices and cumulative volume so source positions can be generated
    // directly, without bounding-box rejection or repeated Inside() checks.
    for (G4int i = 0; i < tetData->GetNumTetrahedron(); ++i) {
        const G4int organID = tetData->GetMaterialIndex(i);
        const auto indexIt = fSourceIndexByID.find(organID);
        if (indexIt == fSourceIndexByID.end()) continue;

        SourceRegion& region = fSourceRegions[indexIt->second];
        G4Tet* tetSolid = tetData->GetTetrahedron(i);

        TetrahedronData tetDataForSampling;
        tetSolid->GetVertices(
            tetDataForSampling.vertices[0],
            tetDataForSampling.vertices[1],
            tetDataForSampling.vertices[2],
            tetDataForSampling.vertices[3]);

        const G4double volume = tetSolid->GetCubicVolume();
        if (!(volume > 0.0)) {
            continue; // Ignore degenerate or invalid zero-volume tetrahedra.
        }

        region.totalVolume += volume;
        tetDataForSampling.cumulativeVolume = region.totalVolume;
        region.tetrahedra.push_back(tetDataForSampling);
    }

    for (const auto& region : fSourceRegions) {
        if (region.tetrahedra.empty() || !(region.totalVolume > 0.0)) {
            std::ostringstream message;
            message << "There are no positive-volume tetrahedra for requested source organ ID "
                    << region.organID << ".";
            G4Exception("TETPrimaryGeneratorAction::TETPrimaryGeneratorAction", "",
                        FatalErrorInArgument, message.str().c_str());
            return;
        }
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

    // Export the planned allocation to a text file instead of the terminal.
    // In multithreaded builds, only worker thread 0 writes the report.
#ifdef G4MULTITHREADED
    if (G4Threading::G4GetThreadId() == 0) {
#endif
        namespace fs = std::filesystem;
        const fs::path outputPath(fOutputFilePath.c_str());
        const fs::path outputDirectory = outputPath.has_parent_path()
                                             ? outputPath.parent_path()
                                             : fs::current_path();
        const fs::path logPath = outputDirectory / "source_particle_distribution.txt";
        std::ofstream logFile(logPath, std::ios::out | std::ios::app);

        if (logFile.is_open()) {
            logFile << "========================================\n";
            logFile << "Source Particle Allocation\n";
            logFile << "Run ID: " << runID << "\n";
            logFile << "Total requested events: " << nEvents << "\n";
            logFile << "----------------------------------------\n";

            G4int reportedTotal = 0;
            for (const auto& region : fSourceRegions) {
                const G4double fraction =
                    (nEvents > 0)
                        ? static_cast<G4double>(region.allocatedParticles) / nEvents
                        : 0.0;

                logFile << "Organ ID " << region.organID
                        << ": " << region.allocatedParticles
                        << " particles (" << std::fixed
                        << std::setprecision(4) << 100.0 * fraction << "%)"
                        << " | mass = " << std::setprecision(10)
                        << region.mass / g << " g\n";
                reportedTotal += region.allocatedParticles;
            }

            logFile << "----------------------------------------\n";
            logFile << "Total allocated: " << reportedTotal << "\n";
            if (reportedTotal != nEvents) {
                logFile << "WARNING: allocated count does not equal requested event count!\n";
            }
            logFile << "========================================\n\n";
            logFile.close();
        }
#ifdef G4MULTITHREADED
    }
#endif

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

    // Select a tetrahedron with probability proportional to its volume.
    // The cached cumulative volumes are sorted, so binary search selects the
    // tetrahedron without rejection sampling or scanning the full mesh.
    const G4double targetVolume =
        G4UniformRand() * selectedRegion->totalVolume;

    const auto tetIt = std::lower_bound(
        selectedRegion->tetrahedra.begin(),
        selectedRegion->tetrahedra.end(),
        targetVolume,
        [](const TetrahedronData& tet, G4double target) {
            return tet.cumulativeVolume < target;
        });

    const TetrahedronData& selectedTet =
        (tetIt == selectedRegion->tetrahedra.end())
            ? selectedRegion->tetrahedra.back()
            : *tetIt;

    // Draw four independent exponential variates and normalize them. The
    // resulting barycentric weights are uniform on the tetrahedral simplex.
    G4double weights[4];
    G4double weightSum = 0.0;
    for (G4int i = 0; i < 4; ++i) {
        G4double u = 0.0;
        do {
            u = G4UniformRand();
        } while (u <= 0.0);

        weights[i] = -std::log(u);
        weightSum += weights[i];
    }

    G4ThreeVector pos(0., 0., 0.);
    for (G4int i = 0; i < 4; ++i) {
        weights[i] /= weightSum;
        pos += weights[i] * selectedTet.vertices[i];
    }

    fParticleGun->SetParticlePosition(pos);
    fParticleGun->SetParticleMomentumDirection(G4RandomDirection());
    fParticleGun->GeneratePrimaryVertex(anEvent);
}
