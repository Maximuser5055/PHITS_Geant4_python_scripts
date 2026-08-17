#include "TETPSPhotonFluence.hh"

#include "G4Gamma.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4THitsMap.hh"
#include "G4MultiFunctionalDetector.hh"

#include <vector>
#include <cmath>
#include <algorithm>


// ============================================================
// Static energy-bin storage
// ============================================================

std::vector<G4double>
TETPSPhotonFluence::energyEdges;


// ============================================================
// Constructor
// ============================================================

TETPSPhotonFluence::TETPSPhotonFluence(
    G4String name,
    TETModelImport* _tetData
)
    : G4VPrimitiveScorer(name),
      tetData(_tetData),
      HCID(-1),
      evtMap(nullptr)
{
    InitializeEnergyBins();
}


// ============================================================
// Destructor
// ============================================================

TETPSPhotonFluence::~TETPSPhotonFluence()
{
}


// ============================================================
// Create logarithmic energy bins
// ============================================================

void TETPSPhotonFluence::InitializeEnergyBins()
{
    if (!energyEdges.empty())
        return;

    energyEdges.resize(nEnergyBins + 1);

    for (G4int i = 0; i <= nEnergyBins; i++)
    {
        G4double fraction =
            static_cast<G4double>(i)
            / static_cast<G4double>(nEnergyBins);

        energyEdges[i] =
            Emin *
            std::pow(
                Emax / Emin,
                fraction
            );
    }
}


// ============================================================
// Lower edge of energy bin
// ============================================================

G4double TETPSPhotonFluence::GetEnergyLow(
    G4int bin
)
{
    if (
        bin < 0 ||
        bin >= nEnergyBins
    )
    {
        return -1.0;
    }


    return energyEdges[bin];
}


// ============================================================
// Upper edge of energy bin
// ============================================================

G4double TETPSPhotonFluence::GetEnergyHigh(
    G4int bin
)
{
    if (
        bin < 0 ||
        bin >= nEnergyBins
    )
    {
        return -1.0;
    }


    return energyEdges[bin + 1];
}


// ============================================================
// Geometric center of energy bin
//
// Because bins are logarithmic, the geometric mean is used.
// ============================================================

G4double TETPSPhotonFluence::GetEnergyCenter(
    G4int bin
)
{
    if (
        bin < 0 ||
        bin >= nEnergyBins
    )
    {
        return -1.0;
    }


    return std::sqrt(
        energyEdges[bin]
        *
        energyEdges[bin + 1]
    );
}


// ============================================================
// Initialize scorer
// ============================================================

void TETPSPhotonFluence::Initialize(
    G4HCofThisEvent* HCE
)
{
    evtMap =
        new G4THitsMap<G4double>(
            GetMultiFunctionalDetector()->GetName(),
            GetName()
        );


    if (HCID < 0)
    {
        HCID =
            G4SDManager::GetSDMpointer()
                ->GetCollectionID(
                    GetMultiFunctionalDetector()
                        ->GetName()
                    + "/" + GetName()
                );
    }


    HCE->AddHitsCollection(
        HCID,
        evtMap
    );
}


// ============================================================
// Process photon steps
// ============================================================

G4bool TETPSPhotonFluence::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
    // --------------------------------------------------------
    // Only score photons
    // --------------------------------------------------------

    G4Track* track =
        aStep->GetTrack();


    if (
        track->GetDefinition()!= G4Gamma::GammaDefinition()
    )
    {
        return false;
    }

    // --------------------------------------------------------
    // Photon kinetic energy before the step
    // --------------------------------------------------------

    G4double energy =
        aStep->GetPreStepPoint()
            ->GetKineticEnergy();


    G4double energyMeV =
        energy / MeV;


    // --------------------------------------------------------
    // Only score between energy min and max
    // --------------------------------------------------------

    if (
        energyMeV < Emin ||
        energyMeV > Emax
    )
    {
        return false;
    }


    // --------------------------------------------------------
    // Find energy bin
    //
    // Find first upper edge > energy
    // --------------------------------------------------------

    auto upper =
        std::upper_bound(
            energyEdges.begin(),
            energyEdges.end(),
            energyMeV
        );


    G4int bin =
        static_cast<G4int>(
            upper - energyEdges.begin()
        ) - 1;


    // Safety check

    if (bin < 0)
    {
        bin = 0;
    }


    if (bin >= nEnergyBins)
    {
        bin = nEnergyBins - 1;
    }


    // --------------------------------------------------------
    // Find tetrahedron copy number
    // --------------------------------------------------------

    G4int copyNo =
        aStep->GetPreStepPoint()
            ->GetTouchable()
            ->GetCopyNumber();


    // --------------------------------------------------------
    // Convert tetrahedron index
    // to organ/material ID
    // --------------------------------------------------------

    G4int organID =
        tetData->GetMaterialIndex(
            copyNo
        );


    // --------------------------------------------------------
    // Get photon track length
    // --------------------------------------------------------

    G4double trackLength =
        aStep->GetStepLength();


    if (trackLength <= 0.0)
    {
        return false;
    }


    // --------------------------------------------------------
    // Encode organ ID + energy bin
    //
    // key = organID * nEnergyBins + bin
    // --------------------------------------------------------

    G4int key =
        organID * nEnergyBins
        + bin;


    // --------------------------------------------------------
    // Add track length
    // --------------------------------------------------------

    evtMap->add(
        key,
        trackLength
    );


    return true;
}


// ============================================================
// End of event
// ============================================================

void TETPSPhotonFluence::EndOfEvent(
    G4HCofThisEvent*
)
{
}


// ============================================================
// Clear
// ============================================================

void TETPSPhotonFluence::clear()
{
}


// ============================================================
// Get hits map
// ============================================================

G4THitsMap<G4double>*
TETPSPhotonFluence::GetHitsMap()
{
    return evtMap;
}