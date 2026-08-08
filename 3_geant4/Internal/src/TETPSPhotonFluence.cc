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
TETPSPhotonFluence::energyBins;


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
//
// 0.001 MeV = 1 keV
// 1.0   MeV
//
// 101 points gives 100 intervals.
// ============================================================

void TETPSPhotonFluence::InitializeEnergyBins()
{
    if (!energyBins.empty())
        return;

    const G4int nBins = 101;

    const G4double Emin = 0.001;  // MeV
    const G4double Emax = 1.0;    // MeV

    energyBins.resize(nBins);

    for (G4int i = 0; i < nBins; i++)
    {
        G4double fraction =
            static_cast<G4double>(i)
            / static_cast<G4double>(nBins - 1);

        energyBins[i] =
            Emin *
            std::pow(
                Emax / Emin,
                fraction
            );
    }
}


// ============================================================
// Return representative energy of a bin
// ============================================================

G4double TETPSPhotonFluence::GetEnergy(
    G4int bin
)
{
    if (bin < 0 ||
        bin >= static_cast<G4int>(energyBins.size()))
    {
        return -1.0;
    }

    return energyBins[bin];
}


// ============================================================
// Initialize
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

G4bool TETPSPhotonFluence::ProcessHits(
    G4Step* aStep,
    G4TouchableHistory*
)
{
    // --------------------------------------------------------
    // Only score photons
    // --------------------------------------------------------

    G4Track* track =
        aStep->GetTrack();

    if (track->GetDefinition()
        != G4Gamma::GammaDefinition())
    {
        return false;
    }


    // --------------------------------------------------------
    // Photon kinetic energy BEFORE the step
    // --------------------------------------------------------

    G4double energy =
        aStep->GetPreStepPoint()
            ->GetKineticEnergy();


    // Only energies within our range
    if (energy < 0.001 * MeV ||
        energy > 1.0 * MeV)
    {
        return false;
    }


    // --------------------------------------------------------
    // Find energy bin
    // --------------------------------------------------------

    G4double energyMeV =
        energy / MeV;

    auto upper =
        std::lower_bound(
            energyBins.begin(),
            energyBins.end(),
            energyMeV
        );

    G4int bin =
        upper - energyBins.begin();

    if (bin >= static_cast<G4int>(
                   energyBins.size()))
    {
        bin =
            energyBins.size() - 1;
    }


    // --------------------------------------------------------
    // Find tetrahedron
    // --------------------------------------------------------

    G4int copyNo =
        aStep->GetPreStepPoint()
            ->GetTouchable()
            ->GetCopyNumber();


    // Convert tetrahedron index
    // to organ/material ID
    G4int organID =
        tetData->GetMaterialIndex(copyNo);


    // --------------------------------------------------------
    // Track length
    // --------------------------------------------------------

    G4double trackLength =
        aStep->GetStepLength();


    if (trackLength <= 0.)
        return false;


    // --------------------------------------------------------
    // Encode:
    //
    // key = organID * number_of_bins + energy_bin
    //
    // Example:
    //
    // organ 1400, bin 20
    //
    // key = 1400 * 101 + 20
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