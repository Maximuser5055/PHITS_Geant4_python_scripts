#ifndef TETPSPhotonFluence_h
#define TETPSPhotonFluence_h 1

#include "G4VPrimitiveScorer.hh"
#include "G4THitsMap.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"

#include "TETModelImport.hh"

#include <vector>


class TETPSPhotonFluence
    : public G4VPrimitiveScorer
{
public:

    TETPSPhotonFluence(
        G4String name,
        TETModelImport* tetData
    );

    virtual ~TETPSPhotonFluence();


    virtual G4bool ProcessHits(
        G4Step* aStep,
        G4TouchableHistory* ROhist
    ) override;


    virtual void Initialize(
        G4HCofThisEvent* HCE
    ) override;


    virtual void clear() override;


    virtual void EndOfEvent(
        G4HCofThisEvent* HCE
    ) override;


    virtual void DrawAll() override
    {
    }


    virtual void PrintAll() override
    {
    }

    G4THitsMap<G4double>*
    GetHitsMap();


    // ========================================================
    // Enable fluence for photons
    // Number of energy bins, energy min and max
    // ========================================================

    static constexpr G4bool ENABLE_PHOTON_FLUENCE = false;
    static constexpr G4int nEnergyBins = 100;
    static constexpr G4double Emin = 0.01;
    static constexpr G4double Emax = 10;

    // ========================================================
    // Get energy-bin boundaries
    // ========================================================

    static G4double GetEnergyLow(
        G4int bin
    );


    static G4double GetEnergyHigh(
        G4int bin
    );


    static G4double GetEnergyCenter(
        G4int bin
    );


private:

    TETModelImport* tetData;

    G4int HCID;

    G4THitsMap<G4double>* evtMap;

    static std::vector<G4double>
        energyEdges;


    static void InitializeEnergyBins();
};

#endif