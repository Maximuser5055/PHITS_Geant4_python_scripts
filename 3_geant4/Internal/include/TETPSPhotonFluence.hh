#ifndef TETPSPhotonFluence_h
#define TETPSPhotonFluence_h 1

#include "G4VPrimitiveScorer.hh"
#include "G4THitsMap.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"
#include "TETModelImport.hh"

#include <vector>

class TETPSPhotonFluence : public G4VPrimitiveScorer
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

    virtual void DrawAll() override {}

    virtual void PrintAll() override {}

    G4THitsMap<G4double>* GetHitsMap();

    // Energy bin information
    static const G4int nEnergyBins = 101;

    static G4double GetEnergy(G4int bin);

private:

    TETModelImport* tetData;

    G4int HCID;

    G4THitsMap<G4double>* evtMap;

    // Log-spaced energy bins:
    // 0.001 MeV -> 1.0 MeV
    static std::vector<G4double> energyBins;

    static void InitializeEnergyBins();
};

#endif