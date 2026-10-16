//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// TETPrimaryGeneratorAction.hh
// \file   MRCP_GEANT4/Internal/include/TETPrimaryGeneratorAction.hh
// \author HUREL
//

#ifndef TETPrimaryGeneratorAction_h
#define TETPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4RandomDirection.hh"
#include "G4ThreeVector.hh"
#include "G4Tet.hh"

#include <cstddef>
#include <map>
#include <vector>

#include "TETModelImport.hh"

class TETPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    TETPrimaryGeneratorAction(TETModelImport* tetData,
                              const std::vector<G4int>& internalSources,
                              const G4String& outputFilePath);
    virtual ~TETPrimaryGeneratorAction();

    virtual void GeneratePrimaries(G4Event* anEvent);

  private:
    struct TetrahedronData {
        G4ThreeVector vertices[4];
        G4double cumulativeVolume;
    };

    struct SourceRegion {
        G4int organID;
        G4double mass;
        std::vector<TetrahedronData> tetrahedra;
        G4double totalVolume;
        G4int allocatedParticles;
        G4int cumulativeParticleEnd;
        G4double fractionalRemainder;
    };

    void PrepareAllocationForCurrentRun();

    std::vector<SourceRegion> fSourceRegions;
    std::map<G4int, std::size_t> fSourceIndexByID;
    G4ParticleGun* fParticleGun;
    G4int fLastPreparedRunID;
    G4int fPreparedEventCount;
    G4String fOutputFilePath;
};

#endif

