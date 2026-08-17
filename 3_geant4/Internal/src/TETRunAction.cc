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
// TETRunAction.cc
// \file   MRCP_GEANT4/Internal/src/TETRunAction.cc
// \author HUREL
//

#include "TETRunAction.hh"
#include "TETPSPhotonFluence.hh"

#include <filesystem>

TETRunAction::TETRunAction(TETModelImport* _tetData, G4String _output)
:tetData(_tetData), fRun(0), numOfEvent(0), runID(0), outputFile(_output)
{}

TETRunAction::~TETRunAction()
{}

G4Run* TETRunAction::GenerateRun()
{
	// generate run
	fRun = new TETRun();
	return fRun;
}

void TETRunAction::BeginOfRunAction(const G4Run* aRun)
{
	// print the progress at the interval of 10%
	numOfEvent=aRun->GetNumberOfEventToBeProcessed();
	G4RunManager::GetRunManager()->SetPrintProgress(int(numOfEvent*0.1));
}

void TETRunAction::EndOfRunAction(const G4Run* aRun)
{
	// print the result only in the Master
	if(!isMaster) return;

	// get the run ID
	runID = aRun->GetRunID();

	// print the run result by G4cout and std::ofstream
	//
	// print by G4cout
	PrintResult(G4cout);

	// print by std::ofstream
	std::ofstream ofs(outputFile.c_str());
	PrintResult(ofs);
	ofs.close();

	// Write a csv file
	std::filesystem::path csvPath(outputFile.c_str());

	csvPath.replace_extension(".csv");

	std::ofstream csv(csvPath);

	if (!csv)
	{
		G4cerr << "Unable to create " << csvPath << G4endl;
	}
	else
	{
		PrintCSV(csv);
		csv.close();
	}

    if (TETPSPhotonFluence::ENABLE_PHOTON_FLUENCE) 
    {
        std::filesystem::path fluxPath(outputFile.c_str());

        fluxPath.replace_extension(".csv");

        fluxPath = fluxPath.parent_path() / (fluxPath.stem().string() + "_photon_fluence.csv");

        std::ofstream fluxFile(fluxPath);
        
        if (!fluxFile) 
        {
		    G4cerr
			<< "Unable to create "
			<< fluxPath
			<< G4endl;
	    }
	    else
	    {
		    PrintPhotonFluence(fluxFile);
		    fluxFile.close();
	    }
    }
}

void TETRunAction::PrintResult(std::ostream &out)
{
	// Print run result
	//
	using namespace std;
	EDEPMAP edepMap = *fRun->GetEdepMap();

	out << G4endl
	    << "=====================================================================" << G4endl
	    << " Run #" << runID << " / Number of event processed : "<< numOfEvent     << G4endl
	    << "=====================================================================" << G4endl
	    << "organ ID| "
		<< setw(19) << "Organ Mass (g)"
		<< setw(19) << "Dose (Gy/source)"
		<< setw(19) << "Relative Error" << G4endl;

	out.precision(3);
	for(auto itr : tetData->GetMassMap()){
		G4double meanDose    = edepMap[itr.first].first  / itr.second / numOfEvent;
		G4double squareDose = edepMap[itr.first].second / (itr.second*itr.second);
		G4double variance    = ((squareDose/numOfEvent) - (meanDose*meanDose))/numOfEvent;
		G4double relativeE   = sqrt(variance)/meanDose;

		out << setw(8)  << itr.first << "| "
			<< setw(19) << fixed      << itr.second/g;
		out	<< setw(19) << scientific << meanDose/(joule/kg);
		out	<< setw(19) << fixed      << relativeE << G4endl;
	}
	out << "=====================================================================" << G4endl << G4endl;
}

void TETRunAction::PrintCSV(std::ostream& out)
{
    EDEPMAP edepMap = *fRun->GetEdepMap();

    out << "Organ ID,"
        << "Organ Mass (g),"
        << "Dose (Gy/source),"
        << "Relative Error\n";

    out.precision(10);

    for (auto itr : tetData->GetMassMap())
    {
        G4double meanDose =
            edepMap[itr.first].first /
            itr.second /
            numOfEvent;

        G4double squareDose =
            edepMap[itr.first].second /
            (itr.second * itr.second);

        G4double variance =
            ((squareDose / numOfEvent) -
             (meanDose * meanDose)) /
            numOfEvent;

        G4double relativeE =
            (meanDose > 0.0)
                ? std::sqrt(variance) / meanDose
                : 0.0;

        out << itr.first << ","
            << itr.second / g << ","
            << meanDose / (joule / kg) << ","
            << relativeE << "\n";
    }
}

void TETRunAction::PrintPhotonFluence(std::ostream& out)
{
    auto fluxMap =
        *fRun->GetPhotonFluenceMap();

    // ========================================================
    // Skeletal IDs
    // ========================================================

    std::vector<G4int> skeletalIDs =
    {
        1400,
        1500,
        1700,
        1800,
        2000,
        2100,
        2300,
        2500,
        2700,
        2900,
        3000,
        3200,
        3300,
        3500,
        3600,
        3800,
        4000,
        4200,
        4400,
        4600,
        4800,
        5000,
        5200,
        5400,
        5600
    };


    // ========================================================
    // CSV header
    // ========================================================

    out
        << "Skeletal ID,"
        << "Energy Low (MeV),"
        << "Energy High (MeV),"
        << "Energy Center (MeV),"
        << "Fluence (photons/m2/source)"
        << "\n";


    out.precision(12);


    // ========================================================
    // Loop over skeletal regions
    // ========================================================

    for (G4int organID : skeletalIDs)
    {
        // ----------------------------------------------------
        // Total volume of this skeletal region
        // ----------------------------------------------------

        G4double volume =
            tetData->GetVolume(
                organID
            );


        // ----------------------------------------------------
        // Loop over energy bins
        // ----------------------------------------------------

        for (
            G4int bin = 0;
            bin < TETPSPhotonFluence::nEnergyBins;
            bin++
        )
        {
            // ------------------------------------------------
            // Construct key
            // ------------------------------------------------

            G4int key =
                organID
                * TETPSPhotonFluence::nEnergyBins
                + bin;


            // ------------------------------------------------
            // Find track length
            // ------------------------------------------------

            auto itr =
                fluxMap.find(key);


            G4double totalTrackLength = 0.0;


            if (itr != fluxMap.end())
            {
                totalTrackLength =
                    itr->second.first;
            }


            // ------------------------------------------------
            // Mean track length per source particle
            // ------------------------------------------------

            G4double meanTrackLength =
                totalTrackLength
                / numOfEvent;


            // ------------------------------------------------
            // Track-length fluence estimator
            //
            // Phi = L / V
            // ------------------------------------------------

            G4double fluence =
                meanTrackLength
                / volume;


            // ------------------------------------------------
            // Convert to photons/m2
            // ------------------------------------------------

            G4double fluence_m2 =
                fluence * m2;


            // ------------------------------------------------
            // Energy-bin information
            // ------------------------------------------------

            G4double energyLow =
                TETPSPhotonFluence::GetEnergyLow(
                    bin
                );


            G4double energyHigh =
                TETPSPhotonFluence::GetEnergyHigh(
                    bin
                );


            G4double energyCenter =
                TETPSPhotonFluence::GetEnergyCenter(
                    bin
                );


            // ------------------------------------------------
            // Write row
            // ------------------------------------------------

            out
                << organID
                << ","
                << energyLow
                << ","
                << energyHigh
                << ","
                << energyCenter
                << ","
                << fluence_m2
                << "\n";
        }
    }
}