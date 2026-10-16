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
// Internal.cc
// \file   MRCP_GEANT4/Internal/Internal.cc
// \author HUREL
//

#include <filesystem>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdexcept>


#include "TETDetectorConstruction.hh"
#include "TETModelImport.hh"
#include "TETPhysicsList.hh"
#include "TETActionInitialization.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "G4UIterminal.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "Randomize.hh"

void PrintUsage(){
	G4cerr << "Usage: ./Internal -i [SOURCE_ID[,SOURCE_ID...]] -m [MACRO] -o [OUTPUT] -f (option for MRCP-AF phantom)" << G4endl;
	G4cerr << "Example (one source ID): ./Internal -i 9500 -m run.mac -o run.out (-f)" << G4endl;
	G4cerr << "Example (multiple source IDs): ./Internal -i 9700,9900 -m run.mac -o run.out (-f)" << G4endl;
}

// Parse a comma-separated list of source-region IDs, e.g. "9700,9900".
std::vector<G4int> ParseSourceIDs(const std::string& sourceIDsText)
{
	std::vector<G4int> sourceIDs;
	std::stringstream stream(sourceIDsText);
	std::string token;

	while (std::getline(stream, token, ',')) {
		const auto first = token.find_first_not_of(" \t\r\n");
		const auto last  = token.find_last_not_of(" \t\r\n");

		if (first == std::string::npos) {
			throw std::invalid_argument("empty source ID in comma-separated list");
		}

		token = token.substr(first, last - first + 1);

		size_t parsedChars = 0;
		int sourceID = std::stoi(token, &parsedChars);
		if (parsedChars != token.size() || sourceID < 0) {
			throw std::invalid_argument("source IDs must be non-negative integers");
		}

		sourceIDs.push_back(static_cast<G4int>(sourceID));
	}

	if (sourceIDs.empty()) {
		throw std::invalid_argument("no source IDs were provided");
	}

	return sourceIDs;
}


int main(int argc,char** argv) 
{
// Read the arguments for batch mode
	//
	G4String macro;
	G4String output;
	std::vector<G4int> internalSources;
	G4bool   isAF(false);
	G4UIExecutive* ui = 0;

	for ( G4int i=1; i<argc; i++ ) {
		// macro file name
		if ( G4String(argv[i]) == "-m" ) {
			macro = argv[i+1];
			i++;
		}
		// output file name
		else if ( G4String(argv[i]) == "-o" ) {
			output = argv[i+1];
			i++;
		}
		// ID of the source organ
 		else if ( G4String(argv[i]) == "-i" ) {
			if (i + 1 >= argc) {
				G4cerr << "ERROR: -i requires one or more comma-separated source IDs." << G4endl;
				PrintUsage();
				return 1;
			}
			try {
				internalSources = ParseSourceIDs(argv[i + 1]);
			}
			catch (const std::exception& e) {
				G4cerr << "ERROR: Invalid source ID list '" << argv[i + 1]
				       << "': " << e.what() << G4endl;
				PrintUsage();
				return 1;
			}
			i++;
		}
		// switch for MRCP-AF phantom
		else if ( G4String(argv[i]) == "-f" ) {
			isAF = true;
		}
		else {
			PrintUsage();
		}
	}

	// Require at least one source-region ID.
	if (internalSources.empty()) {
		PrintUsage();
		return 1;
	}

	// print usage when there are more than eight arguments
	if ( argc>8 ){
		PrintUsage();
		return 1;
	}

	// Detect interactive mode (if no macro file name) and define UI session
	//
	if ( !macro.size() ) {
		ui = new G4UIExecutive(argc, argv, "csh");
	}
	// default output file name
	else if ( !output.size() ) output = macro + ".out";

	// Choose the Random engine
	//
	G4Random::setTheEngine(new CLHEP::RanecuEngine);
	G4Random::setTheSeed(time(0));

	// Construct the default run manager
	//
	#ifdef G4MULTITHREADED
		G4MTRunManager * runManager = new G4MTRunManager;
		// set the default number of threads as one
		runManager->SetNumberOfThreads(1);
	#else
		G4RunManager * runManager = new G4RunManager;
	#endif

	// Set a class to import phantom data
	//
	TETModelImport* tetData = new TETModelImport(isAF, ui);

	// Set mandatory initialisation classes
	//
	// detector construction
	runManager->SetUserInitialization(new TETDetectorConstruction(tetData));
	// physics list
	runManager->SetUserInitialization(new TETPhysicsList());
	// user action initialisation
	auto* actions =
 		new TETActionInitialization(
			tetData,
			internalSources,
			output
		);

	runManager->SetUserInitialization(actions);

	// Visualization manager
	//
	G4VisManager* visManager = new G4VisExecutive;
	visManager->Initialise();

	// Process macro or start UI session
	//
	G4UImanager* UImanager = G4UImanager::GetUIpointer();
/*
	if ( ! ui ){
		// batch mode
		G4String command = "/control/execute ";
		UImanager->ApplyCommand(command+macro);
	}*/

	if (!ui) {
		// batch mode
		namespace fs = std::filesystem;
		fs::path macroPath(macro.data());

		G4String macroDirectory = macroPath.parent_path().string();

		if (!macroDirectory.empty()) {
			UImanager->ApplyCommand("/control/macroPath " + macroDirectory);
		}

		// Execute the .in using its full path
		UImanager->ApplyCommand("/control/execute " + G4String(macroPath.string()));
	}

	else {
		// interactive mode
		UImanager->ApplyCommand("/control/execute init_vis.mac");
		ui->SessionStart();
		delete visManager;
		delete ui;
	}
namespace fs = std::filesystem;

	// Job termination
	delete runManager;
}


