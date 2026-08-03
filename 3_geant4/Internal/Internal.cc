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
#include <chrono>
#include <fstream>

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
	G4cerr<< "Usage: ./Internal -i [SOURCE ORGAN ID] -m [MACRO] -o [OUTPUT] -f (option for MRCP-AF phantom)"  <<G4endl;
	G4cerr<< "Example: ./Internal -i 9500 -m run.mac -o run.out (-f)" <<G4endl;
}

int main(int argc,char** argv) 
{
	auto totalStart = std::chrono::steady_clock::now();
	auto setupStart = std::chrono::steady_clock::now();

	double setupSeconds = 0.0;
	double runSeconds   = 0.0;

	// Read the arguments for batch mode
	//
	G4String macro;
	G4String output;
	G4int    internalSource(-1);
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
			internalSource = G4UIcommand::ConvertToInt(G4String(argv[i+1]));
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

	// print usage when there is no internal source input
	if ( internalSource<0 ) {
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
			internalSource,
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

	auto setupEnd = std::chrono::steady_clock::now();

	setupSeconds = std::chrono::duration<double>(setupEnd - setupStart).count();

	//if ( ! ui ){
		// batch mode
	//	G4String command = "/control/execute ";
	//	UImanager->ApplyCommand(command+macro);
	//}

	auto runStart = std::chrono::steady_clock::now();
	
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

		auto runEnd = std::chrono::steady_clock::now();

		runSeconds = std::chrono::duration<double>(runEnd - runStart).count();
	}

	else {
		// interactive mode
		UImanager->ApplyCommand("/control/execute init_vis.mac");
		ui->SessionStart();
		delete visManager;
		delete ui;
	}

	double transportSeconds =
    actions->GetRunAction()->GetTransportSeconds();

	// Job termination
	//
	delete runManager;

	auto totalEnd = std::chrono::steady_clock::now();

	double totalSeconds =
		std::chrono::duration<double>(totalEnd - totalStart).count();

	namespace fs = std::filesystem;

	fs::path timingPath(output.data());

	std::string filename = timingPath.filename().string();

	const std::string depositPrefix = "Geant4_deposit_";
	const std::string timingPrefix  = "Geant4_timing_";

	if (filename.rfind(depositPrefix, 0) == 0)
	{
		filename.replace(
			0,
			depositPrefix.size(),
			timingPrefix
		);
	}

	timingPath = timingPath.parent_path() / filename;
	timingPath.replace_extension(".txt");

	std::ofstream timing(timingPath);

	if (timing)
	{
		timing << "Geant4 Timing Summary\n";
		timing << "=====================\n\n";

		timing << "Macro      : " << macro << "\n";
		timing << "Output     : " << output << "\n\n";

		timing << "Setup time      : "
			<< setupSeconds
			<< " s\n";

		timing << "BeamOn time  : "
		<< transportSeconds
		<< " s\n";

		timing << "Execution time  : "
			<< runSeconds
			<< " s\n";

		timing << "Total wall time : "
			<< totalSeconds
			<< " s\n";
	}

	G4cout << "\nTotal wall time: "
		<< totalSeconds
		<< " s\n";
}


