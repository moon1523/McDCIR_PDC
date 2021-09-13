// Geant4
#include "G4UImanager.hh"
#include "G4MTRunManager.hh"

#include "DetectorConstruction.hh"
#include "FTFP_BERT.hh"
#include "G4StepLimiterPhysics.hh"
#include "ActionInitialization.hh"

#include "G4Timer.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "Randomize.hh"

//
#include "TETModelImport.hh"
#include "PhantomAnimator.hh"

int main(int argc, char** argv)
{
	G4Timer* initTimer = new G4Timer;
	initTimer->Start();
	G4String macro;
	G4String output;
	G4String phantomName;
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
		// phantom file name
		else if ( G4String(argv[i]) == "-p" ) {
			phantomName = argv[i+1];
			i++;
		}
		else if ( G4String(argv[i]) == "-v" )
		{
			ui = new G4UIExecutive(argc, argv, "Qt");
		}
		else {
			cout << "argument check" << endl;
			return 1;
		}
	}

	// default output file name
	if ( !output.size() ) output = macro + ".out";

	G4MTRunManager* runManager = new G4MTRunManager();

	// Choose the Random engine
	G4Random::setTheEngine(new CLHEP::RanecuEngine);
	G4Random::setTheSeed(time(0));

	// Import phantom
	TETModelImport* tetData = new TETModelImport(phantomName);

	// Set mandatory initialization classes
	runManager->SetUserInitialization(new DetectorConstruction(tetData));
	G4VModularPhysicsList* physicsList = new FTFP_BERT;
	physicsList->RegisterPhysics(new G4StepLimiterPhysics());
	runManager->SetUserInitialization(physicsList);
	runManager->SetUserInitialization(new ActionInitialization(tetData, output, initTimer));

	// Initialize visualization
	G4VisManager* visManager = new G4VisExecutive("Quiet");
	visManager->Initialize();

	// Get the pointer to the User Interface manager
	G4UImanager* UImanager = G4UImanager::GetUIpointer();
	if ( ! ui ){
		// batch mode
		G4String command = "/control/execute ";
		UImanager->ApplyCommand(command+macro);
	}
	else {
		// interactive mode
		UImanager->ApplyCommand("/control/execute ../init_vis.mac");
		ui->SessionStart();
		delete ui;
	}

	delete visManager;
	delete runManager;

    cout << "EXIT_SUCCESS" << endl;
    return EXIT_SUCCESS;
}
