#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <Eigen/SparseCore>

Eigen::Vector3d ROOT;
std::vector<Eigen::Quaterniond,Eigen::aligned_allocator<Eigen::Quaterniond> > VQ;

// Geant4
#include "G4UImanager.hh"
#include "G4RunManagerFactory.hh"

#include "include/DetectorConstruction.hh"
#include "include/ParallelPhantom.hh"
#include "FTFP_BERT.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4ParallelWorldPhysics.hh"
#include "include/ActionInitialization.hh"

#include "G4Timer.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "Randomize.hh"

#include "include/TETModelImport.hh"
#include "include/PhantomAnimator.hh"


int main(int argc, char** argv)
{
	G4Timer* initTimer = new G4Timer;
	initTimer->Start();
	G4String macro;
	G4String output = "output";
	G4String phantomName;
	G4UIExecutive* ui = 0;
	G4bool isBack(false);
	G4bool isSQL(false);

	auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::MT);

	for ( G4int i=1; i<argc; i++ ) {
		// macro file name
		if ( G4String(argv[i]) == "-m" ) {
			macro = argv[++i];
		}
		// output file name
		else if ( G4String(argv[i]) == "-o" ) {
			output = argv[++i];
		}
		// phantom file name
		else if ( G4String(argv[i]) == "-p" ) {
			phantomName = argv[++i];
		}
		else if ( G4String(argv[i]) == "-v" )
		{
			ui = new G4UIExecutive(argc, argv, "Qt");
		}
		// SQL data flag are set to 0
		else if ( G4String(argv[i]) == "-r" ) {
			isBack = true;
			i++;
		}
		else if ( G4String(argv[i]) == "-np" )
		{
			runManager->SetNumberOfThreads(atoi(argv[i+1]));
			i++;
		}
		else if ( G4String(argv[i]) == "-sql" )
		{
			isSQL = true;
			i++;
		}
		else {
			cout << "argument check" << endl;
			return 1;
		}
	}

	// default output file name
	if ( !output.size() ) output = macro + ".out";

	// Choose the Random engine
	G4Random::setTheEngine(new CLHEP::RanecuEngine);
	G4Random::setTheSeed(time(0));

	// Import phantom
	TETModelImport* tetData = new TETModelImport(phantomName);

	// Set mandatory initialization classes
	auto det = new DetectorConstruction();
	det->RegisterParallelWorld(new ParallelPhantom("parallel", tetData));
	runManager->SetUserInitialization(det);
	G4VModularPhysicsList* physicsList = new FTFP_BERT;
	// physicsList->RegisterPhysics(new G4StepLimiterPhysics());
	physicsList->RegisterPhysics(new G4ParallelWorldPhysics("parallel", true));
	runManager->SetUserInitialization(physicsList);
	runManager->SetUserInitialization(new ActionInitialization(tetData, output, initTimer));

	// Initialize visualization
	G4VisManager* visManager = new G4VisExecutive("Quiet");
	visManager->Initialize();

	// Get the pointer to the User Interface manager  virtual
	G4UImanager* UImanager = G4UImanager::GetUIpointer();

	if (ui)
	{
		// interactive mode
		UImanager->ApplyCommand("/control/execute init_vis.mac");
		ui->SessionStart();
		delete ui;

		delete visManager;
		delete runManager;
		cout << "EXIT_SUCCESS" << endl;
		return EXIT_SUCCESS;
	}

	if(isSQL)
	{
		// batch mode
		runManager->Initialize();
		UImanager->ApplyCommand("/machine/c-arm 0 0 119.5");

		// Initialize DB system
		MYSQL* conn = mysql_init(NULL);
		MYSQL_RES* result;
		MYSQL_ROW row;

		string DB_HOST, DB_USER, DB_PASS, DB_NAME, DB_TABLE;
		int  PORT_ID;

		ifstream ifs("./pdc.dat");
		if(!ifs.is_open()) {
			cerr << "pdc.dat file is not opened" << endl; exit(1);
		}
		string dump;
		ifs >> dump >> DB_HOST;
		ifs >> dump >> DB_USER;
		ifs >> dump >> DB_PASS;
		ifs >> dump >> DB_NAME;
		ifs >> dump >> DB_TABLE;
		ifs >> dump >> PORT_ID;
		ifs.close();

		// Connect SQL server
		if( mysql_real_connect(conn, DB_HOST.c_str(), DB_USER.c_str(), DB_PASS.c_str(), DB_NAME.c_str(), PORT_ID, NULL, 0) != NULL ) {
			cout << "=======================================================================" << endl;
			cout << "SQL server is successfully connected" << endl;
			cout << "DB Host IP: " << DB_HOST << endl;
			cout << "DB User   : " << DB_USER << endl;
			cout << "DB Name   : " << DB_NAME << endl;
			cout << "DB_TABLE  : " << DB_TABLE << endl;
			cout << "My IP     : " << GetMyIPAddress() << endl;
			cout << "Port ID   : " << PORT_ID << endl;
			cout << "=======================================================================" << endl;
		}
		else
			finish_with_error(conn);


		// Main loop
		cout << "Start PDC Module !!" << endl;
		int frameNo(0); string sql;
		while(1)
		{
			if (isBack)
			{
				sql = "SELECT * FROM "+DB_TABLE+";";
				if(mysql_query(conn, sql.c_str()) == 1) {
					finish_with_error(conn);
				}
				result = mysql_store_result(conn);
				while ( (row = mysql_fetch_row(result)) != NULL )
				{
					frameNo = atoi(row[0]);
					cout << "Update Frame " << frameNo << " status => 0" << endl;
					sql = "UPDATE "+DB_TABLE+" SET Flag = 0 WHERE FrameNo = " + to_string(frameNo) + ";";
					if(mysql_query(conn, sql.c_str()) == 1) {
						finish_with_error(conn);
					}
				}
				break;
			}
			// Check the Frame No
			sql = "SELECT FrameNo, Flag FROM "+DB_TABLE+" WHERE FrameNo = -1;";
			if(mysql_query(conn, sql.c_str()) == 1) {
				finish_with_error(conn);
			}
			result = mysql_store_result(conn);
			row = mysql_fetch_row(result);
			frameNo = atoi(row[1]);

			if (frameNo == 410) {
				cout << "END" << endl;
				break;
			}

			// Set FrameNo ==> Show Progress
			sql = "UPDATE "+DB_TABLE+" SET Flag = " + to_string(frameNo + 1) + " WHERE FrameNo = -1;";
			if(mysql_query(conn, sql.c_str()) == 1) {
				finish_with_error(conn);
			}

			// Select Frame ----------------------------
			//
			sql = "SELECT * FROM "+DB_TABLE+" WHERE FrameNo = " + to_string(frameNo) + ";";
			if(mysql_query(conn, sql.c_str()) == 1) {
				finish_with_error(conn);
			}
			result = mysql_store_result(conn);
			row = mysql_fetch_row(result);

			int calcFlag = atoi(row[2]);
			if (calcFlag != 0) {
				cerr << "WARNING!! The frame is already allocated in other server." << endl;
				continue;
			}

			// Radiologist
			ROOT = Vector3d( atof(row[4]),atof(row[5]),atof(row[6]) );
			for (int i=0;i<22;i++) {
				double w = atof(row[4*i+0 + 7]);
				double x = atof(row[4*i+1 + 7]);
				double y = atof(row[4*i+2 + 7]);
				double z = atof(row[4*i+3 + 7]);
				VQ.push_back(Quaterniond(w,x,y,z));
			}

			// Operating Table

			// C-arm

			// Pb Glass


			// ---------------------------------------



			// flag On
			cout << "FrameNo (status): " << frameNo << " (" << calcFlag << ") => Update status (1)" << endl;
			sql = "UPDATE "+DB_TABLE+" SET Flag = 1 WHERE FrameNo = " + to_string(frameNo) + ";";
			if(mysql_query(conn, sql.c_str()) == 1) {
				finish_with_error(conn);
			}
			sql = "UPDATE "+DB_TABLE+" SET IP = '"+GetMyIPAddress()+"' WHERE FrameNo = " + to_string(frameNo) + ";";
			if(mysql_query(conn, sql.c_str()) == 1) {
				finish_with_error(conn);
			}
			UImanager->ApplyCommand("/phantom/frameSQL "+to_string(frameNo));
			VQ.clear();

			runManager->BeamOn(100000);
			cout << "Update status (2)" << endl;

			sql = "UPDATE "+DB_TABLE+" SET Flag = 2 WHERE FrameNo = " + to_string(frameNo) + ";";
			if(mysql_query(conn, sql.c_str()) == 1) {
				finish_with_error(conn);
			}
		}
		mysql_close(conn);
	}
	else
	{
		if ( ! ui ){
			// batch mode
			G4String command = "/control/execute ";
			UImanager->ApplyCommand(command+macro);
		}
		else {
			// interactive mode
			UImanager->ApplyCommand("/control/execute init_vis.mac");
			// runManager->Initialize();
			ui->SessionStart();
			delete ui;
		}
	}

	delete visManager;
	delete runManager;

    cout << "EXIT_SUCCESS" << endl;
    return EXIT_SUCCESS;
}
