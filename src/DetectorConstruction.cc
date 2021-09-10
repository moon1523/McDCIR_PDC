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
//

#include "DetectorConstruction.hh"

using namespace std;

DetectorConstruction::DetectorConstruction(TETModelImport* _tetData)
:worldLogical(0) ,worldPhysical(0), container_logic(0), container_phy(0),
 tetData(_tetData), lv_pic(0), pv_pic(0)
{
	// Material
	vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
	water = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
	lead = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb");

	// Radiologist
	G4int frameNo = 130;
	tetData->Deform(frameNo);
	doctor_translation = tetData->GetRootPosition(frameNo) * cm;
	phantomSize        = tetData -> GetPhantomSize();
	phantomBoxMin      = tetData -> GetPhantomBoxMin();
	phantomBoxMax      = tetData -> GetPhantomBoxMax();
	nOfTetrahedrons    = tetData -> GetNumTetrahedron();


	// Operating Table (Patient + Glass)
//	table_ocr = G4ThreeVector(-280,-960,-10); // Initial position
	table_ocr = G4ThreeVector(-80,20-10);     // Operating position
	table_pivot_angle = 0 * deg;              // Z axis rotation
	table_rotation_center = G4ThreeVector(1200, 0, 820);
	table_rotation_matrix = new G4RotationMatrix;
	table_rotation_matrix->rotateZ(table_pivot_angle);
	table_rotation_matrix2 = new G4RotationMatrix;
	table_rotation_matrix2->rotateZ(-table_pivot_angle); // to rotate solid
	table_translation = (*table_rotation_matrix) * table_ocr + table_rotation_center;

	// Pb glass
	glass_translation = G4ThreeVector(950, 500, 1500);
	AngleAxisd aa(80*deg, Vector3d::UnitX());
	Quaterniond glass_quaternion(aa);
	MatrixXd mat = glass_quaternion.matrix();
	G4ThreeVector col1(mat(0,0), mat(1,0), mat(2,0));
	G4ThreeVector col2(mat(0,1), mat(1,1), mat(2,1));
	G4ThreeVector col3(mat(0,2), mat(1,2), mat(2,2));
	glass_rotation_matrix = new G4RotationMatrix(-col1,-col2,-col3);


	// C-arm (visualization)
	carm_isocenter = G4ThreeVector(1120,600,1135);
	carm_primary   = 20 * deg;   // +LAO, -RAO
	carm_secondary = 20 * deg;   // +CAU, -CRA
}

DetectorConstruction::~DetectorConstruction()
{
	delete tetData;
	delete table_rotation_matrix;
	delete table_rotation_matrix2;
	delete glass_rotation_matrix;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{  
	SetupWorldGeometry();
	ConstructPhantom();
	PrintPhantomInformation();

	ConstructOperatingTable();
	ConstructPatient();
	ConstructPbGlass();
	ConstructCarm();

	return worldPhysical;
}

void DetectorConstruction::SetupWorldGeometry()
{
	// Define the world box (size: 10*10*10 m3)
	//
    G4double worldHalfX = 10. * m;
    G4double worldHalfY = 10. * m;
    G4double worldHalfZ = 10. * m;

	G4VSolid* worldSolid = new G4Box("worldSolid", worldHalfX, worldHalfY, worldHalfZ);
	worldLogical = new G4LogicalVolume(worldSolid,vacuum,"worldLogical");
	worldPhysical = new G4PVPlacement(0,G4ThreeVector(), worldLogical,"worldPhysical", 0, false,0,false);
	G4VisAttributes* va_World = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
	va_World->SetForceWireframe(true);
	worldLogical->SetVisAttributes(va_World);

	// Define floor
	G4double floorZ = 0 * cm;
	G4Material* concrete = G4NistManager::Instance()->FindOrBuildMaterial("G4_CONCRETE");
	G4Box* floor = new G4Box("floor", worldHalfX, worldHalfY, (floorZ+worldHalfZ)*0.5);
	G4LogicalVolume* lv_floor = new G4LogicalVolume(floor, concrete, "lv_floor");
	new G4PVPlacement(0, G4ThreeVector(0,0,(-worldHalfZ)+(floorZ+worldHalfZ)*0.5), lv_floor, "pv_floor", worldLogical,0,0);

	G4VisAttributes* vis_floor  = new G4VisAttributes(G4Colour(0.8,0.8,0.8,0.5));
	vis_floor->SetForceAuxEdgeVisible();
	lv_floor->SetVisAttributes(vis_floor);

	// Define the phantom container (10-cm margins from the bounding box of phantom)
	//
	G4ThreeVector max = phantomBoxMax;
	G4ThreeVector min = phantomBoxMin;
	G4double dimX = max.getX()>-min.getX()? max.getX():-min.getX();
	G4double dimY = max.getY()>-min.getY()? max.getY():-min.getY();
	G4double dimZ = max.getZ()>-min.getZ()? max.getZ():-min.getZ();

	doctor_translation.setZ(doctor_translation.getZ() + (dimZ - doctor_translation.getZ()) + floorZ);

	G4Box* containerSolid = new G4Box("phantomBox", dimX + 1.*cm,
										            dimY + 1.*cm,
										            dimZ + 1.*cm);

	container_logic = new G4LogicalVolume(containerSolid, vacuum, "phantomLogical");
	container_phy = new G4PVPlacement(0, doctor_translation, container_logic, "PhantomPhysical", worldLogical, false, 0);

	G4VisAttributes* container_vis = new G4VisAttributes(G4Colour(1.0,1.0,1.0,0.0));
	container_vis->SetForceWireframe(false);
	container_logic->SetVisAttributes(container_vis);
	container_logic->SetOptimisation(TRUE);
	container_logic->SetSmartless( 0.5 ); // for optimization (default=2)
}

void DetectorConstruction::ConstructPhantom()
{
    lv_pic = new G4LogicalVolume(tetData->GetPicTess(), water, "lv_pic");
	pv_pic = new G4PVPlacement(0, G4ThreeVector(), lv_pic, "pv_pic", container_logic, false, 0);
	lv_pic->SetVisAttributes(new G4VisAttributes(G4Colour(1.000000, 0.752941, 0.627451)));
}

void DetectorConstruction::ConstructSDandField()
{

}

void DetectorConstruction::PrintPhantomInformation()
{
	// print brief information on the imported phantom
	G4cout<< G4endl;
	G4cout.precision(3);
	G4cout<<"   Phantom name               "<<tetData->GetPhantomName() << " TET phantom"<<G4endl;
	G4cout<<"   Phantom size               "<<phantomSize.x()<<" * "<<phantomSize.y()<<" * "<<phantomSize.z()<<" mm3"<<G4endl;
	G4cout<<"   Phantom box position (min) "<<phantomBoxMin.x()<<" mm, "<<phantomBoxMin.y()<<" mm, "<<phantomBoxMin.z()<<" mm"<<G4endl;
	G4cout<<"   Phantom box position (max) "<<phantomBoxMax.x()<<" mm, "<<phantomBoxMax.y()<<" mm, "<<phantomBoxMax.z()<<" mm"<<G4endl;
	G4cout<<"   Number of tetrahedrons     "<<nOfTetrahedrons<<G4endl<<G4endl;
}


void DetectorConstruction::ConstructOperatingTable()
{
	// Operating table
	G4Box* table = new G4Box("sol_table", 250, 1600, 65);
	G4LogicalVolume* lv_table = new G4LogicalVolume(table, vacuum, "lv_table");
	new G4PVPlacement(table_rotation_matrix2, table_translation, lv_table, "pv_table", worldLogical, false, 0);
	lv_table->SetVisAttributes( new G4VisAttributes(G4Colour(1.,1.,0.)) );

	// Pb Curtain - 35 x 50 cm lead apron, lead equivalence 0.5 mm Pb
	G4ThreeVector relative_position_to_table = (*table_rotation_matrix) * G4ThreeVector(-250-10, 1600-300-500, -300);
	G4ThreeVector curtain_translation = table_translation + relative_position_to_table;
	G4Box* curtain = new G4Box("sol_curtain", 0.25, 175, 250);
	G4LogicalVolume* lv_curtain = new G4LogicalVolume(curtain, lead, "lv_curtain");
	new G4PVPlacement(table_rotation_matrix2, curtain_translation, lv_curtain, "pv_curtain", worldLogical, false, 0);
	lv_curtain->SetVisAttributes( new G4VisAttributes(G4Colour(0.,0.,1.,0.8)) );


	//	G4Sphere* test = new G4Sphere("solid", 0, 20,0*deg,360*deg,0*deg,360*deg);
	//	G4LogicalVolume* lv_test = new G4LogicalVolume(test, vacuum, "lv_table");
	//	G4VPhysicalVolume* pv_test = new G4PVPlacement(0, table_rotation_center, lv_test, "pv_table", worldLogical, false, 0);
	//	lv_test->SetVisAttributes( new G4VisAttributes(G4Colour(1.,0.,0.)) );
	//
	//	G4Sphere* test2 = new G4Sphere("solid", 0, 20,0*deg,360*deg,0*deg,360*deg);
	//	G4LogicalVolume* lv_test2 = new G4LogicalVolume(test2, vacuum, "lv_table");
	//	G4VPhysicalVolume* pv_test2 = new G4PVPlacement(0, table_rotation_center + table_ocr, lv_test2, "pv_table", worldLogical, false, 0);
	//	lv_test2->SetVisAttributes( new G4VisAttributes(G4Colour(0.,1.,0.)) );
	//
	//	G4Sphere* test4 = new G4Sphere("solid", 0, 20,0*deg,360*deg,0*deg,360*deg);
	//	G4LogicalVolume* lv_test4 = new G4LogicalVolume(test4, vacuum, "lv_table");
	//	G4VPhysicalVolume* pv_test4 = new G4PVPlacement(0, table_translation, lv_test4, "pv_table", worldLogical, false, 0);
	//	lv_test4->SetVisAttributes( new G4VisAttributes(G4Colour(1.,0.,1.)) );
	//
	//
	//	G4Sphere* test3 = new G4Sphere("solid", 0, 50,0*deg,360*deg,0*deg,360*deg);
	//	G4LogicalVolume* lv_test3 = new G4LogicalVolume(test3, vacuum, "lv_table");
	//	G4VPhysicalVolume* pv_test3 = new G4PVPlacement(0, (*table_rotation_matrix) * table_ocr, lv_test3, "pv_table", worldLogical, false, 0);
	//	lv_test3->SetVisAttributes( new G4VisAttributes(G4Colour(0.,0.,1.)) );
	//
	//	G4Sphere* test5 = new G4Sphere("solid", 0, 50,0*deg,360*deg,0*deg,360*deg);
	//	G4LogicalVolume* lv_test5 = new G4LogicalVolume(test5, vacuum, "lv_table");
	//	G4VPhysicalVolume* pv_test5 = new G4PVPlacement(0, table_ocr, lv_test5, "pv_table", worldLogical, false, 0);
	//	lv_test5->SetVisAttributes( new G4VisAttributes(G4Colour(0.,0.,1.)) );
}

void DetectorConstruction::ConstructPatient()
{
	ifstream ifs("../phantoms/patient.ply");
	if (!ifs.is_open()) {
		cerr << "patient.ply is not opened" << endl;
		return;
	}
	G4TessellatedSolid* pTess = new G4TessellatedSolid();

	G4String dump;

	G4int vertNum, faceNum;
	G4ThreeVector vertex, face;
	vector<G4ThreeVector> vertVec;

	G4double minZ( DBL_MAX);
	while (getline(ifs,dump))
	{
		stringstream ss(dump);
		ss >> dump;
		if (dump == "element") {
			ss >> dump;
			if (dump == "vertex") {
				ss >> vertNum;
			}
			if (dump == "face") {
				ss >> faceNum;
			}
		}

		if (dump == "end_header") {
			for (G4int i=0; i<vertNum; i++) {
				G4double x,y,z;
				getline(ifs,dump);
				stringstream sv(dump);
				sv >> x >> y >> z;
				vertVec.push_back(G4ThreeVector(x*cm,y*cm,z*cm));

				if (z*cm < minZ) minZ = z*cm;
			}
			for (G4int i=0; i<faceNum; i++) {
				G4int a,b,c;
				getline(ifs,dump);
				stringstream sf(dump);
				sf >> dump >> a >> b >> c;
				pTess->AddFacet(new G4TriangularFacet(vertVec[a], vertVec[b], vertVec[c], ABSOLUTE));

			}
		}
	}
	pTess->SetSolidClosed(true);

	// Patient
	G4ThreeVector relative_position_to_table = (*table_rotation_matrix) * G4ThreeVector(0,300,65-minZ+1*cm);
	G4ThreeVector patient_translation = table_translation + relative_position_to_table;
	G4LogicalVolume* lv_patient = new G4LogicalVolume(pTess, water, "lv_patient");
	new G4PVPlacement(table_rotation_matrix2, patient_translation, lv_patient, "pv_patient", worldLogical, false, 0);
	lv_patient->SetVisAttributes(new G4VisAttributes(G4Colour(1.0, 1.0, 1.0)));
}

void DetectorConstruction::ConstructCarm()
{
	// This member function didn't considered the set functions of G4VPhyscialVolume.
	G4RotationMatrix* carm_rotation_matrix = new G4RotationMatrix;
	carm_rotation_matrix->rotateX(carm_secondary).rotateY(carm_primary);

	ifstream ifs("../phantoms/carm.ply");
	if (!ifs.is_open()) {
		cerr << "carm.ply is not opened" << endl;
		exit(1);
	}
	G4TessellatedSolid* tess_carm = new G4TessellatedSolid();

	G4String dump;
	G4int vertNum, faceNum;
	G4ThreeVector vertex, face;
	vector<G4ThreeVector> vertVec;

	while (getline(ifs,dump))
	{
		stringstream ss(dump);
		ss >> dump;
		if (dump == "element") {
			ss >> dump;
			if (dump == "vertex") {
				ss >> vertNum;
			}
			if (dump == "face") {
				ss >> faceNum;
			}
		}

		if (dump == "end_header") {
			for (G4int i=0; i<vertNum; i++) {
				G4double x,y,z;
				getline(ifs,dump);
				stringstream sv(dump);
				sv >> x >> y >> z;
				vertVec.push_back((*carm_rotation_matrix)*(G4ThreeVector(x,y,z)));
			}
			for (G4int i=0; i<faceNum; i++) {
				G4int a,b,c;
				getline(ifs,dump);
				stringstream sf(dump);
				sf >> dump >> a >> b >> c;
				tess_carm->AddFacet(new G4TriangularFacet(vertVec[a], vertVec[b], vertVec[c], ABSOLUTE));
			}
		}
	}
	tess_carm->SetSolidClosed(true);

	// C-arm
	G4LogicalVolume* lv_carm = new G4LogicalVolume(tess_carm, vacuum, "lv_carm");
	new G4PVPlacement(0, carm_isocenter, lv_carm, "pv_carm", worldLogical, false, 0);
	lv_carm->SetVisAttributes(new G4VisAttributes(G4Colour(1.0, 1.0, 1.0,0.5)));

	delete carm_rotation_matrix;
}

void DetectorConstruction::ConstructPbGlass()
{
	// Pb Glass - 40 x 50 cm tiltable lead acrylic shield, lead equivalence 0.5 mm Pb
	G4Box* glass = new G4Box("sol_glass", 200, 250, 0.25);
	G4LogicalVolume* lv_glass = new G4LogicalVolume(glass, lead, "lv_glass");
	new G4PVPlacement(glass_rotation_matrix, glass_translation, lv_glass, "pv_glass", worldLogical, false, 0);
	lv_glass->SetVisAttributes( new G4VisAttributes(G4Colour(0.,1.,0.,0.8)) );
}






