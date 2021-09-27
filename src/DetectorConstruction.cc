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
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Tet.hh"
#include "Eigen/Core"

using namespace std;

DetectorConstruction::DetectorConstruction()
:worldLogical(0) ,worldPhysical(0)
{
	// Operating Table (w/ patient, curtain)
//	table_ocr = G4ThreeVector(-280,-960,-10); // Initial position
	// table_trans = G4ThreeVector(-80*mm, 20*mm -10*mm);     // Operating position
	// table_pivot_angle = 0 * deg;              // Z axis rotation
	// table_rotation_center = G4ThreeVector(1200*mm, 0*mm, 820*mm);

	// C-arm det
	carm_isocenter = G4ThreeVector(1120 * mm, 600 * mm,1135 * mm);
	// carm_primary   = 20 * deg;   // +LAO, -RAO
	// carm_secondary = 20 * deg;   // +CAU, -CRA

	messenger = new DetectorMessenger(this);
}

DetectorConstruction::~DetectorConstruction()
{
	delete messenger;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{  
	SetupWorldGeometry();
	ConstructOperatingTable();
	// ConstructPatient();
	ConstructPbGlass();
	ConstructCarmDet();
	return worldPhysical;
}

void DetectorConstruction::SetupWorldGeometry()
{
	// Define the world box (size: 20*20*20 m3)
    G4double worldHalfX = 10. * m;
    G4double worldHalfY = 10. * m;
    G4double worldHalfZ = 10. * m;

	G4VSolid* worldSolid = new G4Box("worldSolid", worldHalfX, worldHalfY, worldHalfZ);
	worldLogical = new G4LogicalVolume(worldSolid,G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"),"worldLogical");
	worldPhysical = new G4PVPlacement(0,G4ThreeVector(), worldLogical,"worldPhysical", 0, false,0,false);
	G4VisAttributes* va_World = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
	va_World->SetForceWireframe(true);
	worldLogical->SetVisAttributes(va_World);

	// Define floor
	// G4double floorZ = 0 * cm;
	// G4Material* concrete = G4NistManager::Instance()->FindOrBuildMaterial("G4_CONCRETE");
	// G4Box* floor = new G4Box("floor", worldHalfX, worldHalfY, (floorZ+worldHalfZ)*0.5);
	// G4LogicalVolume* lv_floor = new G4LogicalVolume(floor, concrete, "lv_floor");
	// new G4PVPlacement(0, G4ThreeVector(0,0,(-worldHalfZ)+(floorZ+worldHalfZ)*0.5), lv_floor, "pv_floor", worldLogical,0,0);

	// G4VisAttributes* vis_floor  = new G4VisAttributes(G4Colour(0.8,0.8,0.8,0.5));
	// vis_floor->SetForceAuxEdgeVisible();
	// lv_floor->SetVisAttributes(vis_floor);
}

void DetectorConstruction::ConstructOperatingTable()
{
	G4ThreeVector tableHalfSize(250*mm, 1600*mm, 1.43*mm*0.5);
	G4ThreeVector curtainHalfSize(0.25*mm, 250*mm, 200*mm);

	G4LogicalVolume* lv_phantomBox = ConstructPatient();

	// frame
	G4double halfZ = tableHalfSize.z()+curtainHalfSize.z()+((G4Box*) lv_phantomBox->GetSolid())->GetZHalfLength();
	G4double halfX = max(tableHalfSize.x()+curtainHalfSize.x(),((G4Box*) lv_phantomBox->GetSolid())->GetXHalfLength());
	G4Box* frame = new G4Box("sol_frame", halfX, tableHalfSize.y(), halfZ);
	G4LogicalVolume* lv_frame = new G4LogicalVolume(frame, G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"), "lv_frame");
	frame_rotation_matrix = new G4RotationMatrix();
	frame_rotation_matrix->setAxis(G4ThreeVector(0,0,-1));
 	pv_frame = 	new G4PVPlacement(frame_rotation_matrix, table_translation, lv_frame, "pv_frame", worldLogical, false, 0);

	// phantom box
	new G4PVPlacement(0, G4ThreeVector(0,frame->GetYHalfLength()-((G4Box*) lv_phantomBox->GetSolid())->GetYHalfLength()-10*cm,halfZ-((G4Box*) lv_phantomBox->GetSolid())->GetZHalfLength()),lv_phantomBox, "phantom box", lv_frame, false, 0);

	// Operating table
	G4Box* table = new G4Box("sol_table", tableHalfSize.x(), tableHalfSize.y(), tableHalfSize.z());
	G4LogicalVolume* lv_table = new G4LogicalVolume(table, G4NistManager::Instance()->FindOrBuildMaterial("G4_Al"), "lv_table");
	auto pv_table = new G4PVPlacement(0, G4ThreeVector(0,0,halfZ-((G4Box*) lv_phantomBox->GetSolid())->GetZHalfLength()*2.-tableHalfSize.z()), lv_table, "pv_table", lv_frame, false, 0);
	lv_table->SetVisAttributes( new G4VisAttributes(G4Colour(1.,1.,0.)) );
	frame_ralative_to_table = -pv_table->GetTranslation();

	// Pb Curtain - ? x ? cm lead apron, lead equivalence 0.5 mm Pb
	G4Box* curtain = new G4Box("sol_curtain", curtainHalfSize.x(), curtainHalfSize.y(), curtainHalfSize.z());
	G4ThreeVector curtain_margin(0*mm,10*cm,0*mm);
	G4ThreeVector relative_position_to_table
		= ( G4ThreeVector(-table->GetXHalfLength()-curtain->GetXHalfLength(), table->GetYHalfLength()-curtain->GetYHalfLength(), -curtain->GetZHalfLength())
						- curtain_margin);
	G4LogicalVolume* lv_curtain = new G4LogicalVolume(curtain, G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb"), "lv_curtain");
	new G4PVPlacement(0, relative_position_to_table + pv_table->GetTranslation(), lv_curtain, "pv_curtain", lv_frame, false, 0);
	lv_curtain->SetVisAttributes( new G4VisAttributes(G4Colour(0.,0.,1.,0.8)) );
}

G4LogicalVolume* DetectorConstruction::ConstructPatient()
{
	Eigen::MatrixXd V;
	ifstream ifs("./phantoms/patient.node");
	if(!ifs.is_open()){
		cout<<"there is no ./phantoms/patient.node!"<<endl;
		exit(100);
	}
	int rowV, tmp;
	ifs>>rowV>>tmp>>tmp>>tmp;
	V.resize(rowV, 3);
	for(int i=0;i<rowV;i++)
	{
		G4double x, y, z;
		ifs>>tmp>>x>>y>>z;
		V.row(i) = Eigen::RowVector3d(x,y,z)*cm;
	}
	ifs.close();

	//move the bbox center to the origin
	Eigen::RowVector3d max= V.colwise().maxCoeff();
	Eigen::RowVector3d min= V.colwise().minCoeff();
	V = V.rowwise()-(max+min)*0.5;

	//set the phantom box with 5-cm-margin
	Eigen::RowVector3d hlafSize = (max-min)*0.5;
	G4VSolid* phantomBox = new G4Box("patientBox", hlafSize.x(), hlafSize.y(), hlafSize.z());
	G4LogicalVolume* lv_phantomBox = new G4LogicalVolume(phantomBox, G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"), "patientBox");
	
	//construc tets
	ifstream ifsEle("./phantoms/patient.ele");
	if(!ifsEle.is_open()){
		cout<<"there is no ./phantoms/patient.ele!"<<endl;
		exit(100);
	}
	int numT;
	ifsEle>>numT>>tmp>>tmp;
	G4Material* bone = G4NistManager::Instance()->FindOrBuildMaterial("G4_BONE_CORTICAL_ICRP");
	G4Material* lung = G4NistManager::Instance()->FindOrBuildMaterial("G4_LUNG_ICRP");
	G4Material* tissue = G4NistManager::Instance()->FindOrBuildMaterial("G4_TISSUE_SOFT_ICRP");
	for(int i=0;i<numT;i++)
	{
		int a, b, c, d, id;
		ifsEle>>tmp>>a>>b>>c>>d>>id;
		G4VSolid* tet = new G4Tet("tet", G4ThreeVector(V(a,0),V(a,1),V(a,2)),
		                                 G4ThreeVector(V(b,0),V(b,1),V(b,2)),
										 G4ThreeVector(V(c,0),V(c,1),V(c,2)),
										 G4ThreeVector(V(d,0),V(d,1),V(d,2)));
		G4Material* mat;
		if(id<100) mat = bone;
		else if(id==125) mat = tissue;
		else if(id==158||id==159) mat=lung;
		else{
			cout<<"wrong ID: "<<id<<endl;
			exit(100);
		}
		G4LogicalVolume* lv_tet = new G4LogicalVolume(tet, mat, "tet");
		new G4PVPlacement(0, G4ThreeVector(), lv_tet, "tet", lv_phantomBox, false, 0);
		// lv_tet->SetVisAttributes(G4VisAttributes::GetInvisible());
	}
	ifsEle.close();
	return lv_phantomBox;
}

void DetectorConstruction::ConstructCarmDet()
{
	// This member function didn't considered the set functions of G4VPhyscialVolume.
	carm_rotation_matrix = new G4RotationMatrix;

	// C-arm	
	G4LogicalVolume* lv_det = new G4LogicalVolume(new G4Box("pv_det", 42*cm*0.5, 52*cm*0.5, 1*cm), G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb"), "lv_det");
	lv_det->SetVisAttributes(new G4VisAttributes(G4Colour(1.0, 1.0, 1.0,0.5)));
	pv_det = new G4PVPlacement(carm_rotation_matrix, G4ThreeVector(), lv_det, "pv_det", worldLogical, false, 0);
}

void DetectorConstruction::ConstructPbGlass()
{
	// Pb Glass - 40 x 50 cm tiltable lead acrylic shield, lead equivalence 0.5 mm Pb
	G4Box* glass = new G4Box("sol_glass", 200*mm, 250*mm, 0.25*mm);
	G4LogicalVolume* lv_glass = new G4LogicalVolume(glass, G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb"), "lv_glass");
	lv_glass->SetVisAttributes( new G4VisAttributes(G4Colour(0.,1.,0.,0.8)) );
	pv_glass  = new G4PVPlacement(new G4RotationMatrix, G4ThreeVector(), lv_glass, "pv_glass", worldLogical, false, 0);
}






