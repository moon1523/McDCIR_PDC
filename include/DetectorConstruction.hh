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

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"

#include <cmath>

#include "globals.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tet.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4Sphere.hh"

#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4GeometryManager.hh"

#include "TETModelImport.hh"
#include "G4PSEnergyDeposit.hh"

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction(TETModelImport* tetData);
    virtual ~DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

    // Operating Table
    void SetTableTrans(G4ThreeVector _table_trans) { table_trans = _table_trans; }
    void SetTablePivotAngle(G4double _table_pivot_angle) { table_pivot_angle = _table_pivot_angle; }
    // Glass
    void SetGlassPose(G4ThreeVector _glass_trans, Quaterniond _glass_quat) 
	{
		glass_trans = _glass_trans;
    	if (glass_rotation_matrix) delete glass_rotation_matrix;
    	MatrixXd mat = _glass_quaternion.matrix();
    	G4ThreeVector col1(mat(0,0), mat(1,0), mat(2,0));
    	G4ThreeVector col2(mat(0,1), mat(1,1), mat(2,1));
    	G4ThreeVector col3(mat(0,2), mat(1,2), mat(2,2));
    	glass_rotation_matrix = new G4RotationMatrix(-col1,-col2,-col3);
    }
    // C-arm will be updated


private:
    void SetupWorldGeometry();
	void ConstructPhantom();
	void PrintPhantomInformation();

	void ConstructOperatingTable();
	void ConstructPatient();
	void ConstructCarm();
	void ConstructPbGlass();

	G4LogicalVolume*   worldLogical;
	G4VPhysicalVolume* worldPhysical;
	G4LogicalVolume*   container_logic;
	G4VPhysicalVolume* container_phy;

	TETModelImport*    tetData;

	// Material
	G4Material* vacuum;
	G4Material* water;
	G4Material* lead;
	G4Material* carbonfiber;

	// Radiologist
	G4ThreeVector doctor_translation;
	G4LogicalVolume* lv_doctor;
	G4VPhysicalVolume* pv_doctor;

	// Operating Table
	G4ThreeVector table_trans;
	G4double table_pivot_angle;
	G4ThreeVector table_rotation_center;
	G4RotationMatrix* table_rotation_matrix;
	G4RotationMatrix* table_rotation_matrix2;
	G4ThreeVector table_translation;

	// Glass
	G4ThreeVector glass_trans;
	G4RotationMatrix* glass_rot;
	G4VPhysicalVolume* pv_glass;

	// C-arm
	G4ThreeVector carm_isocenter;
	G4double carm_primary;
	G4double carm_secondary;
};


#endif

