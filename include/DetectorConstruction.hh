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
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "DetectorMessenger.hh"

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction();
    virtual ~DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();

    // Operating Table
    void SetTablePose(G4ThreeVector _table_trans, G4double _table_pivot_angle) {  ;}
    // Glass
    void SetGlassPose(G4ThreeVector _glass_trans, G4ThreeVector _glass_axis, G4double _glass_theta) 
	{
		pv_glass->SetTranslation(_glass_trans);
		pv_glass->GetRotation()->setAxis(_glass_axis);
		pv_glass->GetRotation()->setTheta(-_glass_theta);
    }
    // C-arm will be updated


private:
    void SetupWorldGeometry();

	void ConstructOperatingTable();
	G4LogicalVolume* ConstructPatient();
	void ConstructCarmDet();
	void ConstructPbGlass();

	G4LogicalVolume*   worldLogical;
	G4VPhysicalVolume* worldPhysical;
	G4LogicalVolume*   container_logic;
	G4VPhysicalVolume* container_phy;

	// Material
	G4Material* vacuum;
	G4Material* water;
	G4Material* lead;
	G4Material* carbonfiber;

	// Operating Table
	G4VPhysicalVolume* pv_frame;
	G4ThreeVector table_rotation_center;
	G4RotationMatrix* table_rotation_matrix;
	G4RotationMatrix* table_rotation_matrix2;
	G4ThreeVector table_translation;

	// Glass
	G4VPhysicalVolume* pv_glass;

	// C-arm
	G4ThreeVector carm_isocenter;
	G4double carm_primary; //rao, lao
	G4double carm_secondary; //cran, caud
	G4VPhysicalVolume* pv_det;

	//messenger
	DetectorMessenger* messenger;
	
};


#endif

