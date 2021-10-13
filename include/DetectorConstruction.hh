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
    void SetTablePose(G4ThreeVector table_trans, G4double table_pivot_angle) 
	//_table_trans is the translation from when rotation center is at the center of the table
	{  
		G4ThreeVector frameOrigin = frame_ralative_to_table + table_default_trans + table_trans; //before rotation
		if(table_pivot_angle==0) 
		{
			pv_frame->SetTranslation(frameOrigin);
			return;
		}

		frame_rotation_matrix->set(G4RotationMatrix::IDENTITY.axisAngle());
		frame_rotation_matrix->rotateZ(-table_pivot_angle);
		pv_frame->SetTranslation(frameOrigin + frame_rotation_matrix->inverse() * (table_rotation_center - frameOrigin));
	}

    // Glass
    void SetGlassPose(G4ThreeVector _glass_trans, G4ThreeVector _glass_axis, G4double _glass_theta) 
	{
		pv_glass->SetTranslation(_glass_trans);
		pv_glass->GetRotation()->setAxis(_glass_axis);
		pv_glass->GetRotation()->setTheta(-_glass_theta);
    }

    // C-arm det
    void SetCarmDetPose(G4double _carm_primary, G4double _carm_secondary, G4double _carm_sid) 
	// 	carm_primary: rao, lao | carm_secondary: cran, caud
	{
		carm_rotation_matrix->setTheta(0); //set identity
		carm_rotation_matrix->rotateY(_carm_primary).rotateX(_carm_secondary);
		carm_rotation_matrix->invert();
		pv_det->SetTranslation(carm_rotation_matrix->inverse()*G4ThreeVector(0,0,_carm_sid-focalLength) + carm_isocenter);
    }


private:
    void SetupWorldGeometry();

	void ConstructOperatingTable();
	G4LogicalVolume* ConstructPatient();
	void ConstructCarmDet();
	void ConstructPbGlass();

	G4LogicalVolume*   worldLogical;
	G4VPhysicalVolume* worldPhysical;

	// Operating Table
	G4VPhysicalVolume* pv_frame;
	G4ThreeVector table_rotation_center, table_default_trans; //program input (relative coordinate to ChArUco)
	G4RotationMatrix* frame_rotation_matrix; //inverse
	G4ThreeVector frame_ralative_to_table;

	// Glass
	G4VPhysicalVolume* pv_glass;

	// C-arm
	G4VPhysicalVolume* pv_det;
	G4ThreeVector carm_isocenter;
	G4RotationMatrix* carm_rotation_matrix; //inverse
	G4double focalLength;

	//messenger
	DetectorMessenger* messenger;
	
};


#endif

