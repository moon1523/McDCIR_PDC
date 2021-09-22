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
// DetectorMessenger.cc
// \author Haegin Han
//

#include "G4UIdirectory.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4RunManager.hh"
#include <sstream>
#include <vector>
#include "RunAction.hh"
#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"

DetectorMessenger::DetectorMessenger(DetectorConstruction* _det)
:G4UImessenger(), fDet(_det)
{
	fMachineDir = new G4UIdirectory("/machine/");
	fTableTransCmd = new G4UIcmdWith3VectorAndUnit("/machine/tableTrans", this);
	fTablePivotCmd = new G4UIcmdWithADoubleAndUnit("/mahcine/pivot", this);
	fDetCmd = new G4UIcmdWith3Vector("/machine/c-arm", this);
	fGlassTransCmd = new G4UIcmdWith3VectorAndUnit("/machine/glassTrans", this);
	fGlassRotCmd = new G4UIcmdWith3Vector("/machine/glassRot", this);
	fCloseCmd = new G4UIcmdWithoutParameter("/machine/close", this);

	fTableTransCmd->SetDefaultUnit("cm");
	fTablePivotCmd->SetDefaultUnit("deg");
	fGlassTransCmd->SetDefaultUnit("cm");

	fDetCmd->SetParameterName("lao[deg]", "caud[deg]", "sid[cm]", true, true);
	fGlassRotCmd->SetParameterName("axisX*[deg]", "axisY*[deg]", "axisZ*[deg]", false);

}

DetectorMessenger::~DetectorMessenger() {
	delete fMachineDir;
	delete fTableTransCmd; //trans
	delete fTablePivotCmd; //pivot
	delete fDetCmd; //primary, secondary, sid
	delete fGlassTransCmd;
	delete fGlassRotCmd; // axis * angle(in deg)
	delete fCloseCmd;
}

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
	if(command == fTableTransCmd){
		tableTrans = fTableTransCmd->GetNew3VectorValue(newValue);
		fDet->SetTablePose(tableTrans, tablePivot);
	}
	else if(command == fTablePivotCmd){
		tablePivot = fTablePivotCmd->GetNewDoubleValue(newValue);
		fDet->SetTablePose(tableTrans, tablePivot);
	}
	else if(command == fDetCmd){
		G4ThreeVector det = fDetCmd->GetNew3VectorValue(newValue);
		fDet->SetCarmDetPose(det.x()*deg, det.y()*deg, det.z()*cm);
	}
	else if(command == fGlassTransCmd){
		glassTrans = fGlassTransCmd->GetNew3VectorValue(newValue);
		fDet->SetGlassPose(glassTrans, glassAxis, glassTheta);
	}
	else if(command == fGlassRotCmd){
		G4ThreeVector rot = fGlassRotCmd->GetNew3VectorValue(newValue);
		glassTheta = rot.mag() * deg;
		glassAxis = rot.unit();
		fDet->SetGlassPose(glassTrans, glassAxis, glassTheta);
	}
	else if(command == fCloseCmd){
		G4RunManager::GetRunManager()->GeometryHasBeenModified();
	}
}

