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
// PrimaryMessenger.cc
// \author Haegin Han
//

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "PrimaryMessenger.hh"

PrimaryMessenger::PrimaryMessenger(PrimaryGeneratorAction* _primary)
:G4UImessenger(), fPrimary(_primary), fd(DetectorZoomField::FD16), sid(119.5 * cm)
{
	fBeamDir = new G4UIdirectory("/beam/");
	
	fFDCmd = new G4UIcmdWithAString("/beam/FD", this);
	fFDCmd->SetCandidates("FD48 FD42 FD37 FD31 FD27 FD23 FD19 FD16");
	
	fBeamCmd = new G4UIcmdWith3Vector("/beam/c-arm", this);
	fBeamCmd->SetParameterName("lao[deg]", "caud[deg]", "sid[cm]", true, true);

	fPeakEnergyCmd = new G4UIcmdWithAnInteger("/beam/peak", this);
	fPeakEnergyCmd->SetParameterName("peakE[keV]", false);
}

PrimaryMessenger::~PrimaryMessenger() {
	delete fBeamCmd;
	delete fFDCmd;
	delete fBeamDir;
}

void PrimaryMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
	if(command == fFDCmd){
		if(newValue=="FD48")	    fd = DetectorZoomField::FD48;
		else if(newValue=="FD42")	fd = DetectorZoomField::FD42;
		else if(newValue=="FD37")	fd = DetectorZoomField::FD37;
		else if(newValue=="FD31")	fd = DetectorZoomField::FD31;
		else if(newValue=="FD27")	fd = DetectorZoomField::FD27;
		else if(newValue=="FD23")	fd = DetectorZoomField::FD23;
		else if(newValue=="FD19")	fd = DetectorZoomField::FD19;
		else if(newValue=="FD16")	fd = DetectorZoomField::FD16;

		fPrimary->FlatDetectorInitialization(fd, sid);
	}
	else if(command == fBeamCmd){
		G4ThreeVector input = fBeamCmd->GetNew3VectorValue(newValue);
		sid = input.z()*cm;
		fPrimary->FlatDetectorInitialization(fd, sid);
		fPrimary->SetCarmAngles(input.x()*deg, input.y()*deg);
	}
	else if(command == fPeakEnergyCmd){
		fPrimary->SetSourceEnergy(fPeakEnergyCmd->GetNewIntValue(newValue));
	}
}

