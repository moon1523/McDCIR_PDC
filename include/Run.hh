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


#ifndef Run_h
#define Run_h 1

#include "G4Run.hh"
#include "G4Event.hh"
#include "G4THitsMap.hh"
#include "G4SDManager.hh"

#include "TETModelImport.hh"

typedef std::map<G4int, std::pair<G4double, G4double>> EDEPMAP;

class Run : public G4Run
{
public:
	Run(TETModelImport* tetData);
	virtual ~Run();

	virtual void RecordEvent(const G4Event*);
    virtual void Merge(const G4Run*);

    EDEPMAP* GetEdepMap() {return &edepMap;}
    G4String GetParticleName() {return primary;}
    G4double GetParticleEnergy()   {return primaryE;}

private:
	EDEPMAP edepMap;
	G4int fCollID;
	G4int fCollID_DRF;

	std::map<G4int, std::vector<G4int>>  organ2dose;
	std::map<G4int, G4double>  rbmFactor;
	std::map<G4int, G4double>  bsFactor;
	G4bool doseOrganized;

	G4String primary;
	G4double primaryE;
};

#endif
