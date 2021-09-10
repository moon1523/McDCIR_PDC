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

#ifndef PrimaryGeneratorAction_hh
#define PrimaryGeneratorAction_hh 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "G4RotationMatrix.hh"
#include "G4RandomDirection.hh"
#include "TETModelImport.hh"
#include "G4AffineTransform.hh"

#include <map>
#include <algorithm>

using namespace std;

class Messenger;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction();
    virtual ~PrimaryGeneratorAction();

    virtual void GeneratePrimaries(G4Event*);

    // Access Function
    //
	G4ParticleGun*  GetParticleGun()            const {return fPrimary;}
	void SetSourceEnergy();

	void FlatDetectorInitialization(G4int FD, G4double SID);
	G4ThreeVector SampleRectangularBeamDirection();



  private:
    G4ParticleGun*    fPrimary;
    G4double carm_primary, carm_secondary;
    G4double angle1, angle2;
    G4RotationMatrix rotate;

    // Energy
    map<G4double, G4double> pdf;
    map<G4double, G4double, greater<G4double>> pdf_sort;
    map<G4double, G4double> cdf_sort;
    G4int peak_energy;

};

#endif
