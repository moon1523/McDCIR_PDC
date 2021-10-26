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
#include <map>
#include <algorithm>
#include "TETModelImport.hh"

using namespace std;

class PrimaryMessenger;

enum DetectorZoomField
{
  FD48,
  FD42,
  FD37,
  FD31,
  FD27,
  FD23,
  FD19,
  FD16
};

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction();
  virtual ~PrimaryGeneratorAction();

  virtual void GeneratePrimaries(G4Event *);

  G4ParticleGun *GetParticleGun() const { return fPrimary; }
  void SetSourceEnergy(G4int peakE); //peakE in keV

  void FlatDetectorInitialization(DetectorZoomField FD, G4double SID);
  void SetCarmAngles(G4double primary, G4double secondary)
  // carm_primary = 20 * deg;   // +LAO, -RAO
  // carm_secondary = 20 * deg; // +CAU, -CRA
  {
    rotate.set(G4RotationMatrix::IDENTITY.axisAngle());
    rotate.rotateY(primary).rotateX(secondary);

    G4ThreeVector focalSpot = rotate * G4ThreeVector(0, 0, -810*mm); //what is 810?
    fPrimary->SetParticlePosition(focalSpot + isocenter);
  }

  G4ThreeVector SampleRectangularBeamDirection();

private:
  G4ParticleGun *fPrimary;
  G4double angle1, angle2;
  G4RotationMatrix rotate;
  G4ThreeVector isocenter;

  // Energy
  map<G4double, G4double> cdf;

  //messenger
  PrimaryMessenger *messenger;
};

#endif
