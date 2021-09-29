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
// \author: Haegin Han
//
#include "ParallelPhantom.hh"
#include "TETParameterisation.hh"

#include "G4Box.hh"
#include "G4Tet.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4RunManager.hh"
#include "G4Material.hh"
#include "G4SystemOfUnits.hh"    
#include "G4MultiFunctionalDetector.hh"
#include "G4SDManager.hh"
#include "TETParameterisation.hh"
#include "TETPSEnergyDeposit.hh"
#include "ParallelMessenger.hh"

ParallelPhantom
::ParallelPhantom(G4String parallelWorldName, TETModelImport* _tetData)
:G4VUserParallelWorld(parallelWorldName),fConstructed(false), tetData(_tetData)
{
  messenger = new ParallelMessenger(this);
}

ParallelPhantom::~ParallelPhantom()
{
  delete messenger;
}

void ParallelPhantom::Construct()
{
  if(fConstructed) return;
  fConstructed = true;

  //
  // World
  //
  G4VPhysicalVolume* ghostWorld = GetWorld();
  G4LogicalVolume* worldLogical = ghostWorld->GetLogicalVolume();

  //
  // parallel world placement box
  //
  G4ThreeVector center = (tetData->GetPhantomBoxMax() + tetData->GetPhantomBoxMin())*0.5;
  G4ThreeVector halfSize = (tetData->GetPhantomBoxMax() - tetData->GetPhantomBoxMin())*0.5 + G4ThreeVector(5., 5., 5.)*cm; //5-cm-margin
  G4VSolid* paraBox = new G4Box("phantomBox",halfSize.x(),halfSize.y(),halfSize.z());
  G4LogicalVolume* lv_phantomBox = new G4LogicalVolume(paraBox,0,"phantomBox");
  pv_doctor = new G4PVPlacement(0,center,lv_phantomBox,"phantomBox",worldLogical,false,0);

  //
  // mother of parallel world parameterized volumes
  //
  lv_tet = new G4LogicalVolume(new G4Tet("tet", G4ThreeVector(), G4ThreeVector(0, 0, 1*cm),
                                         G4ThreeVector(0, 1*cm, 0), G4ThreeVector(1*cm, 0, 0)),0,"tet");
  TETParameterisation* param = new TETParameterisation(tetData);
  // new G4PVParameterised("paraPara",lv_tet, lv_phantomBox, kUndefined, tetData->GetNumTetrahedron(), param);
  new G4PVParameterised("param",lv_tet, lv_phantomBox, kUndefined, 1, param);
}

void ParallelPhantom::ConstructSD()
{
  G4MultiFunctionalDetector* mfd = new G4MultiFunctionalDetector("phantom");
  G4SDManager::GetSDMpointer()->AddNewDetector(mfd);
  G4VPrimitiveScorer* scorer = new TETPSEnergyDeposit("edep", tetData);
  mfd->RegisterPrimitive(scorer);
  SetSensitiveDetector(lv_tet, mfd);
}

void ParallelPhantom::Deform(RotationList vQ, Vector3d root)
{
  tetData->Deform(vQ, root);
  G4ThreeVector center = (tetData->GetPhantomBoxMax() + tetData->GetPhantomBoxMin())*0.5;
  G4ThreeVector halfSize = (tetData->GetPhantomBoxMax() - tetData->GetPhantomBoxMin())*0.5 + G4ThreeVector(5., 5., 5.)*cm; //5-cm-margin
  pv_doctor->SetTranslation(center);
  dynamic_cast<G4Box*>(pv_doctor->GetLogicalVolume()->GetSolid())->SetXHalfLength(halfSize.x());
  dynamic_cast<G4Box*>(pv_doctor->GetLogicalVolume()->GetSolid())->SetYHalfLength(halfSize.y());
  dynamic_cast<G4Box*>(pv_doctor->GetLogicalVolume()->GetSolid())->SetZHalfLength(halfSize.z());
 	G4RunManager::GetRunManager()->GeometryHasBeenModified();
}