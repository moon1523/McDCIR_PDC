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

#include "PrimaryGeneratorAction.hh"
#include "G4ParticleTable.hh"

enum DetectorZoomField
{
	FD48, FD42, FD37, FD31, FD27, FD23, FD19, FD16
};

PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction()
{
	fPrimary = new G4ParticleGun();
	fPrimary->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("gamma"));

	FlatDetectorInitialization(DetectorZoomField::FD48, 119.5 * cm); // FD, SID
	SetSourceEnergy();
	G4ThreeVector translate_from_origin(1120,600,1135);
	carm_primary   = 20 * deg;  // +LAO, -RAO
	carm_secondary = 20 * deg;  // +CAU, -CRA
	rotate.rotateX(carm_secondary).rotateY(carm_primary);
	G4ThreeVector focalSpot = rotate * G4ThreeVector(0,0,-810);
	fPrimary->SetParticlePosition(focalSpot + translate_from_origin);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
	delete fPrimary;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

	G4double rand_energy = G4UniformRand();

	if (rand_energy == 1)	rand_energy = pdf_sort.rbegin()->second;
	else {
		for (auto itr : cdf_sort) {
			if (rand_energy < itr.first) {
				rand_energy = itr.second;
				break;
			}
		}
	}
	fPrimary->SetParticleEnergy(rand_energy*keV);
	fPrimary->SetParticleMomentumDirection(SampleRectangularBeamDirection());
	fPrimary->GeneratePrimaryVertex(anEvent);
}

G4ThreeVector PrimaryGeneratorAction::SampleRectangularBeamDirection()
{
	G4ThreeVector upVec = G4ThreeVector(0,0,0) - G4ThreeVector(0,0,-810);
	G4double a = fabs(upVec.z()) * tan(angle2);
	G4double b = fabs(upVec.z()) * tan(angle1);
	G4double x = 2*a*G4UniformRand() - a;
	G4double y = 2*b*G4UniformRand() - b;
	G4double z = upVec.z();

	return rotate * G4ThreeVector(x,y,z);
}

void PrimaryGeneratorAction::SetSourceEnergy()
{
	pdf.clear();
	pdf_sort.clear();
	cdf_sort.clear();

	peak_energy = 80;
	G4String fileName(to_string(peak_energy) + ".spec");
	G4String spectra("../spectra/" + fileName);

	G4cout << "Read x-ray spectra: " << spectra << G4endl;
	ifstream ifs(spectra);
	if(!ifs.is_open()) { G4cerr << "X-ray spectra file was not opened" << G4endl; exit(1); }

	G4double sum(0);
	G4String dump;
	while(getline(ifs,dump)) {
		stringstream ss(dump);
		ss >> dump;
		if (dump == "Energy[keV]") {
			while (getline(ifs,dump)) {
				G4double energy, intensity;
				stringstream ss2(dump);
				ss2 >> energy >> intensity;
				sum += energy * intensity;
				pdf[energy] = energy*intensity;
			}
		}
	}

	for (auto &itr : pdf) {
		itr.second /= sum;
		pdf_sort[itr.second] = itr.first;
	}

	G4double cdf(0);
	for (auto itr : pdf_sort) {
		cdf += itr.first;
		cdf_sort[cdf] = itr.second;
	}

	ifs.close();
}

void PrimaryGeneratorAction::FlatDetectorInitialization(G4int FD, G4double SID)
{
	switch(FD) {
	case DetectorZoomField::FD48:
		angle1 = atan((38 * cm) * 0.5 / SID);
		angle2 = atan((30 * cm) * 0.5 / SID);
		break;
	case DetectorZoomField::FD42:
		angle1 = atan((30 * cm) * 0.5 / SID);
		angle2 = angle1;
		break;
	case DetectorZoomField::FD37:
		angle1 = atan((26 * cm) * 0.5 / SID);
		angle2 = angle1;
		break;
	case DetectorZoomField::FD31:
		angle1 = atan((22 * cm) * 0.5 / SID);
		angle2 = angle1;
		break;
	case DetectorZoomField::FD27:
		angle1 = atan((19 * cm) * 0.5 / SID);
		angle2 = angle1;
		break;
	case DetectorZoomField::FD23:
		angle1 = atan((16 * cm) * 0.5 / SID);
		angle2 = angle1;
		break;
	case DetectorZoomField::FD19:
		angle1 = atan((13.5 * cm) * 0.5 / SID);
		angle2 = angle1;
		break;
	case DetectorZoomField::FD16:
		angle1 = atan((11 * cm) * 0.5 / SID);
		angle2 = angle1;
		break;
	}
}
