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
#include "PrimaryMessenger.hh" 

PrimaryGeneratorAction::PrimaryGeneratorAction()
	: G4VUserPrimaryGeneratorAction(), isocenter(1120 * mm, 600 * mm, 1135 * mm)
{
	fPrimary = new G4ParticleGun();
	fPrimary->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("gamma"));

	FlatDetectorInitialization(DetectorZoomField::FD48, 119.5 * cm); // FD, SID
	SetSourceEnergy(80);
	G4double carm_primary = 20 * deg;   // +LAO, -RAO
	G4double carm_secondary = 20 * deg; // +CAU, -CRA
	rotate.rotateY(carm_primary).rotateX(carm_secondary);
	G4ThreeVector focalSpot = rotate * G4ThreeVector(0, 0, -810 * mm);
	fPrimary->SetParticlePosition(focalSpot + isocenter);

	messenger = new PrimaryMessenger(this);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
	delete fPrimary;
	delete messenger;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{

	G4double rand = G4UniformRand();

	G4double rand_energy = cdf.rbegin()->second;
	if (rand < 1)
	{
		for (auto itr : cdf)
		{
			if (rand > itr.first) continue;
			rand_energy = itr.second;
			break;
		}
	}
	fPrimary->SetParticleEnergy(rand_energy);
	fPrimary->SetParticleMomentumDirection(SampleRectangularBeamDirection());
	fPrimary->GeneratePrimaryVertex(anEvent);
}

G4ThreeVector PrimaryGeneratorAction::SampleRectangularBeamDirection()
{
	G4ThreeVector upVec = G4ThreeVector(0, 0, 0) - G4ThreeVector(0, 0, -810);
	G4double a = fabs(upVec.z()) * tan(angle2);
	G4double b = fabs(upVec.z()) * tan(angle1);
	G4double x = 2 * a * G4UniformRand() - a;
	G4double y = 2 * b * G4UniformRand() - b;
	G4double z = upVec.z();

	return rotate * G4ThreeVector(x, y, z);
}

void PrimaryGeneratorAction::SetSourceEnergy(G4int peakE)
{
	G4String fileName(to_string(peakE) + ".spec");
	G4String spectra("./spectra/" + fileName);

	G4cout << "Read x-ray spectra: " << spectra << G4endl;
	ifstream ifs(spectra);

	vector<pair<G4double, G4double>> pdf; //there could be duplicated probabilities, which may cause error in 'map'
	if (!ifs.is_open())
	{
		G4cerr << "X-ray spectra file was not opened" << G4endl;
		exit(1);
	}

	G4double sum(0);
	G4String dump;
	while (getline(ifs, dump))
	{
		stringstream ss(dump);
		ss >> dump;
		if (dump == "Energy[keV]")
		{
			while (getline(ifs, dump))
			{
				G4double energy, intensity;
				stringstream ss2(dump);
				ss2 >> energy >> intensity;
				sum += energy * intensity;
				pdf.push_back(make_pair(energy * intensity, energy * keV));
			}
		}
	}
	ifs.close();

	sort(pdf.begin(), pdf.end(), greater<>());
	transform(pdf.begin(), pdf.end(), pdf.begin(),
			  [&sum](pair<G4double, G4double> iter) -> pair<G4double, G4double> { return make_pair(iter.first / sum, iter.second); });

	cdf.clear();
	G4double sumProb(0);
	for (auto itr : pdf)
	{
		sumProb += itr.first;
		cdf[sumProb] = itr.second;
	}
}

void PrimaryGeneratorAction::FlatDetectorInitialization(DetectorZoomField FD, G4double SID)
{
	switch (FD)
	{
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
