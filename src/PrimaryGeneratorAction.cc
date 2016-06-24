#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "G4EmCalculator.hh"
#include "G4NistManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4UImanager.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include "G4Proton.hh"
#include "G4IonTable.hh"
#include <math.h>
#include "G4NavigationHistory.hh"
#include "G4TouchableHistory.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TransportationManager.hh"
#include "TFile.h"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
using namespace std;

PrimaryGeneratorAction* PrimaryGeneratorAction::theGenerator = NULL;

PrimaryGeneratorAction::PrimaryGeneratorAction(G4double energy):ENER(energy)
{ 
  theGenerator = this;
  theDetector = DetectorConstruction::GetInstance();

  particleGun = new G4ParticleGun();
  nProtonsGenerated = 0;
  IrradiatedEnergy  = 0.;

  beamPosY   = linspace(theDetector->Ymin,theDetector->Ymax,0.8);
  beamPosZ   = linspace(theDetector->Zmin,theDetector->Zmax,0.8);

  fieldSizeY = 2*(theDetector->halfY*theDetector->NbinsY);
  fieldSizeZ = 2*(theDetector->halfZ*theDetector->NbinsZ);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  theGenerator = NULL;
  delete particleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  particle = G4IonTable::GetIonTable()->GetIon(2,4,0);
  particleGun->SetParticleDefinition(particle);
  Einit = ENER*MeV;
  
  /*G4double cosTheta = (1./5.)*(G4UniformRand()-0.5);
  phi = (1./16.)*(pi*G4UniformRand()-pi/2);
  G4double sinTheta = std::sqrt(1. - cosTheta*cosTheta);*/
  px0 = 1;//sinTheta*std::cos(phi);
  py0 = 0;//sinTheta*std::sin(phi);
  pz0 = 0;//cosTheta;

  x0 = theDetector->Xmin - 0.5*mm;
  y0 = G4UniformRand()*fieldSizeY-fieldSizeY/2 + theDetector->shift.y();     
  z0 = G4UniformRand()*fieldSizeZ-fieldSizeZ/2 + theDetector->shift.z();

  Position = G4ThreeVector(x0,y0,z0); 
  Momentum = G4ThreeVector(px0,py0,pz0);

  particleGun->SetParticleEnergy(Einit);
  particleGun->SetParticleMomentumDirection(Momentum);
  particleGun->SetParticlePosition(Position);
  particleGun->GeneratePrimaryVertex(anEvent);
  IrradiatedEnergy += Einit;
  nProtonsGenerated++;
  if(nProtonsGenerated%5000==0) cout << nProtonsGenerated << endl;
}

vector<G4double> PrimaryGeneratorAction::linspace(double a, double b, double step) {
  vector<double> array;
  while(a <= b) {
    array.push_back(a);
    a += step;
  }
  return array;
}

