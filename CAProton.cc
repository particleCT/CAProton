#include "EmPhysics_pCT.hh"
#include "HadrontherapyPhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "OrganicMaterial.hh"
#include "G4NistManager.hh"
#include "DetectorConstruction.hh"
#include "SteppingAction.hh"
#include "G4EmCalculator.hh"
#include "EventAction.hh"
#include "G4RunManager.hh"
#include <TROOT.h>
#include "G4Material.hh"
#include "ParallelWorldConstruction.hh"
#include "G4NistManager.hh"
#include "G4Proton.hh"
#include "Analysis.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "OrganicMaterial.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4IonTable.hh"
#ifdef VIS
#include "G4VisExecutive.hh"

#endif
#include <iostream>
#include <fstream>
using namespace std;


void calcRSP(DetectorConstruction*, PrimaryGeneratorAction *);


int main(int argc,char** argv) {
  gROOT->ProcessLine("#include <vector>");
  time_t seed;
  time(&seed);
  CLHEP::RanecuEngine *theRanGenerator = new CLHEP::RanecuEngine;  
  theRanGenerator->setSeed(seed);
  CLHEP::HepRandom::setTheEngine(theRanGenerator);
  G4String paraWorldName = "ParallelWorld";

  G4int nProtons   = atoi(argv[1]);
  G4String phantom = argv[2];
  G4String name    = argv[3];
  G4double energy  = atof(argv[4]);
  G4double angle   = atof(argv[5]);
  G4int thread     = atoi(argv[6]);
  
  G4RunManager* runManager = new G4RunManager;
  runManager->SetUserInitialization(new HadrontherapyPhysicsList(paraWorldName));
  DetectorConstruction* myDC = new DetectorConstruction(phantom,angle);
  myDC->RegisterParallelWorld( new ParallelWorldConstruction(paraWorldName));
  PrimaryGeneratorAction* theGenerator = new PrimaryGeneratorAction(energy);
  runManager->SetUserAction( theGenerator );
  Analysis* theAnalysis    = new Analysis(thread,name);
  runManager->SetUserAction( new SteppingAction() );
  runManager->SetUserAction( new EventAction() );
  runManager->SetUserInitialization( myDC );
  runManager->Initialize();
  
  //G4UImanager * UImanager = G4UImanager::GetUIpointer();
  //UImanager->ApplyCommand("/process/eLoss/fluct false");
  //UImanager->ApplyCommand("/run/verbose 1");
  //UImanager->ApplyCommand("/event/verbose 1");
  //UImanager->ApplyCommand("/tracking/verbose 2");
  //G4UIExecutive * ui = new G4UIExecutive(argc,argv);
  //ui->SessionStart();
#ifdef VIS
  G4UImanager * UImanager = G4UImanager::GetUIpointer();
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
  G4cout << " UI session starts ..." << G4endl;
  G4UIExecutive* ui = new G4UIExecutive(argc, argv);
  UImanager->ApplyCommand("/control/execute vis.mac");
  ui->SessionStart();
#endif
  runManager->BeamOn( nProtons );
  theAnalysis->Save();
  calcRSP(myDC,theGenerator);
  
  delete runManager;
  delete Analysis::GetInstance();
  return 0;
  }

void calcRSP(DetectorConstruction* myDC, PrimaryGeneratorAction* theGenerator){
  G4ParticleDefinition* particle = theGenerator->particle;
  G4EmCalculator* emCal = new G4EmCalculator;

  ofstream myfile;
  myfile.open ("RSP.txt");
  for(size_t i =0;i<myDC->theMaterialList.size();i++){
    G4int HU = myDC->huList[i];
    G4int I = 0;
    G4double tot =0;
    G4Material* water = myDC->water;
    for(int j=10;j<19000;j++){
      G4double dedx_w = emCal->ComputeElectronicDEDX( double(j)/10*MeV,particle,water);
      G4double dedx_b = emCal->ComputeElectronicDEDX( double(j)/10*MeV,particle,myDC->theMaterialList.at(i));
      tot +=dedx_b/dedx_w;
      I+=1;
    }
    G4double RSP = tot/I;
    myfile<<HU<<" "<<myDC->theMaterialList.at(i)->GetDensity()/(g/cm3)<<" "<<RSP<<" "<<myDC->theMaterialList.at(i)->GetName()<<endl;
  }
  myfile.close();

}




