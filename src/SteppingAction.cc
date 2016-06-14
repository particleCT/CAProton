#include "TH3S.h"
#include <vector>
#include <iostream>
#include <iomanip>
#include "G4Step.hh"
#include "G4VProcess.hh"
#include "G4RegularNavigationHelper.hh"
#include "SteppingAction.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

SteppingAction *SteppingAction::theSteppingAction=NULL;
SteppingAction::~SteppingAction() 
{theSteppingAction=NULL;}


SteppingAction::SteppingAction()
{
  theSteppingAction=this;
  Einit = 330.0;
}


void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4Track *tr = aStep->GetTrack();
  if ( tr->GetTrackID()!=1 ){
    tr->SetTrackStatus(fStopAndKill);
    temp_X.clear();
    temp_Y.clear();
    temp_Z.clear();
    temp_name.clear();
    Eloss.clear();
    return;
  }
  if(tr->GetTrackID()==1){
  G4ThreeVector preStepPos = aStep->GetPreStepPoint()->GetPosition();                                                                                                               
  G4double en  = aStep->GetTotalEnergyDeposit();
  Eloss.push_back(en);
  Einit-=en;
  cout<<"Pre  "<<aStep->GetPreStepPoint()->GetKineticEnergy()<<endl;
  cout<<"Post "<<aStep->GetPostStepPoint()->GetKineticEnergy()<<endl;
  cout<<"Energy deposity "<<en<<endl;
  cout<<Einit<<endl;
  temp_X.push_back(preStepPos.x());
  temp_Y.push_back(preStepPos.y());
  temp_Z.push_back(preStepPos.z());
  temp_name.push_back(tr->GetMaterial()->GetName()); 
  }
}

  



