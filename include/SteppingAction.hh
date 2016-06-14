#ifndef SteppingAction_hh
#define SteppingAction_hh

#include "globals.hh"
#include "DetectorConstruction.hh"
#include "G4UserSteppingAction.hh"
#include <map>

class TH3S;
class DetectorConstruction;
class SteppingAction : public G4UserSteppingAction {
public:
  SteppingAction();
  ~SteppingAction();
  void UserSteppingAction(const G4Step*);
  static inline SteppingAction* GetInstance() { return theSteppingAction; }
  std::map<G4int,G4double> Edep;
  G4double Einit;
  std::vector<double> temp_X;
  std::vector<double> temp_Y;
  std::vector<double> temp_Z;
  std::vector<double> Eloss;
  std::vector<double> Length;
  std::vector<TString> temp_name;

private:  
  static SteppingAction* theSteppingAction;
  vector< pair<G4int,G4double> >::iterator it,it_end;
  pair<map<G4int,double>::iterator,bool> ret;
};

#endif



