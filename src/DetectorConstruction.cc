#include "TFile.h"
#include "TH3D.h"
#include "TH3S.h"
#include "G4Proton.hh"
#include "G4ThreeVector.hh"
#include "G4PhantomParameterisation.hh"
#include "G4NistManager.hh"
#include "G4EmCalculator.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4RunManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4Trap.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "OrganicMaterial.hh"
#include "Analysis.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "SensitiveDetector.hh"
#include "DetectorConstruction.hh"
#include <set>
#include <stdlib.h> 

DetectorConstruction* DetectorConstruction::theDetector=NULL;
DetectorConstruction::~DetectorConstruction(){theDetector = NULL; }

DetectorConstruction::DetectorConstruction(G4String theName, G4double angle):theAngle(angle),thePhantom(theName)
{
  
  theDetector = this;   
  theOrganicMaterial = OrganicMaterial::GetInstance();
  f      = new TFile(Form("%s",thePhantom.data()),"update");
  hu     = (TH3D*)f->Get("hu");
  rho    = (TH3D*)f->Get("rho");

  NbinsX =  hu->GetNbinsX();
  NbinsY =  hu->GetNbinsY();
  NbinsZ =  hu->GetNbinsZ();
  Xmax   =  hu->GetXaxis()->GetXmax();
  Ymax   =  hu->GetYaxis()->GetXmax();
  Zmax   =  hu->GetZaxis()->GetXmax();
  Xmin   =  hu->GetXaxis()->GetXmin();
  Ymin   =  hu->GetYaxis()->GetXmin();
  Zmin   =  hu->GetZaxis()->GetXmin();
  midX   =  (Xmax+Xmin)/2.;
  midY   =  (Ymax+Ymin)/2.;
  midZ   =  (Zmax+Zmin)/2.;

  halfX  =  hu->GetXaxis()->GetBinWidth(1)/2.;
  halfY  =  hu->GetYaxis()->GetBinWidth(1)/2.;
  halfZ  =  hu->GetZaxis()->GetBinWidth(1)/2.;
  shift  =  G4ThreeVector(midX,midY,midZ);

  G4cout<<"Centroid shift from (0,0,0) ="<<shift<<G4endl;
  cout<<"NbinsX - Xmax - Xmin :"<<NbinsX<<" "<<Xmax<<" "<<Xmin<<endl;
  cout<<"NbinsY - Ymax - Ymin :"<<NbinsY<<" "<<Ymax<<" "<<Ymin<<endl;
  cout<<"NbinsZ - Zmax - Zmin :"<<NbinsZ<<" "<<Zmax<<" "<<Zmin<<endl;


  vector< G4String > MaterialName  = {"Air","LN300","LN450","AP6Adipose","Breast","Water", "CTsolidwater","Brain","Liver","InnerBone","MineralBone","CB230","CB250","CorticalBone"};
  vector< G4double > Threshold     = {-1000, -727,-522, -83, -41 ,0,1,25,83 ,204 , 225   , 455 , 798   , 1203};
  //vector< G4double > Threshold     = {-1000, -727.9,-522.6, -83.8, -41.9 ,0.0,1.0,25,83.7 ,204.0 , 225.8    , 455.8 , 798   , 1203.1};


  //vector< G4String > MaterialName  = {"SinusCavities_LL","SoftTissue_LL","BrainTissue_LL","SpinalDisc_LL","TrabecularBone_LL","CorticalBone_LL","ToothDentin_LL","ToothEnamel_LL"};
  //vector< G4double > Threshold     = {-998, 24, 52, 92, 197, 923, 1280, 2310};

  map<int,double> ids; //A map of the hu and their density counterpart 
  for (G4int k=0;k<NbinsZ;k++) {
    for (G4int j=0;j<NbinsY;j++) {
      for (G4int i=0;i<NbinsX;i++) ids.insert(pair<double,double>(hu->GetBinContent(i+1,j+1,k+1),rho->GetBinContent(i+1,j+1,k+1)) ) ;
    }
  }

  for(auto it=ids.begin();it!=ids.end();++it){

    G4double huUnit   = it->first;
    G4double density  = it->second;

    // Find the closest material based on the density threshold in the material list
    int lower = lower_bound(Threshold.begin(), Threshold.end(), huUnit)-Threshold.begin()-1;
    int mid   = upper_bound(Threshold.begin(), Threshold.end(), huUnit)-Threshold.begin()-1;
    int upper = lower+1;

    //Construct composite material with the relevant density
    if(lower<0) mat = theOrganicMaterial->ConstructMaterial(MaterialName.at(0),density);                                   // Smaller than the smallest point
    else if(lower== int(Threshold.size())-1) mat = theOrganicMaterial->ConstructMaterial(MaterialName.at(lower),density);  // Larger than the largest point
    else if(mid==upper ) mat = theOrganicMaterial->ConstructMaterial(MaterialName.at(mid),density);                        // Intersection Point

    else // In-between
      {
	G4double ratio  = (huUnit - Threshold[lower])/(Threshold[upper]-Threshold[lower]);
	mat = theOrganicMaterial->ConstructCompositeMaterial(MaterialName.at(lower),MaterialName.at(upper),density,ratio); 
      }

    hu2id[huUnit] = theMaterialList.size();
    huList.push_back(huUnit);
    theMaterialList.push_back(mat);
    
  }

  G4cout << ids.size() << " materials " << G4endl;
}


G4VPhysicalVolume* DetectorConstruction::Construct()
{  
  G4cout << "DetectorConstruction::Construct" << G4endl;
  air   = theOrganicMaterial->ConstructMaterial("Air",0.0001025);
  //Cubic Water world
  G4double world_size          = 3*m;  

  G4Box* hallWorld             = new G4Box("World",world_size,world_size,world_size);
  G4LogicalVolume* logicWorld  = new G4LogicalVolume(hallWorld,air,"World");                                    
  G4VPhysicalVolume* physWorld = new G4PVPlacement(0,G4ThreeVector(),logicWorld,"physWorld",0,false,0);
  G4VisAttributes* hall_att    = new G4VisAttributes(G4Colour(1,0,0));
  hall_att->SetVisibility(true);
  logicWorld->SetVisAttributes(hall_att);

  size_t* materialIDs          = new size_t[NbinsX*NbinsY*NbinsZ];
  unsigned int p               = 0;

  for (G4int k=0;k<NbinsZ;k++) {
  for (G4int j=0;j<NbinsY;j++) {
  for (G4int i=0;i<NbinsX;i++) {
    G4int huUnit = hu->GetBinContent(i+1,j+1,k+1);
    materialIDs[p] = hu2id[huUnit];
    p++;
  }
  }
  }

  // Container volume
  G4RotationMatrix* rot1 = new G4RotationMatrix();
  rot1->rotateZ(theAngle*pi/180.);

  G4Box* cont_vol                      = new G4Box("cont_vol",NbinsX*halfX,NbinsY*halfY,NbinsZ*halfZ);
  G4LogicalVolume * cont_log           = new G4LogicalVolume(cont_vol,air,"cont_log",0,0,0);
  G4VPhysicalVolume* cont_phys         = new G4PVPlacement(rot1,G4ThreeVector( (Xmax+Xmin)/2, (Ymax+Ymin)/2, (Zmax+Zmin)/2 ),cont_log,"cont_phys",logicWorld,false,0);
  G4VisAttributes* cont_att = new G4VisAttributes(G4Colour(0,1,1));
  cont_att->SetVisibility(true);
  cont_log->SetVisAttributes(cont_att);

  // Phantom voxel volume 

  G4Box* phantomVoxel_vol              = new G4Box("phantomVoxel_vol",halfX,halfY,halfZ);
  G4LogicalVolume *phantomVoxel_log    = new G4LogicalVolume(phantomVoxel_vol,air,"phantomVoxel_log",0,0,0);

  // Parameterisation of the human phantom

  G4cout << "Make a parameterised phantom" << G4endl;

  param = new G4PhantomParameterisation();
  param->SetVoxelDimensions(halfX,halfY,halfZ);
  param->SetNoVoxel(NbinsX,NbinsY,NbinsZ);
  param->SetMaterials(theMaterialList);
  param->SetMaterialIndices(materialIDs);
  param->CheckVoxelsFillContainer(cont_vol->GetXHalfLength(),cont_vol->GetYHalfLength(),cont_vol->GetZHalfLength());
  param->SetSkipEqualMaterials(false);
  param->BuildContainerSolid(cont_phys);
  param->CheckVoxelsFillContainer(cont_vol->GetXHalfLength(),cont_vol->GetYHalfLength(),cont_vol->GetZHalfLength());
  phantomPhys = new G4PVParameterised("phantomPhys",phantomVoxel_log,cont_log,kZAxis,param->GetNoVoxel(),param);
  phantomPhys->SetRegularStructureId(1);

  return physWorld;

}


