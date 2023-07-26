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
/// \file optical/LXe/src/LXeMainVolume.cc
/// \brief Implementation of the LXeMainVolume class
//
//

//#define NOTOFU 1

#include "LXeMainVolume.hh"

#include "globals.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalSurface.hh"
#include "G4Sphere.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"

#include "G4NistManager.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4PhysicalConstants.hh"
#include "G4SubtractionSolid.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4UnionSolid.hh"
#include <fstream>
using namespace std;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeMainVolume::LXeMainVolume(G4RotationMatrix* pRot, const G4ThreeVector& tlate,
                             G4LogicalVolume* pMotherLogical, G4bool pMany,
                             G4int pCopyNo, LXeDetectorConstruction* c)
  // Pass info to the G4PVPlacement constructor
  : G4PVPlacement(pRot, tlate,
                  // Temp logical volume must be created here
                  new G4LogicalVolume(new G4Box("temp", 1, 1, 1),
                                      G4Material::GetMaterial("Vacuum"), "temp",
                                      0, 0, 0),
                  "housing", pMotherLogical, pMany, pCopyNo)
  , fConstructor(c)
{


    G4VisAttributes* Red=
            new G4VisAttributes(G4Colour(1.0,0.0,0.0,0.3));  // Red
    Red->SetVisibility(true);
    Red->SetForceSolid(true);

    G4VisAttributes* Green=
            new G4VisAttributes(G4Colour(0.0,1.0,0.0,0.3));  // Green
    Green->SetVisibility(true);
    Green->SetForceSolid(true);

    G4VisAttributes* Blue=
            new G4VisAttributes(G4Colour(0.0,0.0,1.0,0.3));  // Blue
    Blue->SetVisibility(true);
    Blue->SetForceSolid(true);

    G4VisAttributes* BlueL=
            new G4VisAttributes(G4Colour(0.0,0.0,1.0,0.1));  // Blue
    Blue->SetVisibility(true);
    Blue->SetForceSolid(true);
//    Blue->SetForceAuxEdgeVisible(true);
//    Blue->SetForceSolid(false);
//    Blue->SetForceWireframe(true);

    G4VisAttributes* Yellow=
            new G4VisAttributes(G4Colour(1.0,1.0,0.0,0.3));  // Yellow
    Yellow->SetVisibility(true);
    Yellow->SetForceSolid(true);

    G4VisAttributes* Turquoise=
            new G4VisAttributes(G4Colour(0.0,0.8,0.8,0.3));  // Turquoise
    Turquoise->SetVisibility(true);
    Turquoise->SetForceSolid(true);

    G4VisAttributes* Magenta=
            new G4VisAttributes(G4Colour(1.0,0.0,1.0,0.3));  // Magenta
    Magenta->SetVisibility(true);
    Magenta->SetForceSolid(true);

    G4VisAttributes* INVI=
            new G4VisAttributes(G4Colour(1.0,0.0,1.0,0.3));  // INVI
    INVI->SetVisibility(false);

    G4VisAttributes* Grey=
            new G4VisAttributes(G4Colour(0.5,0.5,0.5,0.3));  // Gray


    G4VisAttributes* wireFrame=
            new G4VisAttributes(G4Colour(1.0,0.0,0.0));  // Gray
    wireFrame->SetVisibility(true);
    wireFrame->SetForceSolid(false);
    wireFrame->SetForceWireframe(true);


    Grey->SetVisibility(true);
    Grey->SetForceSolid(true);
    ColorMap.emplace("Red",Red);
    ColorMap.emplace("Green",Green);
    ColorMap.emplace("Blue",Blue);
    ColorMap.emplace("Yellow",Yellow);
    ColorMap.emplace("Turquoise",Turquoise);
    ColorMap.emplace("Magenta",Magenta);
    ColorMap.emplace("Grey",Grey);
    ColorMap.emplace("invi",INVI);
    ColorMap.emplace("BlueL",BlueL);
    ColorMap.emplace("wireFrame",wireFrame);

  CopyValues();

  fExpHall_box =
    new G4Box("ExpHall_box", fExpHall_x / 2., fExpHall_y / 2., fExpHall_z / 2.);
  fExpHall_log = new G4LogicalVolume(
    fExpHall_box, G4Material::GetMaterial("Air"), "ExpHall_log", 0, 0, 0);
  fExpHall_log->SetVisAttributes(ColorMap["invi"]);

  //*************************** housing and scintillator
//  fScint_box =
//    new G4Box("scint_box", fScint_x / 2., fScint_y / 2., fScint_z / 2.);
//  fHousing_box =
//    new G4Box("housing_box", housing_x / 2., housing_y / 2., housing_z / 2.);


  fDetContainer_log = BuildPlaLogic(0);


  //*************** Miscellaneous sphere to demonstrate skin surfaces
  fSphere = new G4Sphere("sphere", 0., 2. * cm, 0. * deg, 360. * deg, 0. * deg,
                         360. * deg);
  fSphere_log =
    new G4LogicalVolume(fSphere, G4Material::GetMaterial("Al"), "sphere_log");
  if(fSphereOn)
    new G4PVPlacement(0, G4ThreeVector(5. * cm, 5. * cm, 5. * cm), fSphere_log,
                      "sphere", fScint_log, false, 0);

  //****************** Build PMTs
  G4double height_pmt = 32.5 * mm;

  G4double fOuterRadius_pmthousing =  30*mm;
  G4double fOuterRadius_glass =  26.2*mm;

  G4double height_photocathode = fD_mtl;

  fPmt = new G4Box("pmt_tube", fOuterRadius_pmthousing/2, fOuterRadius_pmthousing/2, height_pmt/2);
  G4Box* fPmt_glass = new G4Box("pmt_tube_glasswindow", fOuterRadius_glass/2, fOuterRadius_glass/2, height_photocathode);
  // the "photocathode" is a metal slab at the back of the glass that
  // is only a very rough approximation of the real thing since it only
  // absorbs or detects the photons based on the efficiency set below
  fPhotocath = new G4Box("photocath_tube", fOuterRadius_pmt/2, fOuterRadius_pmt/2, height_photocathode/2);

  fPmt_log =
    new G4LogicalVolume(fPmt, G4Material::GetMaterial("Al"), "pmt_log");

  G4LogicalVolume* Pmt_log_glass =
    new G4LogicalVolume(fPmt_glass, G4Material::GetMaterial("Borosilicate_Glass"), "pmt_log");
//  fPmt_log =
//    new G4LogicalVolume(fPmt, G4Material::GetMaterial("Glass"), "pmt_log");
  fPhotocath_log = new G4LogicalVolume(
    fPhotocath, G4Material::GetMaterial("Bialkali"), "photocath_log");


//  fPhotocath_log = new G4LogicalVolume(
//    fPhotocath, G4Material::GetMaterial("Al"), "photocath_log");

  new G4PVPlacement(0, G4ThreeVector(0., 0., -height_photocathode/2), fPhotocath_log,
                    "photocath", Pmt_log_glass, false, 0);

  new G4PVPlacement(0, G4ThreeVector(0., 0., +height_pmt / 2.-height_photocathode), Pmt_log_glass,
                    "photocath_glass", fPmt_log, false, 0);

  Pmt_log_glass->SetVisAttributes(ColorMap["wireFrame"]);
  fPmt_log->SetVisAttributes(ColorMap["Blue"]);
  fPhotocath_log->SetVisAttributes(ColorMap["Grey"]);
  //***********Arrange scintillator and Pmt**********


  // Scintilator
  nofTubes = 0;
  G4int tubeGroup[1000];
  G4double x[1000];
  G4double y[1000];
  G4double z[1000];
  G4double rot[1000];
  ifstream infile("detmap.txt");
  while (!infile.eof()){
      G4int myid;
      infile>>myid;
      infile>>x[myid]>>y[myid]>>z[myid]>>rot[myid]>>tubeGroup[myid];
      nofTubes++;
//      G4cout<<myid<<G4endl;
  }
  G4cout<<"Read "<<nofTubes<<G4endl;
  infile.close();

  for (G4int i=0;i<nofTubes;i++)
  {
      //--- Define Detector Group --
      //--------------------------------------------------------------
      G4RotationMatrix* rotmat;
      rotmat = new G4RotationMatrix;
      rotmat->set(0,0,0);
      rotmat->rotateZ(90*degree-rot[i]*degree);
#ifndef NOTOFU
      new G4PVPlacement(rotmat,
                        G4ThreeVector(x[i]*mm,y[i]*mm,z[i]*mm),             //position
                        fDetContainer_log,            //its logical volume
                        "PhysicsDet",             //its name
                        fExpHall_log,             //its mother  volume
                        false,                 //no boolean operation
                        i,                 //copy number
                        fCheckOverlaps);       // checking overlaps
      fScintPositions.push_back(G4ThreeVector(x[i]*mm,y[i]*mm,z[i]*mm));
#endif
  }



  // PMT

  G4double detZ = 800*mm;

  G4int k = 0;
  for (G4int i=0;i<nofTubes;i++)
  {
      G4RotationMatrix* rotmat;
      rotmat = new G4RotationMatrix;
      rotmat->set(0,0,0);
      rotmat->rotateZ(90*degree-rot[i]*degree);
      G4double pmt_posX = x[i]*mm;
      G4double pmt_posY = y[i]*mm;
      G4double pmt_posZ_down = z[i]*mm-detZ/2-height_pmt/2;

#ifndef NOTOFU
      new G4PVPlacement(rotmat, G4ThreeVector(pmt_posX, pmt_posY, pmt_posZ_down), fPmt_log, "pmt",
                        fExpHall_log, false, k,fCheckOverlaps);
      fPmtPositions.push_back(G4ThreeVector(pmt_posX, pmt_posY, pmt_posZ_down));
#endif
      k++;

      G4RotationMatrix* rotmat2;
      rotmat2 = new G4RotationMatrix;
      rotmat2->set(0,0,0);

      rotmat2->rotateZ(90*degree-rot[i]*degree);
      rotmat2->rotateY(180*degree);

      G4double pmt_posZ_up = z[i]*mm+detZ/2+height_pmt/2;

#ifndef NOTOFU
      new G4PVPlacement(rotmat2, G4ThreeVector(pmt_posX, pmt_posY, pmt_posZ_up), fPmt_log, "pmt",
                        fExpHall_log, false, k,fCheckOverlaps);
      fPmtPositions.push_back(G4ThreeVector(pmt_posX, pmt_posY, pmt_posZ_up));
#endif
      k++;

  }


  // Single detector
//  new G4PVPlacement(0, G4ThreeVector(-10*cm,-10*cm,0*cm), fDetContainer_log, "scintillator",
//                    fExpHall_log, false, 0);

//  new G4PVPlacement(0, G4ThreeVector(-10*cm,-10*cm, -detZ/2-height_pmt/2), fPmt_log, "pmt",
//                    fExpHall_log, false, 0,fCheckOverlaps);
//  fPmtPositions.push_back(G4ThreeVector(0, 0, -detZ/2-height_pmt/2));
//  G4RotationMatrix* rotmat2;
//  rotmat2 = new G4RotationMatrix;
//  rotmat2->set(0,0,0);
//  rotmat2->rotateZ(90*degree);
//  rotmat2->rotateY(180*degree);
//  new G4PVPlacement(rotmat2, G4ThreeVector(-10*cm,-10*cm, detZ/2+height_pmt/2), fPmt_log, "pmt",
//                    fExpHall_log, false, 1,fCheckOverlaps);
//  fPmtPositions.push_back(G4ThreeVector(0, 0, detZ/2+height_pmt/2));



  VisAttributes();
  SurfaceProperties();

  SetLogicalVolume(fExpHall_log);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeMainVolume::CopyValues()
{
    fExpHall_x         = fConstructor->GetExpHallX();
    fExpHall_y         = fConstructor->GetExpHallY();
    fExpHall_z         = fConstructor->GetExpHallZ();

  fScint_x         = fConstructor->GetScintX();
  fScint_y         = fConstructor->GetScintY();
  fScint_z         = fConstructor->GetScintZ();
  fD_mtl           = fConstructor->GetHousingThickness();
  fNx              = fConstructor->GetNX();
  fNy              = fConstructor->GetNY();
  fNz              = fConstructor->GetNZ();
  fOuterRadius_pmt = fConstructor->GetPMTRadius();
  fSphereOn        = fConstructor->GetSphereOn();
  fRefl            = fConstructor->GetHousingReflectivity();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeMainVolume::PlacePMTs(G4LogicalVolume* pmt_log, G4RotationMatrix* rot,
                              G4double& a, G4double& b, G4double da,
                              G4double db, G4double amin, G4double bmin,
                              G4int na, G4int nb, G4double& x, G4double& y,
                              G4double& z, G4int& k)
{
  /*  PlacePMTs : a different way to parameterize placement that does not depend
   * on calculating the position from the copy number
   *
   *  pmt_log = logical volume for pmts to be placed
   *  rot = rotation matrix to apply
   *  a,b = coordinates to vary(ie. if varying in the xy plane then pass x,y)
   *  da,db = value to increment a,b by
   *  amin,bmin = start values for a,b
   *  na,nb = number of repitions in a and b
   *  x,y,z = just pass x,y, and z by reference (the same ones passed for a,b)
   *  k = copy number to start with
   *  sd = sensitive detector for pmts
   */
  a = amin;
  for(G4int j = 1; j <= na; ++j)
  {
    a += da;
    b = bmin;
    for(G4int i = 1; i <= nb; ++i)
    {
      b += db;
      new G4PVPlacement(rot, G4ThreeVector(x, y, z), pmt_log, "pmt",
                        fHousing_log, false, k);
      fPmtPositions.push_back(G4ThreeVector(x, y, z));
      ++k;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeMainVolume::VisAttributes()
{
  G4VisAttributes* housing_va = new G4VisAttributes(G4Colour(0.8, 0.8, 0.8));
  fHousing_log->SetVisAttributes(housing_va);

  G4VisAttributes* sphere_va = new G4VisAttributes();
  sphere_va->SetForceSolid(true);
  fSphere_log->SetVisAttributes(sphere_va);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeMainVolume::SurfaceProperties()
{
  std::vector<G4double> ephoton = { 2.004*eV, 2.615*eV, 2.760*eV, 2.922*eV, 3.191*eV };

  //**Scintillator housing properties
  std::vector<G4double> reflectivity     = { fRefl, fRefl,fRefl, fRefl,fRefl };
  std::vector<G4double> efficiency       = { 0.0, 0.0,0.0, 0.0,0.0 };
  G4MaterialPropertiesTable* scintHsngPT = new G4MaterialPropertiesTable();
  scintHsngPT->AddProperty("REFLECTIVITY", ephoton, reflectivity);
  scintHsngPT->AddProperty("EFFICIENCY", ephoton, efficiency);
  G4OpticalSurface* OpScintHousingSurface =
    new G4OpticalSurface("HousingSurface", unified, polished, dielectric_metal);
  OpScintHousingSurface->SetMaterialPropertiesTable(scintHsngPT);

  //**Sphere surface properties
  std::vector<G4double> sphereReflectivity = { 1.0, 1.0,1.0, 1.0,1.0 };
  std::vector<G4double> sphereEfficiency   = { 0.0, 0.0,0.0, 0.0,0.0 };
  G4MaterialPropertiesTable* spherePT      = new G4MaterialPropertiesTable();
  spherePT->AddProperty("REFLECTIVITY", ephoton, sphereReflectivity);
  spherePT->AddProperty("EFFICIENCY", ephoton, sphereEfficiency);
  G4OpticalSurface* OpSphereSurface =
    new G4OpticalSurface("SphereSurface", unified, polished, dielectric_metal);
  OpSphereSurface->SetMaterialPropertiesTable(spherePT);

  //**Photocathode surface properties
  std::vector<G4double> photocath_EFF     = { 1., 1.,1., 1.,1. };
  std::vector<G4double> photocath_ReR     = { 1.92, 1.92,1.92, 1.92,1.92 };
  std::vector<G4double> photocath_ImR     = { 1.69, 1.69,1.69, 1.69,1.69 };
  std::vector<G4double>  photocath_Reflectivity = { 0., 0.,0., 0.,0.};
  G4MaterialPropertiesTable* photocath_mt = new G4MaterialPropertiesTable();
  photocath_mt->AddProperty("EFFICIENCY", ephoton, photocath_EFF);
  photocath_mt->AddProperty("REFLECTIVITY", ephoton, photocath_Reflectivity);
  photocath_mt->AddProperty("REALRINDEX", ephoton, photocath_ReR);
  photocath_mt->AddProperty("IMAGINARYRINDEX", ephoton, photocath_ImR);
  G4OpticalSurface* photocath_opsurf = new G4OpticalSurface(
    "photocath_opsurf", glisur, polished, dielectric_metal);
  photocath_opsurf->SetMaterialPropertiesTable(photocath_mt);

  //**Create logical skin surfaces
  new G4LogicalSkinSurface("housing_surf", fHousing_log, OpScintHousingSurface);
  new G4LogicalSkinSurface("sphere_surface", fSphere_log, OpSphereSurface);
  new G4LogicalSkinSurface("photocath_surf", fPhotocath_log, photocath_opsurf);
}



G4VSolid* LXeMainVolume::BuildPlaSolid(G4double X,G4double Y,G4double Z,G4double XYPMT,G4double tpLength){
    G4double tolerance = 0.001*mm;
    G4Box* solidDet = new G4Box("solidDet",X/2.,Y/2.,Z/2.);
    G4Trd* trapedge1 = new G4Trd("trd1",0.*mm,(X/2-XYPMT/2.)*mm,Y/2.+tolerance,Y/2.+tolerance,tpLength/2+tolerance);
    G4Trd* trapedge2 = new G4Trd("trd2",(X/2-XYPMT/2.)*mm,0.*mm,Y/2.+tolerance,Y/2.+tolerance,tpLength/2+tolerance);
    G4Trd* trapedge3 = new G4Trd("trd3",X/2,X/2,0*mm,(Y-XYPMT)+tolerance,tpLength/2+tolerance);
    G4Trd* trapedge4 = new G4Trd("trd4",X/2,X/2,(Y-XYPMT)+tolerance,0*mm,tpLength/2+tolerance);
    G4SubtractionSolid* solidDetSub1 = new G4SubtractionSolid("solidDetSub1",solidDet,trapedge1,0,G4ThreeVector(X/2,0*mm,Z/2-tpLength/2));
    G4SubtractionSolid* solidDetSub2 = new G4SubtractionSolid("solidDetSub2",solidDetSub1,trapedge1,0,G4ThreeVector(-X/2,0*mm,Z/2-tpLength/2));
    G4SubtractionSolid* solidDetSub3 = new G4SubtractionSolid("solidDetSub3",solidDetSub2,trapedge2,0,G4ThreeVector(X/2,0*mm,-(Z/2-tpLength/2)));
    G4SubtractionSolid* solidDetSub4 = new G4SubtractionSolid("solidDetSub4",solidDetSub3,trapedge2,0,G4ThreeVector(-X/2,0*mm,-(Z/2-tpLength/2)));
    G4SubtractionSolid* solidDetSub5 = new G4SubtractionSolid("solidDetSub5",solidDetSub4,trapedge3,0,G4ThreeVector(0*mm,Y/2,Z/2-tpLength/2));
    G4SubtractionSolid* solidDetSub6 = new G4SubtractionSolid("solidDetSub6",solidDetSub5,trapedge4,0,G4ThreeVector(0*mm,Y/2,-(Z/2-tpLength/2)));
    return solidDetSub6;
}



G4LogicalVolume* LXeMainVolume::BuildPlaLogic(G4int detID)
{
    char tempname[50];
    sprintf(tempname,"_%d",detID);
    G4String identifier=(G4String) tempname;
    //
    // Define detector
    //
    G4double outerTal = 0.63*mm;//0.63*mm;
    G4double outerTtape = outerTal+0.18*mm;//0.18*mm;
    G4double detX = 60*mm;
    G4double detY = 25*mm;
    G4double detZ = 800*mm;
    G4double PMTXY  = 23*mm;
    G4double tpLength = 50*mm;

    G4double headX = 60*mm;
    G4double headY = (4+27+4)*mm;
    G4double headZ = 48*mm;
    G4double headYin = 27*mm;
    G4double headPMTXY = 24*mm;

    G4Box* solidDetContainer = new G4Box("solidDetContainer",(detX+outerTtape)/2,headY/2,detZ/2);
    G4LogicalVolume* logicDetContainer =
            new G4LogicalVolume(solidDetContainer,          //its solid
                                G4Material::GetMaterial("Air"),         //its material
                                "DetContainer"+identifier);        //its name
    logicDetContainer->SetVisAttributes(ColorMap["invi"]);


    G4Box* solidHeadOut = new G4Box("headOut",headX/2,headY/2,headZ/2);
    G4Box* solidHeadIn = new G4Box("solidHeadIn",headX/2+0.002*mm,headYin/2,headZ/2+0.002*mm);
    G4Box* solidHeadTop = new G4Box("solidHeadBot",headX/2+0.002*mm,1*mm,headZ/2+0.002*mm);
    G4SubtractionSolid* solidHead1 = new G4SubtractionSolid("solidHead1",solidHeadOut,solidHeadIn,0,G4ThreeVector(0,0,0));
    G4SubtractionSolid* solidHead = new G4SubtractionSolid("solidHead",solidHead1,solidHeadTop,0,G4ThreeVector(0,headY-1*mm,0));

    G4double trapWX = (headX/2-headPMTXY/2);
    G4double trapWX2 = trapWX*2;

    G4Trd* trapleft =new G4Trd("trapleft",0.00001*mm,trapWX2/2,headYin/2,headYin/2,headZ/2);
    G4double torr = 0.001*mm;
    G4Box* traplefts = new G4Box("traplefts",trapWX2/2,headYin/2+torr,headZ/2+torr);
    G4SubtractionSolid* trapl = new G4SubtractionSolid("trapl",trapleft,traplefts,0,G4ThreeVector(trapWX2/2,0));
    G4Trd* trapright =new G4Trd("trapright",0.00001*mm,trapWX2/2,headYin/2,headYin/2,headZ/2);
    G4Box* traprights = new G4Box("traprights",trapWX2/2,headYin/2+torr,headZ/2+torr);
    G4SubtractionSolid* trapr = new G4SubtractionSolid("trapr",trapright,traprights,0,G4ThreeVector(trapWX2/2,0,0));
    G4RotationMatrix* rotmat = new G4RotationMatrix;
    rotmat->set(0,0,0);
    rotmat->rotateZ(180*degree);
    G4UnionSolid* trap = new G4UnionSolid("trap",trapl,trapr,rotmat,G4ThreeVector(-trapWX2-headPMTXY,0,0));

    G4NistManager* nist = G4NistManager::Instance();
    G4Element* elC  = nist->FindOrBuildElement("C");
    G4Element* elH  = nist->FindOrBuildElement("H");
    G4Element* elO  = nist->FindOrBuildElement("O");
    G4Material* PLAmat= new G4Material( "PLA" ,1.26*g/cm3 ,3);//, kStateUndefined, temperature);
    PLAmat -> AddElement(elC,3);
    PLAmat -> AddElement(elH,4);
    PLAmat -> AddElement(elO,2);
    G4LogicalVolume* logicHead1 =
            new G4LogicalVolume(solidHead,          //its solid
                                PLAmat,         //its material
                                "Head1"+identifier);        //its name
    G4LogicalVolume* logicHead2 =
            new G4LogicalVolume(trap,          //its solid
                                PLAmat,         //its material
                                "Head2"+identifier);        //its name

    logicHead1->SetVisAttributes(ColorMap["Yellow"]);
    logicHead2->SetVisAttributes(ColorMap["Yellow"]);

    G4VSolid* solidDetOuterTape = BuildPlaSolid(detX+outerTtape,detY+outerTtape,detZ,PMTXY+outerTtape,tpLength);


    G4LogicalVolume* logicDetOuterTape =
            new G4LogicalVolume(solidDetOuterTape,          //its solid
                                G4Material::GetMaterial("TAPE"),         //its material
                                "OuterDet2"+identifier);        //its name

    logicDetOuterTape->SetVisAttributes(ColorMap["Magenta"]);
    G4VSolid* solidDetOuterAl = BuildPlaSolid(detX+outerTal,detY+outerTal,detZ,PMTXY+outerTal,tpLength);


    G4VSolid* solidDet = BuildPlaSolid(detX,detY,detZ,PMTXY,tpLength);


    fHousing_log = new G4LogicalVolume(solidDetOuterAl,          //its solid
                                G4Material::GetMaterial("Al"),         //its material
                                "housing_log"+identifier);        //its name
    fHousing_log->SetVisAttributes(ColorMap["Magenta"]);


    fScint_log =  new G4LogicalVolume(solidDet,          //its solid
                                G4Material::GetMaterial("EJ200"),         //its material
                                "scint_log"+identifier);        //its name


    fScint_log->SetVisAttributes(ColorMap["Red"]);

    new G4PVPlacement(0,
                      G4ThreeVector(0.*mm,0.*mm,0*mm),             //position
                      fHousing_log,            //its logical volume
                      "InactivePV"+identifier,             //its name
                      logicDetOuterTape,             //its mother  volume
                      false,                 //no boolean operation
                      0,                 //copy number
                      fCheckOverlaps);       // checking overlaps

    new G4PVPlacement(0,
                      G4ThreeVector(0.*mm,0.*mm,0*mm),             //position
                      fScint_log,            //its logical volume
                      "InactivePV"+identifier,             //its name
                      fHousing_log,             //its mother  volume
                      false,                 //no boolean operation
                      0,                 //copy number
                      fCheckOverlaps);       // checking overlaps

    new G4PVPlacement(0,
                      G4ThreeVector(0.*mm,0.*mm,0*mm),             //position
                      logicDetOuterTape,            //its logical volume
                      "InactivePV"+identifier,             //its name
                      logicDetContainer,             //its mother  volume
                      false,                 //no boolean operation
                      0,                 //copy number
                      fCheckOverlaps);       // checking overlaps
    new G4PVPlacement(0,
                      G4ThreeVector(0.*mm,0.*mm,detZ/2-headZ/2),             //position
                      logicHead1,            //its logical volume
                      "InactivePV"+identifier,             //its name
                      logicDetContainer,             //its mother  volume
                      false,                 //no boolean operation
                      0,                 //copy number
                      fCheckOverlaps);       // checking overlaps
    new G4PVPlacement(0,
                      G4ThreeVector(0.*mm,0.*mm,-(detZ/2-headZ/2)),             //position
                      logicHead1,            //its logical volume
                      "InactivePV"+identifier,             //its name
                      logicDetContainer,             //its mother  volume
                      false,                 //no boolean operation
                      1,                 //copy number
                      fCheckOverlaps);       // checking overlaps
    new G4PVPlacement(0,
                      G4ThreeVector(headX/2,0.*mm,detZ/2-headZ/2),             //position
                      logicHead2,            //its logical volume
                      "InactivePV"+identifier,             //its name
                      logicDetContainer,             //its mother  volume
                      false,                 //no boolean operation
                      0,                 //copy number
                      fCheckOverlaps);       // checking overlaps


    G4RotationMatrix* rotmat2 = new G4RotationMatrix;
    rotmat2->set(0,0,0);
    rotmat2->rotateY(180*degree);
    new G4PVPlacement(rotmat2,
                      G4ThreeVector(-headX/2,0.*mm,-(detZ/2-headZ/2)),             //position
                      logicHead2,            //its logical volume
                      "InactivePV"+identifier,             //its name
                      logicDetContainer,             //its mother  volume
                      false,                 //no boolean operation
                      1,                 //copy number
                      fCheckOverlaps);       // checking overlaps

    return logicDetContainer;
}
