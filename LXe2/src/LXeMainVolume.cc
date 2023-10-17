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
#include "G4Cons.hh"
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
            new G4VisAttributes(G4Colour(0.0,0.0,1.0,0.5));  // Blue
    Blue->SetVisibility(true);
    Blue->SetForceSolid(true);

    G4VisAttributes* BlueL=
            new G4VisAttributes(G4Colour(0.0,0.0,1.0,1.));  // Blue
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
            new G4VisAttributes(G4Colour(0.91,0.91,0.91,1.));  // Gray
    Grey->SetVisibility(true);
    Grey->SetForceSolid(true);

    G4VisAttributes* Grey2=
            new G4VisAttributes(G4Colour(0.,0.,0.,1.));  // Gray2
    Grey2->SetVisibility(true);
    Grey2->SetForceSolid(true);



    G4VisAttributes* wireFrame=
            new G4VisAttributes(G4Colour(1.0,0.0,0.0));  // Gray
    wireFrame->SetVisibility(true);
    wireFrame->SetForceSolid(false);
    wireFrame->SetForceWireframe(true);


    ColorMap.emplace("Red",Red);
    ColorMap.emplace("Green",Green);
    ColorMap.emplace("Blue",Blue);
    ColorMap.emplace("Yellow",Yellow);
    ColorMap.emplace("Turquoise",Turquoise);
    ColorMap.emplace("Magenta",Magenta);
    ColorMap.emplace("Grey",Grey);
    ColorMap.emplace("Grey2",Grey2);
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

    G4double fOuterRadius_pmthousing =  fOuterRadius_pmt;
    G4double fOuterRadius_glass =  fOuterRadius_pmthousing - 0.1*mm;
    G4double height_photocathode = fD_mtl;

    fPmt = new G4Box("pmt_tube", fOuterRadius_pmthousing/2, fOuterRadius_pmthousing/4, height_pmt/2);
    G4Box* fPmt_glass = new G4Box("pmt_tube_glasswindow", fOuterRadius_glass/2, fOuterRadius_glass/4, height_photocathode);
    // the "photocathode" is a metal slab at the back of the glass that
    // is only a very rough approximation of the real thing since it only
    // absorbs or detects the photons based on the efficiency set below
    fPhotocath = new G4Box("photocath_tube", (fOuterRadius_glass-0.1*mm)/2, (fOuterRadius_glass-0.1*mm)/4, height_photocathode/2);

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
    fPmt_log->SetVisAttributes(ColorMap["wireFrame"]);
    fPhotocath_log->SetVisAttributes(ColorMap["wireFrame"]);
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
    }
    nofTubes--;
    G4cout<<"Read "<<nofTubes<<G4endl;
    infile.close();

    G4VPhysicalVolume* physScint[1000];
    for (G4int i=0;i<nofTubes;i++)
    {
        //--- Define Detector Group --
        //--------------------------------------------------------------
        G4RotationMatrix* rotmat;
        rotmat = new G4RotationMatrix;
        rotmat->set(0,0,0);
        rotmat->rotateZ(rot[i]*degree);

#ifndef NOTOFU
        physScint[i]=new G4PVPlacement(rotmat,
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
    //! Reflectors
    G4int nRefl = 0;
    G4double refZ[1000];
    ifstream infile2("refzpos.txt");
    G4cout<<"reading refzpos.txt"<<G4endl;
    while (!infile2.eof()){
        G4int myid;
        infile2>>myid;
        infile2>>refZ[myid];
        nRefl++;
    }
    nRefl--;
    G4cout<<"Read reflector "<<nRefl<<G4endl;
    infile2.close();

    G4double al_x = 96*mm;
    G4double al_y = 96*mm;
    G4double al_z = 0.01*mm;

    G4Box* alwrapSolid = new G4Box("alwrapSolid",al_x/2.,al_y/2.,al_z/2.);
    G4LogicalVolume* alwrapLog = new G4LogicalVolume(alwrapSolid,          //its solid
                                G4Material::GetMaterial("Al"),         //its material
                                "alwrapLogic");        //its name
    alwrapLog->SetVisAttributes(ColorMap["wrireFrame"]);

    G4VPhysicalVolume* physAlwrap[100];
    for (G4int i=0;i<nRefl;i++)
    {
        //--- Define Detector Group --
        //--------------------------------------------------------------
        G4RotationMatrix* rotmat;
        rotmat = new G4RotationMatrix;
        rotmat->set(0,0,0);
        physAlwrap[i] =
            new G4PVPlacement(rotmat,
                              G4ThreeVector(0,0,refZ[i]*mm),             //position
                                       alwrapLog,            //its logical volume
                                       "PhysicsAlWrap",             //its name
                                       fExpHall_log,             //its mother  volume
                                       false,                 //no boolean operation
                                       i,                 //copy number
                                       fCheckOverlaps);       // checking overlaps
    }


    //! PMT
    G4double detxy = 96*mm;
    G4double xpos = -detxy/2-height_pmt/2;
    G4RotationMatrix* rotmatpmt;
    rotmatpmt = new G4RotationMatrix;
    rotmatpmt->set(0,0,0);
    rotmatpmt->rotateZ(90*degree);
    rotmatpmt->rotateX(90*degree);
    new G4PVPlacement(rotmatpmt, G4ThreeVector(xpos, 0., 0.), fPmt_log, "pmt0phys",
                      fExpHall_log, false, 0,fCheckOverlaps);
    fPmtPositions.push_back(G4ThreeVector(xpos, 0., 0.));

    xpos = detxy/2+height_pmt/2;
    G4RotationMatrix* rotmatpmt1 = new G4RotationMatrix;
    rotmatpmt1->set(0,0,0);
    rotmatpmt1->rotateZ(90*degree);
    rotmatpmt1->rotateX(90*degree+180*degree);
    new G4PVPlacement(rotmatpmt1, G4ThreeVector(xpos, 0., 0.), fPmt_log, "pmt1phys",
                      fExpHall_log, false, 1,fCheckOverlaps);
    fPmtPositions.push_back(G4ThreeVector(xpos, 0., 0.));


    xpos = -detxy/2-height_pmt/2;
    G4RotationMatrix* rotmatpmt2 = new G4RotationMatrix;
    rotmatpmt2->set(0,0,0);
    rotmatpmt2->rotateZ(90*degree);
    rotmatpmt2->rotateX(90*degree);
    rotmatpmt2->rotateY(90*degree);
    new G4PVPlacement(rotmatpmt2, G4ThreeVector(0., xpos, 0.), fPmt_log, "pmt2phys",
                      fExpHall_log, false, 2,fCheckOverlaps);
    fPmtPositions.push_back(G4ThreeVector(0., xpos, 0.));

    xpos = detxy/2+height_pmt/2;
    G4RotationMatrix* rotmatpmt3 = new G4RotationMatrix;
    rotmatpmt3->set(0,0,0);
    rotmatpmt3->rotateZ(90*degree);
    rotmatpmt3->rotateX(90*degree);
    rotmatpmt3->rotateY(90*degree+180*degree);
    new G4PVPlacement(rotmatpmt3, G4ThreeVector(0., xpos, 0.), fPmt_log, "pmt3phys",
                      fExpHall_log, false, 3,fCheckOverlaps);
    fPmtPositions.push_back(G4ThreeVector(0., xpos, 0.));

    VisAttributes();
    SetLogicalVolume(fExpHall_log);

    //! al reflector
    std::vector<G4double> ephoton = { 2.004*eV, 2.615*eV, 2.760*eV, 2.922*eV, 3.191*eV };
    std::vector<G4double> reflectivity_al     = { fRefl, fRefl,fRefl, fRefl,fRefl };
    std::vector<G4double> efficiency_al       = { 0.0, 0.0,0.0, 0.0,0.0 };
    G4MaterialPropertiesTable* scintHsngPT = new G4MaterialPropertiesTable();
    scintHsngPT->AddProperty("REFLECTIVITY", ephoton, reflectivity_al);
    scintHsngPT->AddProperty("EFFICIENCY", ephoton, efficiency_al);
    G4OpticalSurface* OpScintHousingSurface =
      new G4OpticalSurface("HousingSurface", unified, polished, dielectric_metal);
    OpScintHousingSurface->SetMaterialPropertiesTable(scintHsngPT);
    new G4LogicalSkinSurface("housing_surf", alwrapLog, OpScintHousingSurface);

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
    new G4LogicalSkinSurface("photocath_surf", fPhotocath_log, photocath_opsurf);

    //! ej200 plastics
    G4OpticalSurface* OpWaterSurface = new G4OpticalSurface("WaterSurface");
    OpWaterSurface->SetType(dielectric_LUTDAVIS);
    OpWaterSurface->SetFinish(Polished_LUT);
    OpWaterSurface->SetModel(DAVIS);

    for (G4int i=0;i<nofTubes;i++)
    {
        G4LogicalBorderSurface* WaterSurface =
                new G4LogicalBorderSurface("WaterSurface",
                                           physScint[i],this,OpWaterSurface);
    }
    //OpticalWaterSurface
    const G4int num = 2;
    G4double Ephoton[num] = {2.038*eV, 4.144*eV};
    G4double RefractiveIndex[num] = {1.35, 1.40};
    G4double SpecularLobe[num]    = {0.3, 0.3};
    G4double SpecularSpike[num]   = {0.2, 0.2};
    G4double Backscatter[num]     = {0.2, 0.2};

    G4double reflectivity[num]     = {0.2, 0.2};
    G4double transmittance[num]     = {0.8, 0.8};
    G4double efficiency[num]     = {0., 0.};


    G4MaterialPropertiesTable* myST1 = new G4MaterialPropertiesTable();

    myST1->AddProperty("RINDEX",                Ephoton, RefractiveIndex, num);
    myST1->AddProperty("SPECULARLOBECONSTANT",  Ephoton, SpecularLobe,    num);
    myST1->AddProperty("SPECULARSPIKECONSTANT", Ephoton, SpecularSpike,   num);
    myST1->AddProperty("BACKSCATTERCONSTANT",   Ephoton, Backscatter,     num);
    myST1->AddProperty("REFLECTIVITY",   Ephoton, reflectivity,     num);
    myST1->AddProperty("TRANSMITTANCE",   Ephoton, transmittance,     num);
    myST1->AddProperty("EFFICIENCY",   Ephoton, efficiency,     num);

    OpWaterSurface->SetMaterialPropertiesTable(myST1);



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

void LXeMainVolume::VisAttributes()
{
    G4VisAttributes* sphere_va = new G4VisAttributes();
    sphere_va->SetForceSolid(true);
    fSphere_log->SetVisAttributes(sphere_va);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



G4VSolid* LXeMainVolume::BuildPlaSolid(G4double X,G4double Y,G4double Z){
    return new G4Box("solidDet",X/2.,Y/2.,Z/2.);
}

G4LogicalVolume* LXeMainVolume::BuildPlaLogic(G4int detID)
{
    char tempname[50];
    sprintf(tempname,"_%d",detID);
    G4String identifier=(G4String) tempname;

    G4VSolid* solidDet = BuildPlaSolid(fScint_x,fScint_y,fScint_z);
    fScint_log =  new G4LogicalVolume(solidDet,          //its solid
                                      G4Material::GetMaterial("EJ200"),      //its material
                                      "scint_log"+identifier);        //its name

    fScint_log->SetVisAttributes(ColorMap["Blue"]);

    return fScint_log;
}
