//---------------------------------------------------------------------
// Create the solids defining an Eurogam Phase-II detector
//---------------------------------------------------------------------
#include "GARi.hh"
//#include "MyMaterials.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Trap.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4LogicalVolume.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "globals.hh"

#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"

#include <stdio.h>
//using namespace std;
GARI::GARI()
{
    GARIsegZ = 2*mm;
    GARIbaseZ = 1*mm;
    GARIsegX = 1.3*mm;
    GARIsegY = 1.3*mm;
    segGap = 0.2*mm;
    fnGARIseg = 34;

    GARIBaseX = fnGARIseg*(GARIsegX+segGap);//51*mm;
    GARIBaseY = fnGARIseg*(GARIsegY+segGap);//51*mm;

    GARIBaseGap = 2*mm;

    GARIpmtX = 48.5*mm;
    GARIpmtY = 48.5*mm;

    height_pmt = 32.5 * mm;

    fOuterRadius_pmthousing =  51.7 * mm;


    fOuterRadius_glass =  51.0 * mm;
    height_photocathode = 0.0635 * cm;
    fOuterRadius_pmt = 48.5 * mm;

    GARIpmtZ = height_pmt;

    GlassWindowX = GARIBaseX*3+GARIBaseGap*2+10*mm;
    GlassWindowY = GARIBaseY*3+GARIBaseGap*2+10*mm;
    GlassWindowZ = 6*mm;

    STLwindowX = 6.*cm;
    STLwindowY = 11*cm;


    GARIbackTapeZ = 0.2*mm;

    GARIX = GlassWindowX+STLwindowX*2;
    GARIY = GlassWindowY+STLwindowY*2;
    GARIZ = GARIbaseZ+GARIsegZ+GlassWindowZ+GARIpmtZ+GARIbackTapeZ;


    STLwindowZ = (GARIZ/2-GetGlassCenter())*2;
    GARITapeoutZ = GARIZ-STLwindowZ-GARIbackTapeZ;

    OffSetSurfaceGlass = -(GARIZ/2-GARIbaseZ-GARIsegZ-GlassWindowZ) * cos(45*pi/180.);
//    OffSetSurfaceGlass = (GARIZ/2-GARIbaseZ-GARIsegZ) * cos(45*pi/180.);

    G4VisAttributes* Red=
            new G4VisAttributes(G4Colour(1.0,0.0,0.0,0.3));  // Red
    Red->SetVisibility(true);
    Red->SetForceSolid(true);

    G4VisAttributes* RedL=
            new G4VisAttributes(G4Colour(1.0,0.0,0.0,0.1));  // Red
    RedL->SetVisibility(true);
    RedL->SetForceSolid(true);

    G4VisAttributes* Green=
            new G4VisAttributes(G4Colour(0.0,1.0,0.0,0.3));  // Green
    Green->SetVisibility(true);
    Green->SetForceSolid(true);

    G4VisAttributes* Blue=
            new G4VisAttributes(G4Colour(0.0,0.0,1.0,0.3));  // Blue
    Blue->SetVisibility(true);
    Blue->SetForceSolid(true);

    G4VisAttributes* BlueL=
            new G4VisAttributes(G4Colour(0.0,0.0,1.0,0.1));  // Blue Light
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
    ColorMap.emplace("RedL",RedL);
    ColorMap.emplace("wireFrame",wireFrame);


    //create the solids.....
    CreateSolids();
    SetOpticalProperties();

}

//Destructor
GARI::~GARI() {
    //
}

void GARI::SetPosition(G4ThreeVector thisPos) {
    position = thisPos*mm;
}

void GARI::SetRotation(G4RotationMatrix thisRot) {
    rotation = thisRot;
}



//---------------------------------------------------------------------
// Create the solids defining Phase-II GARIs
//---------------------------------------------------------------------
void  GARI::CreateSolids()
{
    G4Box* solidDetContainer = new G4Box("solidDetContainer",
                                         GARIX/2,GARIY/2,GARIZ/2);
    logicGARI = new G4LogicalVolume(solidDetContainer,          //its solid
                                    G4Material::GetMaterial("Vacuum"),         //its material
                                    "GARIOut");        //its name
    logicGARI->SetVisAttributes(ColorMap["wireFrame"]);

    G4Box* solidGlassWindow = new G4Box("solidGlassWindow",
                                         GlassWindowX/2,GlassWindowY/2,GlassWindowZ/2);
    logicGARIGlass = new G4LogicalVolume(solidGlassWindow,          //its solid
                                    G4Material::GetMaterial("Borosilicate_Glass"),         //its material
                                    "logicGARIGlass");        //its name
    logicGARIGlass->SetVisAttributes(ColorMap["Green"]);


    G4Box* solidbase = new G4Box("solidbase",GARIBaseX/2,GARIBaseX/2,GARIbaseZ/2);
    logicGARIbase =
            new G4LogicalVolume(solidbase,          //its solid
                                G4Material::GetMaterial("EJ200"),         //its material
                                "GARIBase");        //its name
    logicGARIbase->SetVisAttributes(ColorMap["Yellow"]);
    G4Box* solidSeg = new G4Box("solidSeg",GARIsegX/2,GARIsegY/2,GARIsegZ/2);
    logicGARIseg =
            new G4LogicalVolume(solidSeg,          //its solid
                                G4Material::GetMaterial("EJ200"),         //its material
                                "logigGARIseg");        //its name
    logicGARIseg->SetVisAttributes(ColorMap["Red"]);

//    G4Box* solidPMT = new G4Box("solidpmt",GARIpmtX/2,GARIpmtY/2,GARIpmtZ/2);
//    logicGARIpmt =
//            new G4LogicalVolume(solidPMT,          //its solid
//                                G4Material::GetMaterial("Al"),         //its material
//                                "logigGARIpmt");        //its name
//    logicGARIpmt->SetVisAttributes(ColorMap["Grey"]);

    //****************** Build PMTs

    G4Box* fPmt = new G4Box("pmt_tube", fOuterRadius_pmthousing/2, fOuterRadius_pmthousing/2, height_pmt/2);
    G4Box* fPmt_glass = new G4Box("pmt_tube_glasswindow", fOuterRadius_glass/2, fOuterRadius_glass/2, height_photocathode);
    // the "photocathode" is a metal slab at the back of the glass that
    // is only a very rough approximation of the real thing since it only
    // absorbs or detects the photons based on the efficiency set below
    G4Box* fPhotocath = new G4Box("photocath_tube", fOuterRadius_pmt/2, fOuterRadius_pmt/2, height_photocathode/2);

    fPmt_log =  new G4LogicalVolume(fPmt, G4Material::GetMaterial("Al"), "pmt_log");

    G4LogicalVolume* Pmt_log_glass =
            new G4LogicalVolume(fPmt_glass, G4Material::GetMaterial("Borosilicate_Glass"), "pmt_log");
    fPhotocath_log = new G4LogicalVolume( fPhotocath, G4Material::GetMaterial("Bialkali"), "photocath_log");

    new G4PVPlacement(0, G4ThreeVector(0., 0., -height_photocathode/2), fPhotocath_log,
                      "photocath", Pmt_log_glass, false, 0);

    new G4PVPlacement(0, G4ThreeVector(0., 0., +height_pmt / 2.-height_photocathode), Pmt_log_glass,
                      "photocath_glass", fPmt_log, false, 0);

    Pmt_log_glass->SetVisAttributes(ColorMap["wireFrame"]);
//    fPmt_log->SetVisAttributes(ColorMap["Bluein"]);
//    fPhotocath_log->SetVisAttributes(ColorMap["Grey"]);
    fPmt_log->SetVisAttributes(ColorMap["invi"]);
    fPhotocath_log->SetVisAttributes(ColorMap["invi"]);
    //***********Arrange scintillator and Pmt**********


    //! outer window
    // STL outer window

    G4Box* solidOuterBox = new G4Box("solidOuterBox",
                                         GlassWindowX/2+STLwindowX,GlassWindowY/2+STLwindowY,STLwindowZ/2);
    G4Box* solidInnerBox = new G4Box("solidInnerBox",
                                         GlassWindowX/2,GlassWindowY/2,STLwindowZ/2+0.1*mm);
    G4SubtractionSolid* solidSTLWindow = new G4SubtractionSolid("solidSTLWindow",solidOuterBox,solidInnerBox,0,G4ThreeVector(0,0,0));

    logicSTLWindow = new G4LogicalVolume(solidSTLWindow,          //its solid
                                    G4Material::GetMaterial("STL_NIST"),         //its material
                                    "logicSTLWindow");        //its name
    logicSTLWindow->SetVisAttributes(ColorMap["invi"]);
//    logicSTLWindow->SetVisAttributes(ColorMap["Red"]);


    G4Box* solidTapeBack = new G4Box("solidTapeBack",
                                         GlassWindowX/2+STLwindowX,GlassWindowY/2+STLwindowY,GARIbackTapeZ/2);
    logicTapeBack = new G4LogicalVolume(solidTapeBack,          //its solid
                                    G4Material::GetMaterial("TAPE"),         //its material
                                    "logicTapeBack");        //its name
    logicTapeBack->SetVisAttributes(ColorMap["invi"]);
//    logicTapeBack->SetVisAttributes(ColorMap["RedL"]);


    G4Box* solidOuterBoxAround = new G4Box("solidOuterBoxAround",
                                         GlassWindowX/2+STLwindowX,GlassWindowY/2+STLwindowY,GARITapeoutZ/2);
    G4Box* solidInnerBoxAround = new G4Box("solidInnerBoxAround",
                                         GlassWindowX/2+STLwindowX-GARIbackTapeZ/2,GlassWindowY/2+STLwindowY-GARIbackTapeZ/2,GARITapeoutZ/2+0.1*mm);
    G4SubtractionSolid* solidTapeOut = new G4SubtractionSolid("solidTapeOut",solidOuterBoxAround,solidInnerBoxAround,0,G4ThreeVector(0,0,0));

    logicTapeOut = new G4LogicalVolume(solidTapeOut,          //its solid
                                    G4Material::GetMaterial("TAPE"),         //its material
                                    "logicTapeOut");        //its name
    logicTapeOut->SetVisAttributes(ColorMap["invi"]);
//    logicTapeOut->SetVisAttributes(ColorMap["RedL"]);


}

//------------------------------------------------------------------
void GARI::SetOpticalProperties()
{
    //**Photocathode surface properties
    std::vector<G4double> ephoton = { 2.484*eV, 2.615*eV, 2.760*eV, 2.922*eV, 3.105*eV };
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
}
//------------------------------------------------------------------
void GARI::Placement(G4int copyNo, G4VPhysicalVolume* physiMother)
{

    G4double stepXY[] = {-fOuterRadius_pmthousing-GARIBaseGap,0*mm,fOuterRadius_pmthousing+GARIBaseGap};
    G4int k = 0;
    for (G4int i = 0 ;i<3;i++){
        for (G4int j = 0 ;j<3;j++){
            //! placing base
            new G4PVPlacement(0,
                              G4ThreeVector(stepXY[i],stepXY[j],GARIZ/2-GARIbaseZ/2),             //position
                              logicGARIbase,            //its logical volume
                              "PhysicsGARIBase",             //its name
                              logicGARI,             //its mother  volume
                              false,                 //no boolean operation
                              fCopyNo+k,                 //copy number
                              fCheckOverlaps);       // checking overlap

            k++;
        }
    }

    k = 0;

    for (G4int i = 0 ;i<3;i++){
        for (G4int j = 0 ;j<3;j++){
//            if ((i>1)&&(j>1))
//                continue;
            //! placing seg
            G4double offsetX = stepXY[i]-GARIBaseX/2;
            G4double offsetY = stepXY[i]-GARIBaseY/2;
//            G4double offsetX = GARIBaseX/2;
//            G4double offsetY = GARIBaseY/2;
            G4double xx = offsetX+GARIsegX/2;
            for (G4int l = 0 ;l<fnGARIseg;l++){
                G4double yy = offsetY+GARIsegY/2;
                for (G4int n = 0 ;n<fnGARIseg;n++){
                    new G4PVPlacement(0,
                                      G4ThreeVector(xx,yy,GARIZ/2-GARIbaseZ-GARIsegZ/2),             //position
                                      logicGARIseg,            //its logical volume
                                      "PhysicsGARISeg",             //its name
                                      logicGARI,             //its mother  volume
                                      false,                 //no boolean operation
                                      fCopyNo+k,                 //copy number
                                      0);       // checking overlap
                    yy+=GARIsegY+segGap;
                }
                xx+=GARIsegX+segGap;
            }
            k++;
        }
    }


    //! placing glass
    new G4PVPlacement(0,
                      G4ThreeVector(0,0,GARIZ/2-GARIbaseZ-GARIsegZ-GlassWindowZ/2),             //position
                      logicGARIGlass,            //its logical volume
                      "PhysicsGARIBase",             //its name
                      logicGARI,             //its mother  volume
                      false,                 //no boolean operation
                      fCopyNo,                 //copy number
                      fCheckOverlaps);       // checking overlap

    //! STL structure outside of the glass
    new G4PVPlacement(0,
                      G4ThreeVector(0,0,GARIZ/2-GARIbaseZ-GARIsegZ-GlassWindowZ/2),             //position
                      logicSTLWindow,            //its logical volume
                      "physSTLWindow",             //its name
                      logicGARI,             //its mother  volume
                      false,                 //no boolean operation
                      0,                 //copy number
                      fCheckOverlaps);       // checking overlaps

    //! placing pmts
//    new G4PVPlacement(0,
//                      G4ThreeVector(0,0,GARIZ/2-GARIbaseZ-GARIsegZ-GlassWindowZ-GARIpmtZ/2),             //position
//                      fPmt_log,            //its logical volume
//                      "PhysicsGARIPMT",             //its name
//                      logicGARI,             //its mother  volume
//                      false,                 //no boolean operation
//                      fCopyNo,                 //copy number
//                      fCheckOverlaps);

    k = 0;
    for (G4int i = 0 ;i<3;i++){
        for (G4int j = 0 ;j<3;j++){
            //! placing base
            new G4PVPlacement(0,
                              G4ThreeVector(stepXY[i],stepXY[j],GARIZ/2-GARIbaseZ-GARIsegZ-GlassWindowZ-GARIpmtZ/2),             //position
                              fPmt_log,            //its logical volume
                              "PhysicsGARIPMT",             //its name
                              logicGARI,             //its mother  volume
                              false,                 //no boolean operation
                              fCopyNo+k,                 //copy number
                              fCheckOverlaps);
            k++;
        }
    }

    //! backing tape
    new G4PVPlacement(0,
                      G4ThreeVector(0,0,GARIZ/2-GARIbaseZ-GARIsegZ-GlassWindowZ-GARIpmtZ-GARIbackTapeZ/2),             //position
                      logicTapeBack,            //its logical volume
                      "physTapeBack",             //its name
                      logicGARI,             //its mother  volume
                      false,                 //no boolean operation
                      0,                 //copy number
                      fCheckOverlaps);       // checking overlaps

    //! sourring tape

    G4double STLPos = GARIZ/2-GARIbaseZ-GARIsegZ-GlassWindowZ/2;
    new G4PVPlacement(0,
                      G4ThreeVector(0,0,STLPos-STLwindowZ/2-GARITapeoutZ/2),             //position
                      logicTapeOut,            //its logical volume
                      "physSTLWindow",             //its name
                      logicGARI,             //its mother  volume
                      false,                 //no boolean operation
                      0,                 //copy number
                      fCheckOverlaps);       // checking overlaps


    //! Place whole things
    new G4PVPlacement(&rotation,
                      position,             //position
                      "PhysicsGARI",             //its name
                      logicGARI,            //its logical volume
                      physiMother,             //its mother  volume
                      false,                 //no boolean operation
                      copyNo,                 //copy number
                      fCheckOverlaps);       // checking overlaps

}
