#ifndef GARI_H
#define GARI_H 1

#include "globals.hh"
#include "G4Material.hh"
#include "G4VUserDetectorConstruction.hh"

#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include <map>
//class G4Box;
//class G4Tubs;
//class G4Polycone;
class G4UnionSolid;
class G4SubtractionSolid;
//class G4IntersectionSolid;
//class G4Polyhedra;
class G4LogicalVolume;
class G4VPhysicalVolume;
class MyDetectorMessenger;
//class MyMaterials;

class GARI {

public:
  GARI();
  GARI(G4String);
  ~GARI();

public:
  void SetPosition( G4ThreeVector );
  void SetRotation( G4RotationMatrix );
  void Placement(G4int, G4VPhysicalVolume*);

  G4double GetSizeX(){return GARIX;};
  G4double GetSizeY(){return GARIY;};
  G4double GetSizeZ(){return GARIZ;};
  G4double GetGlassCenter(){return GARIZ/2-GARIbaseZ-GARIsegZ-GlassWindowZ/2;};

  G4double GetOffsetSurfaceGlass(){return OffSetSurfaceGlass;};

  G4double GetGlassMaxRadius(){return sqrt(GlassWindowX/2*GlassWindowX/2+GlassWindowY/2*GlassWindowY/2);};

  void SetCopyNo(G4int copyNo){
      fCopyNo = copyNo;
  };

  G4int GetCopyNo(){
      return fCopyNo;
  };

   G4LogicalVolume* GetLogPhotoCath() { return fPhotocath_log; }

private:


  G4ThreeVector        position;
  G4RotationMatrix     rotation;


  G4double GARIX;
  G4double GARIY;

  G4double GARIBaseX;
  G4double GARIBaseY;
  G4double GARIsegZ;
  G4double GARIbaseZ;
  G4double GARIBaseGap;
  G4double GARIpmtZ;
  G4double GARIZ;
  G4double GARIsegX;
  G4double GARIsegY;
  G4double GARIpmtX;
  G4double GARIpmtY;
  G4double segGap;
  G4int fnGARIseg;

  G4double GlassWindowX;
  G4double GlassWindowY;
  G4double GlassWindowZ;

  G4double STLwindowX;
  G4double STLwindowY;
  G4double STLwindowZ;

  G4double GARIbackTapeZ;

  G4double GARITapeoutZ;

  G4double height_pmt;

  G4double fOuterRadius_pmthousing;
  G4double fOuterRadius_glass;

  G4double height_photocathode;
  G4double fOuterRadius_pmt;

  G4double OffSetSurfaceGlass;

  G4LogicalVolume* logicGARI;
  G4LogicalVolume* logicGARIbase;
  G4LogicalVolume* logicGARIseg;
  G4LogicalVolume* logicSTLWindow;
  G4LogicalVolume* logicTapeBack;
  G4LogicalVolume* logicTapeOut;



  G4LogicalVolume* logicGARIGlass;

//  G4LogicalVolume* logicGARIpmt;

  G4LogicalVolume* fPmt_log;
  G4LogicalVolume* fPhotocath_log;


  G4int fCheckOverlaps;
  G4int fCopyNo;

  std::map<std::string,G4VisAttributes*> ColorMap;



  void CreateSolids();
  void SetOpticalProperties();
//  void MakeMaterials();

};

#endif
