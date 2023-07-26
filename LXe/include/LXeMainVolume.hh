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
/// \file optical/LXe/include/LXeMainVolume.hh
/// \brief Definition of the LXeMainVolume class
//
#ifndef LXeMainVolume_h
#define LXeMainVolume_h 1

#include "LXeDetectorConstruction.hh"

#include "G4PVPlacement.hh"

#include <string>
#include "G4VSolid.hh"
#include "G4VisAttributes.hh"


#include <map>
#include <string>


class G4Box;
class G4LogicalVolume;
class G4Sphere;
class G4Tubs;

class LXeMainVolume : public G4PVPlacement
{
 public:
  LXeMainVolume(G4RotationMatrix* pRot, const G4ThreeVector& tlate,
                G4LogicalVolume* pMotherLogical, G4bool pMany, G4int pCopyNo,
                LXeDetectorConstruction* c);

  G4LogicalVolume* GetLogPhotoCath() { return fPhotocath_log; }
  G4LogicalVolume* GetLogScint() { return fScint_log; }

  G4int GetNScint(){return nofTubes;}

  std::vector<G4ThreeVector> GetPmtPositions() { return fPmtPositions; }
  std::vector<G4ThreeVector> GetScintPositions() { return fScintPositions; }

 private:
  void VisAttributes();
  void SurfaceProperties();

  void PlacePMTs(G4LogicalVolume* pmt_Log, G4RotationMatrix* rot, G4double& a,
                 G4double& b, G4double da, G4double db, G4double amin,
                 G4double bmin, G4int na, G4int nb, G4double& x, G4double& y,
                 G4double& z, G4int& k);

  void CopyValues();

  LXeDetectorConstruction* fConstructor;

  G4double fExpHall_x;
  G4double fExpHall_y;
  G4double fExpHall_z;

  G4double fScint_x;
  G4double fScint_y;
  G4double fScint_z;
  G4double fD_mtl;
  G4int fNx;
  G4int fNy;
  G4int fNz;
  G4double fOuterRadius_pmt;
  G4bool fSphereOn;
  G4double fRefl;

  // Basic Volumes
  //
  G4Box* fExpHall_box;
//  G4Box* fScint_box;
//  G4Box* fHousing_box;
  G4Box* fPmt;
  G4Box* fPhotocath;

  G4Sphere* fSphere;

  // Logical volumes
  //

  G4int nofTubes;

  G4LogicalVolume* fExpHall_log;
  G4LogicalVolume* fScint_log;
  G4LogicalVolume* fHousing_log;
  G4LogicalVolume* fPmt_log;
  G4LogicalVolume* fPhotocath_log;
  G4LogicalVolume* fSphere_log;

  G4LogicalVolume* fDetContainer_log;

  std::map<std::string,G4VisAttributes*> ColorMap;
  G4bool fCheckOverlaps;

  G4LogicalVolume* BuildPlaLogic(G4int);
  G4VSolid* BuildPlaSolid(G4double X, G4double Y, G4double Z, G4double XYPMT, G4double tpLength);

  // Sensitive Detectors positions
  std::vector<G4ThreeVector> fPmtPositions;


  std::vector<G4ThreeVector> fScintPositions;


};

#endif