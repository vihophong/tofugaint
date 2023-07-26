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
/// \file optical/LXe/src/LXeScintSD.cc
/// \brief Implementation of the LXeScintSD class
//
//
#include "LXeScintSD.hh"

#include "LXeScintHit.hh"

#include "G4ios.hh"
#include "G4LogicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"
#include "G4VTouchable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeScintSD::LXeScintSD(G4String name)
  : G4VSensitiveDetector(name)
  , fScintPositionsX(nullptr)
  , fScintPositionsY(nullptr)
  , fScintPositionsZ(nullptr)
  , fHitsCID(-1)
{
  fScintCollection = nullptr;
  collectionName.insert("scintCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeScintSD::~LXeScintSD()
{
    delete fScintPositionsX;
    delete fScintPositionsY;
    delete fScintPositionsZ;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeScintSD::SetScintPositions(const std::vector<G4ThreeVector>& positions)
{
  for(size_t i = 0; i < positions.size(); ++i)
  {
    if(fScintPositionsX)
      fScintPositionsX->push_back(positions[i].x());
    if(fScintPositionsY)
      fScintPositionsY->push_back(positions[i].y());
    if(fScintPositionsZ)
      fScintPositionsZ->push_back(positions[i].z());
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeScintSD::Initialize(G4HCofThisEvent* hitsCE)
{
  fScintCollection =
    new LXeScintHitsCollection(SensitiveDetectorName, collectionName[0]);

  if(fHitsCID < 0)
  {
    fHitsCID = G4SDManager::GetSDMpointer()->GetCollectionID(fScintCollection);
  }
  hitsCE->AddHitsCollection(fHitsCID, fScintCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool LXeScintSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  if(edep == 0.)
    return false;  // No edep so don't count as hit

  G4StepPoint* thePrePoint = aStep->GetPreStepPoint();
  G4TouchableHistory* theTouchable =
    (G4TouchableHistory*) (aStep->GetPreStepPoint()->GetTouchable());
  G4VPhysicalVolume* thePrePV = theTouchable->GetVolume();

  G4StepPoint* thePostPoint = aStep->GetPostStepPoint();

  G4int scintNumber =
    thePostPoint->GetTouchable()->GetReplicaNumber(3);

  G4double timeAve= thePrePoint->GetGlobalTime()/2+thePostPoint->GetGlobalTime()/2;

  G4double timePre= thePrePoint->GetGlobalTime();

  // Get the average position of the hit
  G4ThreeVector pos = thePrePoint->GetPosition() + thePostPoint->GetPosition();
  pos /= 2.;


  LXeScintHit* scintHit = new LXeScintHit(thePrePV);
  scintHit->SetEdep(edep);
  scintHit->SetT(timeAve);
  scintHit->SetTpre(timePre);
  scintHit->SetPos(pos);
  scintHit->SetScintNumber(scintNumber);
  scintHit->SetScintPos((*fScintPositionsX)[scintNumber], (*fScintPositionsY)[scintNumber],
                        (*fScintPositionsZ)[scintNumber]);
  fScintCollection->insert(scintHit);

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
