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
/// \file optical/LXe/src/LXeHistoManager.cc
/// \brief Implementation of the LXeHistoManager class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



#include "LXeHistoManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeHistoManager::LXeHistoManager()
  : fFileName("lxe")
{
  Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeHistoManager::~LXeHistoManager() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeHistoManager::Book()
{
  // Create or get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("root");
  analysisManager->SetFileName(fFileName);
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(true);  // enable inactivation of histograms

  // Define histogram indices, titles

  // Default values (to be reset via /analysis/h1/set command)
  G4int nbins   = 100;
  G4double vmin = 0.;
  G4double vmax = 100.;

  // 0
  analysisManager->CreateH1("0", "dummy", nbins, vmin, vmax);
  // 1
  analysisManager->CreateH1("eE", "Electron energy",500, 0,
                            10000);
  // 2
  analysisManager->CreateH1("hits above threshold",
                            "hits per event above threshold", nbins, vmin,
                            vmax);
  // 3
  analysisManager->CreateH1("scintillation", "scintillation photons per event",
                            nbins, vmin, vmax);
  // 4
  analysisManager->CreateH1("Cerenkov", "Cerenkov photons per event", nbins,
                            vmin, vmax);
  // 5
  analysisManager->CreateH1("absorbed", "absorbed photons per event", nbins,
                            vmin, vmax);
  // 6
  analysisManager->CreateH1("boundary absorbed",
                            "photons absorbed at boundary per event", nbins,
                            vmin, vmax);
  // 7
  analysisManager->CreateH1(
    "E dep", "energy deposition in scintillator per event", nbins, vmin, vmax);

  // 8
  analysisManager->CreateH1(
    "pulse", "Time distribution", 200, 0, 1000);

  // 0
  analysisManager->CreateH2(
    "hxy0", "xy distribution", 200, -100, 100,200, -100, 100);
  analysisManager->CreateH2(
    "hxy1", "xy distribution", 200, -100, 100,200, -100, 100);
  analysisManager->CreateH2(
    "hxy2", "xy distribution", 200, -100, 100,200, -100, 100);
  analysisManager->CreateH2(
    "hxy3", "xy distribution", 200, -100, 100,200, -100, 100);

  analysisManager->CreateNtuple("tofu","tofu");
  //! neutron hits
  analysisManager->CreateNtupleIColumn("evt");//0
  analysisManager->CreateNtupleDColumn("primE");//1
  analysisManager->CreateNtupleIColumn("scint_mult");//2
  analysisManager->CreateNtupleIColumn("pmt_mult");//3

  for (G4int i=0;i<MAX_SCINT_MULT;i++){
      char tmpStr[500];
      sprintf(tmpStr,"sid%i",i);
      analysisManager->CreateNtupleIColumn(tmpStr);
      sprintf(tmpStr,"se%i",i);
      analysisManager->CreateNtupleDColumn(tmpStr);
      sprintf(tmpStr,"st%i",i);
      analysisManager->CreateNtupleDColumn(tmpStr);
      sprintf(tmpStr,"stmin%i",i);
      analysisManager->CreateNtupleDColumn(tmpStr);
      sprintf(tmpStr,"sm%i",i);
      analysisManager->CreateNtupleIColumn(tmpStr);
      sprintf(tmpStr,"sx%i",i);
      analysisManager->CreateNtupleDColumn(tmpStr);
      sprintf(tmpStr,"sy%i",i);
      analysisManager->CreateNtupleDColumn(tmpStr);
      sprintf(tmpStr,"sz%i",i);
      analysisManager->CreateNtupleDColumn(tmpStr);
  }

  for (G4int i=0;i<MAX_PMT_MULT;i++){
      char tmpStr[500];
      sprintf(tmpStr,"pid%i",i);
      analysisManager->CreateNtupleIColumn(tmpStr);
      sprintf(tmpStr,"pm%i",i);
      analysisManager->CreateNtupleIColumn(tmpStr);
      sprintf(tmpStr,"pt%i",i);
      analysisManager->CreateNtupleDColumn(tmpStr);
      sprintf(tmpStr,"ptmin%i",i);
      analysisManager->CreateNtupleDColumn(tmpStr);
      sprintf(tmpStr,"px%i",i);
      analysisManager->CreateNtupleDColumn(tmpStr);
      sprintf(tmpStr,"py%i",i);
      analysisManager->CreateNtupleDColumn(tmpStr);
      sprintf(tmpStr,"pz%i",i);
      analysisManager->CreateNtupleDColumn(tmpStr);
  }
  analysisManager->FinishNtuple();

  // Create all histograms as inactivated
  for(G4int i = 0; i < analysisManager->GetNofH1s(); ++i)
  {
    analysisManager->SetH1Activation(i, true);
  }
}
