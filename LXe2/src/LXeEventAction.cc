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
/// \file optical/LXe/src/LXeEventAction.cc
/// \brief Implementation of the LXeEventAction class
//
//
#include "LXeEventAction.hh"

#include "LXeDetectorConstruction.hh"
#include "LXeHistoManager.hh"
#include "LXePMTHit.hh"
#include "LXeRun.hh"
#include "LXeScintHit.hh"
#include "LXeTrajectory.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ios.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Trajectory.hh"
#include "G4TrajectoryContainer.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"

#include <fstream>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef struct
{
    G4int smult;
    G4int sid[MAX_SCINT_MULT];
    G4double se[MAX_SCINT_MULT];
    G4double st[MAX_SCINT_MULT];
    G4double stmin[MAX_SCINT_MULT];
    G4int sm[MAX_SCINT_MULT];
    G4double sx[MAX_SCINT_MULT];
    G4double sy[MAX_SCINT_MULT];
    G4double sz[MAX_SCINT_MULT];

    G4int pmult;
    G4int pid[MAX_PMT_MULT];
    G4int pm[MAX_PMT_MULT];
    G4double pt[MAX_PMT_MULT];
    G4double ptmin[MAX_PMT_MULT];
    G4double px[MAX_PMT_MULT];
    G4double py[MAX_PMT_MULT];
    G4double pz[MAX_PMT_MULT];
}evt_t;

void ClearEvt(evt_t* evt){
    evt->smult=0;
    memset(evt->sid,0,sizeof(evt->sid));
    memset(evt->se,0,sizeof(evt->se));
    memset(evt->st,0,sizeof(evt->st));
    memset(evt->stmin,0,sizeof(evt->stmin));
    memset(evt->sm,0,sizeof(evt->sm));
    memset(evt->sx,0,sizeof(evt->sx));
    memset(evt->sy,0,sizeof(evt->sy));
    memset(evt->sz,0,sizeof(evt->sz));

    evt->pmult=0;
    memset(evt->pid,0,sizeof(evt->pid));
    memset(evt->pm,0,sizeof(evt->pm));
    memset(evt->pt,0,sizeof(evt->pt));
    memset(evt->ptmin,0,sizeof(evt->ptmin));
    memset(evt->px,0,sizeof(evt->px));
    memset(evt->py,0,sizeof(evt->py));
    memset(evt->pz,0,sizeof(evt->pz));
}

void FillEvt(evt_t* evt, G4AnalysisManager* analysisManager)
{
    G4int iCol = 2;
    analysisManager->FillNtupleIColumn(iCol,evt->smult);
    iCol++;
    analysisManager->FillNtupleIColumn(iCol,evt->pmult);
    iCol++;
    for (G4int i=0;i<MAX_SCINT_MULT;i++){
        analysisManager->FillNtupleIColumn(iCol,evt->sid[i]);
        iCol++;
        analysisManager->FillNtupleDColumn(iCol,evt->se[i]);
        iCol++;
        analysisManager->FillNtupleDColumn(iCol,evt->st[i]);
        iCol++;
        analysisManager->FillNtupleDColumn(iCol,evt->stmin[i]);
        iCol++;
        analysisManager->FillNtupleIColumn(iCol,evt->sm[i]);
        iCol++;
        analysisManager->FillNtupleDColumn(iCol,evt->sx[i]);
        iCol++;
        analysisManager->FillNtupleDColumn(iCol,evt->sy[i]);
        iCol++;
        analysisManager->FillNtupleDColumn(iCol,evt->sz[i]);
        iCol++;
    }

    for (G4int i=0;i<MAX_PMT_MULT;i++){
        analysisManager->FillNtupleIColumn(iCol,evt->pid[i]);
        iCol++;
        analysisManager->FillNtupleIColumn(iCol,evt->pm[i]);
        iCol++;
        analysisManager->FillNtupleDColumn(iCol,evt->pt[i]);
        iCol++;
        analysisManager->FillNtupleDColumn(iCol,evt->ptmin[i]);
        iCol++;
        analysisManager->FillNtupleDColumn(iCol,evt->px[i]);
        iCol++;
        analysisManager->FillNtupleDColumn(iCol,evt->py[i]);
        iCol++;
        analysisManager->FillNtupleDColumn(iCol,evt->pz[i]);
        iCol++;
    }
}

LXeEventAction::LXeEventAction(const LXeDetectorConstruction* det)
    : fDetector(det)
    , fScintCollID(-1)
    , fPMTCollID(-1)
    , fVerbose(0)
    , fPMTThreshold(1)
    , fForcedrawphotons(false)
    , fForcenophotons(false)
{
    fEventMessenger = new LXeEventMessenger(this);

    fHitCount                = 0;
    fPhotonCount_Scint       = 0;
    fPhotonCount_Ceren       = 0;
    fAbsorptionCount         = 0;
    fBoundaryAbsorptionCount = 0;
    fTotE                    = 0.0;

    fConvPosSet = false;
    fEdepMax    = 0.0;

    fPMTsAboveThreshold = 0;
    fPrimeE = -1.;

    // get number of scintillator
    fNScint = 0;
    G4int tubeGroup[1000];
    G4double x[1000];
    G4double y[1000];
    G4double z[1000];
    G4double rot[1000];
    std::ifstream infile("detmap.txt");
    while (!infile.eof()){
        G4int myid;
        infile>>myid;
        infile>>x[myid]>>y[myid]>>z[myid]>>rot[myid]>>tubeGroup[myid];
        fNScint++;
    }
    infile.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeEventAction::~LXeEventAction() { delete fEventMessenger; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeEventAction::BeginOfEventAction(const G4Event* anEvent)
{
    fHitCount                = 0;
    fPhotonCount_Scint       = 0;
    fPhotonCount_Ceren       = 0;
    fAbsorptionCount         = 0;
    fBoundaryAbsorptionCount = 0;
    fTotE                    = 0.0;

    fPrimeE = -1.;

    fConvPosSet = false;
    fEdepMax    = 0.0;

    fPMTsAboveThreshold = 0;

    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    if(fScintCollID < 0)
        fScintCollID = SDman->GetCollectionID("scintCollection");
    if(fPMTCollID < 0)
        fPMTCollID = SDman->GetCollectionID("pmtHitCollection");
    //    G4cout<<"beg"<<anEvent->GetEventID()<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeEventAction::EndOfEventAction(const G4Event* anEvent)
{
    G4TrajectoryContainer* trajectoryContainer =
            anEvent->GetTrajectoryContainer();

    G4int n_trajectories = 0;
    if(trajectoryContainer)
        n_trajectories = trajectoryContainer->entries();

    // extract the trajectories and draw them
    if(G4VVisManager::GetConcreteInstance())
    {
        for(G4int i = 0; i < n_trajectories; ++i)
        {
            LXeTrajectory* trj =
                    (LXeTrajectory*) ((*(anEvent->GetTrajectoryContainer()))[i]);
            if(trj->GetParticleName() == "opticalphoton")
            {
                trj->SetForceDrawTrajectory(fForcedrawphotons);
                trj->SetForceNoDrawTrajectory(fForcenophotons);
            }
            trj->DrawTrajectory();
        }
    }

    auto analysisManager = G4AnalysisManager::Instance();

    LXeScintHitsCollection* scintHC = nullptr;
    LXePMTHitsCollection* pmtHC     = nullptr;
    G4HCofThisEvent* hitsCE         = anEvent->GetHCofThisEvent();



    analysisManager->FillH1(1,fPrimeE);

    // Get the hit collections
    if(hitsCE)
    {
        if(fScintCollID >= 0)
        {
            scintHC = (LXeScintHitsCollection*) (hitsCE->GetHC(fScintCollID));
        }
        if(fPMTCollID >= 0)
        {
            pmtHC = (LXePMTHitsCollection*) (hitsCE->GetHC(fPMTCollID));
        }
    }

    //    analysisManager->FillNtupleIColumn(0,anEvent->GetEventID());
    //    analysisManager->FillNtupleDColumn(1,anEvent->GetPrimaryVertex()->GetPrimary()->GetKineticEnergy());
    evt_t* evtdata = new evt_t;
    ClearEvt(evtdata);

    G4int mult[MAX_N_SCINT];
    G4double edep_Scint[MAX_N_SCINT];
    G4double t_Scint[MAX_N_SCINT];
    G4double tmin_Scint[MAX_N_SCINT];
    G4ThreeVector eWPos_Scint[MAX_N_SCINT];
    for (size_t i = 0; i < MAX_N_SCINT; ++i){
        eWPos_Scint->set(0.,0.,0.);
        edep_Scint[i] = 0;
        t_Scint[i] = 0;
        tmin_Scint[i] = 1000000.;
        mult[i] = 0;
    }

    if(scintHC)
    {
        size_t n_hit = scintHC->entries();
        G4ThreeVector eWeightPos(0.);
        G4int scintNo;
        G4ThreeVector scintPos(0.);
        G4double edep;
        G4double timeAve;
        G4double edepMax = 0;

        for(size_t i = 0; i < n_hit; ++i)
        {  // gather info on hits in scintillator
            scintNo = (*scintHC)[i]->GetScintNumber();
            scintPos = (*scintHC)[i]->GetScintPos();
            edep = (*scintHC)[i]->GetEdep();
            timeAve = (*scintHC)[i]->GetT();
            edep_Scint[scintNo] += edep;
            t_Scint[scintNo] += timeAve;

            if (timeAve<tmin_Scint[scintNo])
                tmin_Scint[scintNo] = timeAve;

            eWPos_Scint[scintNo] += (*scintHC)[i]->GetPos() * edep;
            mult[scintNo]++;

            fTotE += edep;

            eWeightPos +=
                    (*scintHC)[i]->GetPos() * edep;  // calculate energy weighted pos

            if(edep > edepMax)
            {
                edepMax = edep;  // store max energy deposit
                G4ThreeVector posMax = (*scintHC)[i]->GetPos();
                fPosMax              = posMax;
                fEdepMax             = edep;
            }
        }
        analysisManager->FillH1(7, fTotE);
        if(fTotE == 0.)
        {
            if(fVerbose > 0)
                G4cout << "No hits in the scintillator this event." << G4endl;
        }
        else
        {
            // Finish calculation of energy weighted position
            eWeightPos /= fTotE;
            fEWeightPos = eWeightPos;
            if(fVerbose > 0)
            {
                G4cout << "\tEnergy weighted position of hits in LXe : "
                       << eWeightPos / mm << G4endl;
            }
        }
        if(fVerbose > 0)
        {
            G4cout << "\tTotal energy deposition in scintillator : " << fTotE / keV
                   << " (keV)" << G4endl;
        }
    }


    G4int scintMult = 0;
    for (G4int i = 0; i < fNScint; ++i){
        if (edep_Scint[i]>0){
            if (scintMult<MAX_SCINT_MULT){
                evtdata->sid[scintMult] = i;
                evtdata->se[scintMult] = edep_Scint[i];
                evtdata->st[scintMult] = t_Scint[i]/mult[i];
                evtdata->stmin[scintMult] = tmin_Scint[i];
                evtdata->sm[scintMult] = mult[i];
                if (edep_Scint[i]>0)
                    eWPos_Scint[i]/=edep_Scint[i];
                evtdata->sx[scintMult] = eWPos_Scint[i].getX();
                evtdata->sy[scintMult] = eWPos_Scint[i].getY();
                evtdata->sz[scintMult] = eWPos_Scint[i].getZ();
            }
            scintMult++;
        }
    }
    evtdata->smult = scintMult;


    if(pmtHC)
    {
        G4ThreeVector reconPos(0., 0., 0.);

        G4int pmtNo;
        G4ThreeVector pmtPost(0., 0., 0.);

        size_t pmts = pmtHC->entries();
        evtdata->pmult = pmts;
        // Gather info from all PMTs
        for(size_t i = 0; i < pmts; ++i)
        {
            pmtNo = (*pmtHC)[i]->GetPMTNumber();
            pmtPost = (*pmtHC)[i]->GetPMTPos();
            fHitCount += (*pmtHC)[i]->GetPhotonCount();
            reconPos += (*pmtHC)[i]->GetPMTPos() * (*pmtHC)[i]->GetPhotonCount();
            if((*pmtHC)[i]->GetPhotonCount() >= fPMTThreshold)
            {
                ++fPMTsAboveThreshold;
            }
            else
            {  // wasn't above the threshold, turn it back off
                (*pmtHC)[i]->SetDrawit(false);
            }
            if (i<MAX_PMT_MULT){
                evtdata->pid[i] = pmtNo;
                evtdata->pm[i] = (*pmtHC)[i]->GetPhotonCount();
                evtdata->pt[i] = (*pmtHC)[i]->GetT();
                evtdata->ptmin[i] = (*pmtHC)[i]->GetTmin();
                evtdata->px[i] = (*pmtHC)[i]->GetPos().getX();
                evtdata->py[i] = (*pmtHC)[i]->GetPos().getY();
                evtdata->pz[i] = (*pmtHC)[i]->GetPos().getZ();
            }
            auto hitT = (*pmtHC)[i]->GetTi();
            for (auto it=hitT.begin();it!=hitT.end();it++){
                analysisManager->FillH1(8, *it);
            }
        }


        //        analysisManager->FillH1(1, fHitCount);
        analysisManager->FillH1(2, fPMTsAboveThreshold);

        if(fHitCount > 0)
        {  // don't bother unless there were hits
            reconPos /= fHitCount;
            if(fVerbose > 0)
            {
                G4cout << "\tReconstructed position of hits in LXe : " << reconPos / mm
                       << G4endl;
            }
            fReconPos = reconPos;
        }
        pmtHC->DrawAllHits();
    }

    FillEvt(evtdata,analysisManager);
    analysisManager->AddNtupleRow();


    delete evtdata;


    analysisManager->FillH1(3, fPhotonCount_Scint);
    analysisManager->FillH1(4, fPhotonCount_Ceren);
    analysisManager->FillH1(5, fAbsorptionCount);
    analysisManager->FillH1(6, fBoundaryAbsorptionCount);

    if(fVerbose > 0)
    {
        // End of event output. later to be controlled by a verbose level
        G4cout << "\tNumber of photons that hit PMTs in this event : " << fHitCount
               << G4endl;
        G4cout << "\tNumber of PMTs above threshold(" << fPMTThreshold
               << ") : " << fPMTsAboveThreshold << G4endl;
        G4cout << "\tNumber of photons produced by scintillation in this event : "
               << fPhotonCount_Scint << G4endl;
        G4cout << "\tNumber of photons produced by cerenkov in this event : "
               << fPhotonCount_Ceren << G4endl;
        G4cout << "\tNumber of photons absorbed (OpAbsorption) in this event : "
               << fAbsorptionCount << G4endl;
        G4cout << "\tNumber of photons absorbed at boundaries (OpBoundary) in "
               << "this event : " << fBoundaryAbsorptionCount << G4endl;
        G4cout << "Unaccounted for photons in this event : "
               << (fPhotonCount_Scint + fPhotonCount_Ceren - fAbsorptionCount -
                   fHitCount - fBoundaryAbsorptionCount)
               << G4endl;
    }

    // update the run statistics
    LXeRun* run = static_cast<LXeRun*>(
                G4RunManager::GetRunManager()->GetNonConstCurrentRun());

    run->IncHitCount(fHitCount);
    run->IncPhotonCount_Scint(fPhotonCount_Scint);
    run->IncPhotonCount_Ceren(fPhotonCount_Ceren);
    run->IncEDep(fTotE);
    run->IncAbsorption(fAbsorptionCount);
    run->IncBoundaryAbsorption(fBoundaryAbsorptionCount);
    run->IncHitsAboveThreshold(fPMTsAboveThreshold);


    //    G4cout<<"end"<<anEvent->GetEventID()<<G4endl;
    // If we have set the flag to save 'special' events, save here
    //  if(fPhotonCount_Scint + fPhotonCount_Ceren < fDetector->GetSaveThreshold())
    //  {
    //    G4RunManager::GetRunManager()->rndmSaveThisEvent();
    //  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
