#include "TChain.h"
#include "TLatex.h"
#include "TTree.h"
#include "TFile.h"

#include "TCanvas.h"
#include "TLegend.h"
#include "TPad.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom.h"
#include "TString.h"
#include "TLine.h"
#include "TSpectrum.h"
#include "TMarker.h"
#include "TString.h"
#include "TF1.h"
#include "TGraph.h"
#include <fstream>
#include <iostream>
#define MAX_SCINT_MULT 10
#define MAX_PMT_MULT 20

#define MAX_N_SCINT 150

#define INNER_RING_ID 36

#define MAX_N_HISTO 100


#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif



using namespace std;


TFile* outfile;
TTree* treeData[100];


void getH1d(char* infile,char* treename,Int_t treeid)
{
    TTree* tofu = 0;
    TFile* f = new TFile(infile);
    f->GetObject("tofu",tofu);

    //Declaration of leaves types
    Int_t           evt;
    Double_t        primE;
    Int_t           scint_mult;
    Int_t           pmt_mult;

    Int_t           sid[MAX_SCINT_MULT];
    Double_t        se[MAX_SCINT_MULT];
    Double_t        st[MAX_SCINT_MULT];
    Double_t        stmin[MAX_SCINT_MULT];
    Int_t           sm[MAX_SCINT_MULT];
    Double_t        sx[MAX_SCINT_MULT];
    Double_t        sy[MAX_SCINT_MULT];
    Double_t        sz[MAX_SCINT_MULT];

    Int_t           pid[MAX_PMT_MULT];
    Int_t           pm[MAX_PMT_MULT];
    Double_t        pt[MAX_PMT_MULT];
    Double_t        ptmin[MAX_PMT_MULT];
    Double_t        px[MAX_PMT_MULT];
    Double_t        py[MAX_PMT_MULT];
    Double_t        pz[MAX_PMT_MULT];


    // Set branch addresses.
    tofu->SetBranchAddress("evt",&evt);
    tofu->SetBranchAddress("primE",&primE);
    tofu->SetBranchAddress("scint_mult",&scint_mult);
    tofu->SetBranchAddress("pmt_mult",&pmt_mult);
    for (Int_t i=0;i<MAX_SCINT_MULT;i++){
        tofu->SetBranchAddress(Form("sid%d",i),&sid[i]);
        tofu->SetBranchAddress(Form("se%d",i),&se[i]);
        tofu->SetBranchAddress(Form("st%d",i),&st[i]);
        tofu->SetBranchAddress(Form("stmin%d",i),&stmin[i]);
        tofu->SetBranchAddress(Form("sm%d",i),&sm[i]);
        tofu->SetBranchAddress(Form("sx%d",i),&sx[i]);
        tofu->SetBranchAddress(Form("sy%d",i),&sy[i]);
        tofu->SetBranchAddress(Form("sz%d",i),&sz[i]);
    }

    for (Int_t i=0;i<MAX_PMT_MULT;i++){
        tofu->SetBranchAddress(Form("pid%d",i),&pid[i]);
        tofu->SetBranchAddress(Form("pm%d",i),&pm[i]);
        tofu->SetBranchAddress(Form("pt%d",i),&pt[i]);
        tofu->SetBranchAddress(Form("ptmin%d",i),&ptmin[i]);
        tofu->SetBranchAddress(Form("px%d",i),&px[i]);
        tofu->SetBranchAddress(Form("py%d",i),&py[i]);
        tofu->SetBranchAddress(Form("pz%d",i),&pz[i]);
    }

    Long64_t nentries = tofu->GetEntries();
    Long64_t nbytes = 0;



    outfile->cd();
    treeData[treeid] = new TTree(treename,treename);
    Double_t e[4];
    treeData[treeid]->Branch("e[4]",e,"e[4]/D");

    Int_t tmp[MAX_N_SCINT*2];
    for (int i=0;i<MAX_N_SCINT*2;i++) tmp[i]=-1;
    for (Long64_t ientry=0; ientry<nentries;ientry++) {
        nbytes += tofu->GetEntry(ientry);
//        if (scint_mult!=1) continue;
        for (Int_t j = 0 ; j < 4; j++) e[j]=-9999;
        for (int i=0;i<pmt_mult;i++){
            if (pid[i]<4){
                e[pid[i]] = pm[i];
                treeData[treeid]->Fill();
            }
        }
    }
}


void getDataVaryE()
{
    outfile = new TFile("results/calibdata.root","recreate");
    Double_t estep = 200;
    Double_t emin = 200;
    Double_t emax = 2000;
    Int_t ntrees = 0;

    while (true){
        if (emin>emax)
            break;
        getH1d(Form("build-LXe2/g_vary_e/LXe_gamma_%d.root",(int)round(emin)),Form("tree_%d",(int)round(emin)),ntrees);
        emin+=estep;
        ntrees++;
    }

    outfile->cd();

    for (Int_t i=0;i<ntrees;i++){//true (y)
        treeData[i]->Write();
    }

    outfile->Close();
}
void getSingle(char* filename="build-LXe2/g_vary_e/LXe_gamma_Co60.root",char* treename="tree_Co60")
{
    outfile = new TFile("results/calibdataCo60.root","recreate");
    Double_t estep = 200;
    Double_t emin = 200;
    Double_t emax = 2000;
    Int_t ntrees = 0;

    getH1d(filename,treename,0);

    outfile->cd();

    treeData[0]->Write();

    outfile->Close();
}
