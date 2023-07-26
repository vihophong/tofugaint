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


TH1F* htof[MAX_N_HISTO];
TH1F* he_hires[MAX_N_HISTO];

using namespace std;


Double_t GetCorrectionFactor(TH2F* h2in,Double_t low,Double_t step, Int_t nstep,Int_t opt=0){

    TCanvas* cc1=new TCanvas("cc1","cc1",900,700);
    cc1->Divide(2,2);
    cc1->cd(1);
    h2in->Draw("colz");
    Double_t lowi = low;
    TLine* l1=new TLine();
    l1->SetLineColor(2);
    TLine* l2=new TLine();
    l2->SetLineColor(3);

    Double_t xx[200];
    Double_t yy[200];
    for (Int_t i=0;i<nstep;i++){
        cc1->cd(1);
        Int_t startbin=h2in->GetXaxis()->FindBin(lowi);
        Int_t stopbin=h2in->GetXaxis()->FindBin(lowi+step);
        l1->DrawLine(lowi,h2in->GetYaxis()->GetXmin(),lowi,h2in->GetYaxis()->GetXmax());
        l2->DrawLine(lowi+step,h2in->GetYaxis()->GetXmin(),lowi+step,h2in->GetYaxis()->GetXmax());
        cout<<startbin<<"\tsss"<<stopbin<<endl;
        cc1->cd(3);
        TH1F* hproj=(TH1F*) h2in->ProjectionY(Form("prj%d",i),startbin,stopbin);
        hproj->SetLineColor(i);
        hproj->Draw("same");
        hproj->Fit("gaus","LQE","same goff");
        xx[i] = lowi+step/2;
        yy[i] = hproj->GetFunction("gaus")->GetParameter(1);
        lowi+=step;
    }
    cc1->cd(2);
    TGraph* gr = new TGraph(nstep,xx,yy);
    if (opt==0){
        gr->Fit("pol1");
        gr->SetMarkerStyle(20);
        gr->Draw("AP");
        return gr->GetFunction("pol1")->GetParameter(1);
    }else if (opt==1){
        gr->Fit("pol2");
        gr->SetMarkerStyle(20);
        gr->Draw("AP");
        return gr->GetFunction("pol2")->GetParameter(1);
    }else if (opt==2){
        gr->Fit("pol3");
        gr->SetMarkerStyle(20);
        gr->Draw("AP");
        return gr->GetFunction("pol3")->GetParameter(1);
    }
}


void ana(char* infile,char* outfile = "data_Cf252.root")
{
    Double_t dtSum = 22.899045;
    Double_t l0 = 1086;//mm
    Double_t zcorr[] = {-0.076769/0.00770464,1/0.00770464};
    Double_t cc = 0.005227037589038097;

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

    TH1F* h1 = new TH1F("h1","h1",200,0,100);
    TH1F* h2 = new TH1F("h2","h2",200,0,100);
    TH2F* h3 = new TH2F("h3","h3",200,-600,600,200,-30,30);
    TH2F* h4 = new TH2F("h4","h4",200,-600,600,200,-600,600);

    TH1F* h5 = new TH1F("h5","h5",200,0,100);
    TH1F* h6 = new TH1F("h6","h6",200,0,100);


//    TH1F* h7 = new TH1F("h7","h7",200,0,10);

//    TH1F* h7 = new TH1F("h7","h7",61,0,12200);
//    TH1F* h7 = new TH1F("h7","h7",22,0,4400);
    TH1F* h7 = new TH1F("h7","h7",52,0,200*52+100);

    htof[0] = new TH1F("htof","htof",500,0,2000);

    Long64_t nentries = tofu->GetEntries();

    he_hires[0] = new TH1F("hehires","hehires",500,0,2000);

    TH2F* hlvsE = new TH2F("hlvsE","hlvsE",200,0,1000,200,-5000,5000);
    TH2F* htofvsE = new TH2F("htofvsE","htofvsE",200,0,1000,200,0,1000);
    TH2F* h8 = new TH2F("h8","h8",500,0,500,200,-5000,5000);


    TH2F* h9 = new TH2F("h9","h9",500,0,1000,500,0,500);
    TH1F* h10 = new TH1F("h10","h10",1000,0,1000);
    TH2F* h11 = new TH2F("h11","h11",500,0,5.,500,0,1000);

    std::ofstream ofs("tmp.txt");


    Long64_t nbytes = 0;

    Int_t tmp[MAX_N_SCINT*2];
    for (int i=0;i<MAX_N_SCINT*2;i++) tmp[i]=-1;
    for (Long64_t ientry=0; ientry<nentries;ientry++) {
        nbytes += tofu->GetEntry(ientry);
        if (scint_mult!=1) continue;
        Int_t tofupmtMult = 0;
        Int_t TOFUid = 0;

        Int_t TOFUmultID[MAX_N_SCINT*2];
        memcpy(TOFUmultID,tmp,sizeof(tmp));

        for (int i=0;i<MAX_N_SCINT*2;i++) TOFUmultID[i]=-1;

        for (int i=0;i<pmt_mult;i++){
            if (pid[i]<500){
                if (pid[i]>TOFUid)
                    TOFUid=pid[i];
                TOFUmultID[pid[i]] = i;
                tofupmtMult++;
            }
        }
        if (tofupmtMult!=2) continue;
        if (TOFUid%2!=1) continue;
        TOFUid = (TOFUid-1)/2;
        if (TOFUid!=sid[0]){
            cout<<"pmt and scint No. missmatch at Event  No. "<<evt<<endl;
            continue;
        }

        Int_t TOFUmultIDL = TOFUmultID[TOFUid*2];
        Int_t TOFUmultIDR = TOFUmultID[TOFUid*2+1];
        if (TOFUmultIDL<0||TOFUmultIDR<0) cout<<"sth wrong"<<endl;


//        //! find start detector
//        Int_t multIDstart = -1;
//        for (int i=0;i<pmt_mult;i++){
//            if (pid[i]==504){
//                multIDstart=i;
//            }
//        }
//        if (multIDstart<0)
//            continue;

//        if (TOFUid<INNER_RING_ID){
//            h1->Fill((pt[TOFUmultIDL]+pt[TOFUmultIDR])/2-pt[multIDstart]);
//        }

//        if (!(abs(sz[0]+200)<40)) continue;
        if (TOFUid<INNER_RING_ID){
            Double_t tdiff = pt[TOFUmultIDL]-pt[TOFUmultIDR];
            Double_t zposCal = zcorr[1]+tdiff*zcorr[1];
            h3->Fill(sz[0],tdiff);
            h4->Fill(sz[0],zposCal);

//            h1->Fill((tdiff));
//            h2->Fill(sz[0]);
//            h1->Fill(pt[TOFUmultIDL]+pt[TOFUmultIDR]-2*stmin[0]);

            h2->Fill(stmin[0]);
            h1->Fill(pt[TOFUmultIDL]+pt[TOFUmultIDR]-2*stmin[0]);
            if (abs(sz[0]-300)<40){
                h5->Fill(pt[TOFUmultIDL]+pt[TOFUmultIDR]-2*stmin[0]);
            }
            if (abs(sz[0]+300)<40){
                h6->Fill(pt[TOFUmultIDL]+pt[TOFUmultIDR]-2*stmin[0]);
            }

            Double_t tof = (pt[TOFUmultIDL]+pt[TOFUmultIDR])/2-dtSum/2;
            Double_t l = sqrt(zposCal*zposCal+l0*l0);
            Double_t EMeV = cc* (l/tof)* (l/tof);
            h7->Fill(EMeV*1000);
            htof[0]->Fill(tof);
            hlvsE->Fill(EMeV*1000,zposCal);
            he_hires[0]->Fill(EMeV*1000);
            htofvsE->Fill(EMeV*1000,tof);

            Double_t Qsum = pm[TOFUmultIDL]+pm[TOFUmultIDR];
            h10->Fill(Qsum);
//            if (tof>20)
            h11->Fill(se[0],Qsum);
            ofs<<se[0]<<"\t"<<Qsum<<endl;
            h8->Fill(Qsum,zposCal);
//            if (!(abs(sz[0])<40))
                h9->Fill(tof,Qsum);


        }
    }

    TCanvas* c1 =new TCanvas("c1","c1",900,700);
    c1->Divide(2,2);
    c1->cd(1);
//    h1->Draw();
    htof[0]->Draw();
    h5->SetLineColor(2);
    h6->SetLineColor(3);
//    h5->Draw("same");
//    h6->Draw("same");

    c1->cd(2);
//    h7->Draw();
//    he_hires[0]->Draw();
    h10->Draw();
    c1->cd(3);
//    h3->Draw("colz");
//    htofvsE->Draw("colz");
//    h8->Draw("colz");
    h9->Draw("colz");
    c1->cd(4);
//    h4->Draw("colz");
//    hlvsE->Draw("colz");
    h11->Draw("colz");
    c1->Draw();

//    GetCorrectionFactor(h3,-400,160,5);

    h7->SaveAs(outfile);
}


TH1F* getH1d(char* infile,char* hname,Int_t nbins,Double_t emax,Int_t nn)
{
    Double_t dtSum = 22.899045;
    Double_t l0 = 1086;//mm
    Double_t zcorr[] = {-0.076769/0.00770464,1/0.00770464};
    Double_t cc = 0.005227037589038097;

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

    TH1F* h1 = new TH1F("h1","h1",200,0,100);
    TH1F* h2 = new TH1F("h2","h2",200,0,100);
    TH2F* h3 = new TH2F("h3","h3",200,-600,600,200,-15,15);
    TH2F* h4 = new TH2F("h4","h4",200,-600,600,200,-600,600);

    TH1F* h5 = new TH1F("h5","h5",200,0,100);
    TH1F* h6 = new TH1F("h6","h6",200,0,100);

    htof[nn] = new TH1F(Form("htof%s",hname),Form("htof%s",hname),500,0,2000);


    TH1F* h7 = new TH1F("h7","h7",nbins,0,emax);
    he_hires[nn] = new TH1F(Form("hEhires%s",hname),Form("hEhires%s",hname),500,0,12000);
//    cout<<"EEEE"<<h7->GetBinCenter(1)<<endl;
    Long64_t nentries = tofu->GetEntries();
    Long64_t nbytes = 0;

//    TH1F* h11 = new TH1F("h11","h11",61,0,12200);
//    for (Int_t i=0;i<61;i++)
//        cout<<h11->GetBinCenter(i+1)<<endl;

    Int_t tmp[MAX_N_SCINT*2];
    for (int i=0;i<MAX_N_SCINT*2;i++) tmp[i]=-1;
    for (Long64_t ientry=0; ientry<nentries;ientry++) {
        nbytes += tofu->GetEntry(ientry);
        if (scint_mult!=1) continue;
        Int_t tofupmtMult = 0;
        Int_t TOFUid = 0;

        Int_t TOFUmultID[MAX_N_SCINT*2];
        memcpy(TOFUmultID,tmp,sizeof(tmp));

        for (int i=0;i<MAX_N_SCINT*2;i++) TOFUmultID[i]=-1;

        for (int i=0;i<pmt_mult;i++){
            if (pid[i]<500){
                if (pid[i]>TOFUid)
                    TOFUid=pid[i];
                TOFUmultID[pid[i]] = i;
                tofupmtMult++;
            }
        }
        if (tofupmtMult!=2) continue;
        if (TOFUid%2!=1) continue;
        TOFUid = (TOFUid-1)/2;
        if (TOFUid!=sid[0]){
            cout<<"pmt and scint No. missmatch at Event  No. "<<evt<<endl;
            continue;
        }

        Int_t TOFUmultIDL = TOFUmultID[TOFUid*2];
        Int_t TOFUmultIDR = TOFUmultID[TOFUid*2+1];
        if (TOFUmultIDL<0||TOFUmultIDR<0) cout<<"sth wrong"<<endl;

//        //! find start detector
//        Int_t multIDstart = -1;
//        for (int i=0;i<pmt_mult;i++){
//            if (pid[i]==504){
//                multIDstart=i;
//            }
//        }
//        if (multIDstart<0)
//            continue;

//        if (TOFUid<INNER_RING_ID){
//            h1->Fill((ptmin[TOFUmultIDL]+ptmin[TOFUmultIDR])/2-ptmin[multIDstart]);
//        }

//        if (!(abs(sz[0]+200)<40)) continue;
        if (TOFUid<INNER_RING_ID){
            Double_t tdiff = pt[TOFUmultIDL]-pt[TOFUmultIDR];
            Double_t zposCal = zcorr[1]+tdiff*zcorr[1];
            h3->Fill(sz[0],tdiff);
            h4->Fill(sz[0],zposCal);

//            h1->Fill((tdiff));
//            h2->Fill(sz[0]);
//            h1->Fill(pt[TOFUmultIDL]+pt[TOFUmultIDR]-2*stmin[0]);

            h2->Fill(stmin[0]);
            h1->Fill(pt[TOFUmultIDL]+pt[TOFUmultIDR]-2*stmin[0]);
            if (abs(sz[0]-300)<40){
                h5->Fill(pt[TOFUmultIDL]+pt[TOFUmultIDR]-2*stmin[0]);
            }
            if (abs(sz[0]+300)<40){
                h6->Fill(pt[TOFUmultIDL]+pt[TOFUmultIDR]-2*stmin[0]);
            }

            Double_t tof = (pt[TOFUmultIDL]+pt[TOFUmultIDR])/2-dtSum/2;
            htof[nn]->Fill(tof);
            Double_t l = sqrt(zposCal*zposCal+l0*l0);
            Double_t EkeV = cc* (l/tof)* (l/tof)*1000;
            h7->Fill(EkeV);
            he_hires[nn]->Fill(EkeV);
        }
    }

//    TCanvas* c1 =new TCanvas("c1","c1",900,700);
//    c1->Divide(2,2);
//    c1->cd(1);
//    h1->Draw();
//    h5->SetLineColor(2);
//    h6->SetLineColor(3);
////    h5->Draw("same");
////    h6->Draw("same");

//    c1->cd(2);
//    h7->Draw();
//    c1->cd(3);
//    h3->Draw("colz");
//    c1->cd(4);
//    h4->Draw("colz");
//    c1->Draw();

//    GetCorrectionFactor(h3,-400,160,5);
    h7->SetName(Form("hE%s",hname));
//    f->Close();
    return h7;
}


void det_response()
{
    TH1F* hh[MAX_N_HISTO];
    Int_t nbins = 52;
    Double_t emax = 200*nbins;

    TFile* outfile = new TFile("response_matrix.root","recreate");

    TCanvas *c1 = new TCanvas("c1","c1",900,700);
    c1->Divide(2,1);
    c1->cd(1);
    TH2F* hmatrix = new TH2F("hmatrix","hmatrix",nbins,0,emax,nbins,0,emax);
    for (Int_t i=0;i<nbins;i++){//true (y)
        hh[i] = getH1d(Form("vary_ene/LXe_%d.root",200*i+100),Form("%d",200*i+100),nbins,emax,i);
        hh[i]->SetLineColor(i);
        if (i==1)
            hh[i]->Draw();
        else
            hh[i]->Draw("same");
//        hmatrix->SetBinContent(i+1,)
        for (Int_t j=0;j<nbins;j++){//observe (x)
            hmatrix->SetBinContent(j+1,i+1,hh[i]->GetBinContent(j+1));
        }

        cout<<200*i+100<<"\t"<<hh[i]->GetEntries()<<"\t"<<htof[i]->GetEntries()<<endl;

    }
    c1->cd(2);
    hmatrix->Draw("colz");


    outfile->cd();

    for (Int_t i=0;i<nbins;i++){//true (y)
        hh[i]->Write();
        htof[i]->Write();
        he_hires[i]->Write();
    }
    hmatrix->Write();

    outfile->Close();
}
