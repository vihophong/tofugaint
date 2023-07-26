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
//GARI.hh
//
/// \file optical/LXe/src/LXeDetectorConstruction.cc
/// \brief Implementation of the LXeDetectorConstruction class
//
//
#include "LXeDetectorConstruction.hh"
#include "LXeDetectorMessenger.hh"
#include "LXeMainVolume.hh"
#include "LXePMTSD.hh"
#include "LXeScintSD.hh"
#include "LXeWLSSlab.hh"

#include "globals.hh"
#include "G4Box.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4OpticalSurface.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4PVPlacement.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SolidStore.hh"
#include "G4Sphere.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4UImanager.hh"
#include "G4VisAttributes.hh"
#include "G4UnionSolid.hh"
#include "G4NistManager.hh"
#include "G4CutTubs.hh"

#include "G4SubtractionSolid.hh"




#include "Clover.hh"

#include "G4UnitsTable.hh"

#include <fstream>

#define NOGDML 1
#define SUSFRAMEONLY 1

#define WGARI 1
//#define WCLOVER 1

G4bool LXeDetectorConstruction::fSphereOn = true;

using namespace std;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeDetectorConstruction::LXeDetectorConstruction()
    : fParser(), fLXe_mt(nullptr)
    , fMPTPStyrene(nullptr)
{
    fWorld_box  = nullptr;
    fWorld_log  = nullptr;
    fWorld_phys = nullptr;

    fLXe = fAl = fAir = fVacuum = fGlass = nullptr;
    fPstyrene = fPMMA = fPethylene1 = fPethylene2 = nullptr;

    fN = fO = fC = fH = nullptr;

    fSaveThreshold = 0;
    SetDefaults();

    DefineMaterials();
    fDetectorMessenger = new LXeDetectorMessenger(this);


    //! VIS attributes
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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeDetectorConstruction::~LXeDetectorConstruction()
{
    if(fMainVolume)
    {
        delete fMainVolume;
    }
    delete fLXe_mt;
    delete fDetectorMessenger;
    delete fMPTPStyrene;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::DefineMaterials()
{

    G4NistManager* nist = G4NistManager::Instance();
    G4double a;  // atomic mass
    G4double z;  // atomic number
    G4double density;

    G4int polyPMMA = 1;
    G4int nC_PMMA  = 3 + 2 * polyPMMA;
    G4int nH_PMMA  = 6 + 2 * polyPMMA;

    G4int polyeth = 1;
    G4int nC_eth  = 2 * polyeth;
    G4int nH_eth  = 4 * polyeth;

    //***Elements
    fH = nist->FindOrBuildElement("H");
    fC = nist->FindOrBuildElement("C");
    fN = nist->FindOrBuildElement("N");
    fO = nist->FindOrBuildElement("O");

    //***Materials

    // EJ200 scintillator
    fLXe = new G4Material("EJ200", density = 1.023*g/cm3, 2);
    fLXe->AddElement(nist->FindOrBuildElement("H"), 10);
    fLXe->AddElement(nist->FindOrBuildElement("C"), 9);
    //

    //    fLXe->AddElement(nist->FindOrBuildElement("H"), 0.08457);//H10 mass fraction
    //    fLXe->AddElement(nist->FindOrBuildElement("C"), 0.91543);//C9 mass fraction


    // STL
    STL_NIST = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
    STL_NIST->SetName("STL_NIST");

    // Bialkali K2CsSb Photocathode Material
    fBialkali = new G4Material("Bialkali", density = 2.*g/cm3, 3);
    G4Element* K = nist->FindOrBuildElement("K");
    G4Element* Cs = nist->FindOrBuildElement("Cs");
    G4Element* Sb = nist->FindOrBuildElement("Sb");
    fBialkali->AddElement(K,  2);
    fBialkali->AddElement(Cs, 1);
    fBialkali->AddElement(Sb, 1);

    // Borosilicate grass window
    G4Material* Borosilicate = new G4Material("Borosilicate_Glass",density= 2.23*g/cm3, 5);
    G4Material* Al2O3      = nist->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");
    G4Material* B2O3       = nist->FindOrBuildMaterial("G4_BORON_OXIDE");
    G4Material* SiO2       = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
    G4Material* Na2O       = nist->FindOrBuildMaterial("G4_SODIUM_MONOXIDE");
    G4Material* K2O        = nist->FindOrBuildMaterial("G4_POTASSIUM_OXIDE");
    Borosilicate->AddMaterial(SiO2,   80.6 * perCent);
    Borosilicate->AddMaterial(B2O3,  13.0 * perCent);
    Borosilicate->AddMaterial(Na2O,  2.   * perCent); // 1/2 of wt% for (Na20+K20)
    Borosilicate->AddMaterial(K2O,   2.   * perCent); // 1/2 of wt% for (Na20+K20)
    Borosilicate->AddMaterial(Al2O3, 2.31  * perCent);

    G4MaterialPropertiesTable* BorosilicateMPT = new G4MaterialPropertiesTable();
    // Assign refractive index  of quazt glass to the Borosilicate grass window of the PMT
    //"Refractiveindex.info"
    std::vector<G4double> Ephoton_quartz = {
      1.77120301*eV, 1.796872619*eV, 1.823297216*eV, 1.850510608*eV, 1.878548647*eV,
      1.907449396*eV, 1.937253292*eV, 1.968003345*eV, 1.999745334*eV, 2.032528044*eV,
      2.066403512*eV, 2.1014273*eV, 2.137658805*eV, 2.175161591*eV, 2.214003763*eV,
      2.254258377*eV, 2.296003902*eV, 2.33932473*eV, 2.384311744*eV, 2.431062955*eV,
      2.479684214*eV, 2.530290015*eV, 2.58300439*eV, 2.63796193*eV, 2.695308928*eV,
      2.755204682*eV, 2.817822971*eV, 2.883353737*eV, 2.952005017*eV, 3.024005139*eV,
      3.099605268*eV, 3.179082326*eV, 3.262742387*eV, 3.350924614*eV, 3.444005853*eV,
      3.54240602*eV, 3.646594433*eV, 3.757097294*eV, 3.874506585*eV, 3.999490668*eV,
      4.132807024*eV, 4.275317611*eV, 4.428007525*eV, 4.592007804*eV, 4.768623489*eV,
      4.959368428*eV, 5.16600878*eV, 5.390617857*eV, 5.635645941*eV, 5.904010034*eV,
      6.199210536*eV
    };

    std::vector<G4double> Rindex_quartz = {
      1.455292466, 1.455524071, 1.455763571, 1.456011496, 1.456268423, 1.456534974, 1.456811819, 1.457099689,
      1.457399374, 1.457711733, 1.458037702, 1.4583783, 1.458734641, 1.459107942, 1.459499536, 1.459910886,
      1.460343603, 1.460799458, 1.461280408, 1.461788618, 1.462326487, 1.462896682, 1.463502175, 1.464146283,
      1.464832722, 1.465565665, 1.466349815, 1.467190482, 1.46809369, 1.469066293, 1.470116119, 1.471252144,
      1.472484709, 1.473825777, 1.475289258, 1.476891413, 1.478651361, 1.48059172, 1.482739429, 1.485126813,
      1.487792976, 1.490785646, 1.494163661, 1.498000361, 1.502388312, 1.507446007, 1.513327606, 1.520237459,
      1.528452449, 1.53835762, 1.550505538
    };
    BorosilicateMPT->AddProperty("RINDEX", Ephoton_quartz, Rindex_quartz);
    // Absorption according to NEXTSim
    std::vector<G4double> PhotonEnergy_Borosilicate = { 2.484*eV, 2.615*eV, 2.760*eV, 2.922*eV, 3.105*eV };
    std::vector<G4double> Absorption_Borosilicate =  {125.*cm, 123.5*cm, 122.*cm, 121.*cm, 120.*cm};
    BorosilicateMPT->AddProperty("ABSLENGTH", PhotonEnergy_Borosilicate, Absorption_Borosilicate);
    Borosilicate->SetMaterialPropertiesTable(BorosilicateMPT);


    // Aluminum
    fAl = new G4Material("Al", z = 13., a = 26.98 * g / mole,
                         density = 2.7 * g / cm3);

    //! Set Al as a reflector
    std::vector<G4double> al_ephoton = { 2.004*eV, 2.615*eV, 2.760*eV, 2.922*eV, 3.191*eV };
    std::vector<G4double> reflectivity     = { fRefl, fRefl,fRefl, fRefl,fRefl };
    std::vector<G4double> efficiency       = { 0.0, 0.0,0.0, 0.0,0.0 };
    G4MaterialPropertiesTable* alMPT = new G4MaterialPropertiesTable();
    alMPT->AddProperty("REFLECTIVITY", al_ephoton, reflectivity);
    alMPT->AddProperty("EFFICIENCY", al_ephoton, efficiency);
    fAl->SetMaterialPropertiesTable(alMPT);
    STL_NIST->SetMaterialPropertiesTable(alMPT);


    G4Material* Tapemat = nist->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");
    Tapemat->SetName("TAPE");
    //! outer tape as an absorber
    std::vector<G4double> ephoton = { 2.484*eV, 2.615*eV, 2.760*eV, 2.922*eV, 3.105*eV };
    std::vector<G4double> tape_EFF     = { 1., 1.,1., 1.,1. };
    std::vector<G4double> tape_ReR     = { 1.92, 1.92,1.92, 1.92,1.92 };
    std::vector<G4double> tape_ImR     = { 1.69, 1.69,1.69, 1.69,1.69 };
    std::vector<G4double>  tape_Reflectivity = { 0., 0.,0., 0.,0.};
    G4MaterialPropertiesTable* tape_mt = new G4MaterialPropertiesTable();
    tape_mt->AddProperty("EFFICIENCY", ephoton, tape_EFF);
    tape_mt->AddProperty("REFLECTIVITY", ephoton, tape_Reflectivity);
    tape_mt->AddProperty("REALRINDEX", ephoton, tape_ReR);
    tape_mt->AddProperty("IMAGINARYRINDEX", ephoton, tape_ImR);
    Tapemat->SetMaterialPropertiesTable(tape_mt);


    // Vacuum
    fVacuum = new G4Material("Vacuum", z = 1., a = 1.01 * g / mole,
                             density = universe_mean_density, kStateGas,
                             0.1 * kelvin, 1.e-19 * pascal);
    // Air
    fAir = new G4Material("Air", density = 1.29 * mg / cm3, 2);
    fAir->AddElement(fN, 70 * perCent);
    fAir->AddElement(fO, 30 * perCent);


    // Glass
    fGlass = new G4Material("Glass", density = 1.032 * g / cm3, 2);
    fGlass->AddElement(fC, 91.533 * perCent);
    fGlass->AddElement(fH, 8.467 * perCent);
    // Polystyrene
    fPstyrene = new G4Material("Polystyrene", density = 1.03 * g / cm3, 2);
    fPstyrene->AddElement(fC, 8);
    fPstyrene->AddElement(fH, 8);
    // Fiber(PMMA)
    fPMMA = new G4Material("PMMA", density = 1190. * kg / m3, 3);
    fPMMA->AddElement(fH, nH_PMMA);
    fPMMA->AddElement(fC, nC_PMMA);
    fPMMA->AddElement(fO, 2);
    // Cladding(polyethylene)
    fPethylene1 = new G4Material("Pethylene1", density = 1200. * kg / m3, 2);
    fPethylene1->AddElement(fH, nH_eth);
    fPethylene1->AddElement(fC, nC_eth);
    // Double cladding(flourinated polyethylene)
    fPethylene2 = new G4Material("Pethylene2", density = 1400. * kg / m3, 2);
    fPethylene2->AddElement(fH, nH_eth);
    fPethylene2->AddElement(fC, nC_eth);

    //***Material properties tables
    fLXe_mt = new G4MaterialPropertiesTable();

    G4double photonEnergy_Ej200[44] = {2.004*eV, 2.058*eV, 2.112*eV, 2.166*eV, 2.220*eV, 2.274*eV, 2.328*eV, 2.382*eV, 2.436*eV, 2.490*eV,
                                       2.517*eV, 2.552*eV, 2.585*eV, 2.613*eV, 2.635*eV, 2.656*eV, 2.686*eV, 2.720*eV, 2.749*eV, 2.772*eV,
                                       2.791*eV, 2.809*eV, 2.826*eV, 2.842*eV, 2.861*eV, 2.884*eV, 2.919*eV, 2.946*eV, 2.954*eV, 2.961*eV,
                                       2.967*eV, 2.974*eV, 2.981*eV, 2.987*eV, 2.994*eV, 3.001*eV, 3.009*eV, 3.018*eV, 3.029*eV, 3.041*eV,
                                       3.056*eV, 3.083*eV, 3.137*eV, 3.191*eV};

    G4double ScintilFast_EJ200[44] = {0.000, 0.001, 0.001, 0.002, 0.003, 0.006, 0.010, 0.018, 0.033, 0.060,
                                      0.084, 0.122, 0.175, 0.234, 0.294, 0.356, 0.416, 0.473, 0.533, 0.594,
                                      0.657, 0.720, 0.784, 0.846, 0.903, 0.962, 1.000, 0.917, 0.857, 0.798,
                                      0.732, 0.669, 0.604, 0.542, 0.480, 0.422, 0.359, 0.297, 0.237, 0.170,
                                      0.105, 0.028, 0.004, 0.000};

    G4double photonEnergy_Ej200_2[2] = {2.004*eV, 3.191*eV};
    G4double RefIndex_EJ200[2] = {1.580, 1.580};
    G4double Absorption_EJ200[2] = {400*cm, 400*cm};


    fLXe_mt->AddProperty("SCINTILLATIONCOMPONENT1", photonEnergy_Ej200, ScintilFast_EJ200, 44,true,true);
//    fLXe_mt->AddProperty("SCINTILLATIONCOMPONENT2", photonEnergy_Ej200, ScintilFast_EJ200, 44,true,true);
    fLXe_mt->AddConstProperty("SCINTILLATIONRISETIME1", 0.9 * ns,true);
//    fLXe_mt->AddConstProperty("SCINTILLATIONRISETIME2", 0.9 * ns,true);
    fLXe_mt->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 2.1 * ns,true);
//    fLXe_mt->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 2.1 * ns,true);
//    fLXe_mt->AddConstProperty("SCINTILLATIONYIELD1", 1.0,true);
//    fLXe_mt->AddConstProperty("SCINTILLATIONYIELD2", 0.0,true);

    fLXe_mt->AddProperty("RINDEX", photonEnergy_Ej200_2, RefIndex_EJ200, 2,true,false);
    fLXe_mt->AddProperty("ABSLENGTH", photonEnergy_Ej200_2, Absorption_EJ200, 2,true,false);


    fLXe_mt->AddConstProperty("SCINTILLATIONYIELD", 10000/MeV,true); // Scintillation efficiency as per Eljen specs
    fLXe_mt->AddConstProperty("RESOLUTIONSCALE", 1.0,true); // Intrinsic resolution


    G4double fE = 1.;

    G4double fP = 1.;

    std::vector<G4double> pEne = {0*MeV, 0.1*MeV, 0.13*MeV, 0.17*MeV, 0.2*MeV, 0.24*MeV, 0.3*MeV, 0.34*MeV, 0.4*MeV, 0.48*MeV, 0.6*MeV, 0.72*MeV, 0.84*MeV, 1*MeV, 1.3*MeV, 1.7*MeV, 2*MeV, 2.4*MeV, 3*MeV, 3.4*MeV, 4*MeV, 4.8*MeV, 6*MeV, 7.2*MeV, 8.4*MeV, 10*MeV, 13*MeV, 17*MeV, 20*MeV, 24*MeV, 30*MeV, 34*MeV, 40*MeV};

    std::vector<G4double> eYields = {0*fE, 1000*fE, 1300*fE, 1700*fE, 2000*fE, 2400*fE, 3000*fE, 3400*fE, 4000*fE, 4800*fE, 6000*fE, 7200*fE, 8400*fE, 10000*fE, 13000*fE, 17000*fE, 20000*fE, 24000*fE, 30000*fE, 34000*fE, 40000*fE, 48000*fE, 60000*fE, 72000*fE, 84000*fE, 100000*fE, 130000*fE, 170000*fE, 200000*fE, 240000*fE, 300000*fE, 340000*fE, 400000*fE};

    std::vector<G4double> pYields = {0*fP, 71.2407*fP, 94.0675*fP, 128.148*fP, 155.541*fP, 195.142*fP, 261.181*fP, 307.896*fP, 387.524*fP, 512.806*fP, 719.839*fP, 966.156*fP, 1247.51*fP, 1658.39*fP, 2532.18*fP, 3885.86*fP, 5016.58*fP, 6635.69*fP, 9194.41*fP, 11063*fP, 14088.9*fP, 18240.2*fP, 24525.5*fP, 31320.4*fP, 38433.9*fP, 48307.8*fP, 67524.8*fP, 93749*fP, 114665*fP, 143331*fP, 187923*fP, 217651*fP, 263304*fP};

    std::vector<G4double> aYields = {0*fP, 17.412*fP, 22.1897*fP, 28.8785*fP, 33.9747*fP, 40.982*fP, 52.0238*fP, 59.8804*fP, 71.6654*fP, 88.1219*fP, 114.665*fP, 143.331*fP, 175.819*fP, 222.959*fP, 320.636*fP, 468.214*fP, 596.681*fP, 796.282*fP, 1167.88*fP, 1449.23*fP, 1927*fP, 2712.67*fP, 4321.16*fP, 6444.58*fP, 9236.88*fP, 14014.6*fP, 24950.2*fP, 42786.9*fP, 57757*fP, 78672.7*fP, 110630*fP, 132077*fP, 164565*fP};

    std::vector<G4double> ionYields = {0*fP, 11.0205*fP, 13.4837*fP, 16.7007*fP, 18.9834*fP, 22.0411*fP, 26.6065*fP, 29.6536*fP, 33.8792*fP, 39.0285*fP, 46.3012*fP, 53.3297*fP, 60.3688*fP, 69.7437*fP, 86.2958*fP, 107.838*fP, 123.657*fP, 144.754*fP, 176.403*fP, 198.678*fP, 232.079*fP, 276.618*fP, 343.431*fP, 411.413*fP, 479.405*fP, 573.175*fP, 757.488*fP, 1049.05*fP, 1289.34*fP, 1629.26*fP, 2191.88*fP, 2613.85*fP, 3317.12*fP};

    fLXe_mt->AddProperty("ELECTRONSCINTILLATIONYIELD", pEne, eYields,true,true);

    fLXe_mt->AddProperty("PROTONSCINTILLATIONYIELD", pEne, pYields, true,true);
    fLXe_mt->AddProperty("DEUTERONSCINTILLATIONYIELD", pEne, pYields, true,true);
    fLXe_mt->AddProperty("TRITONSCINTILLATIONYIELD", pEne, pYields, true,true);

    fLXe_mt->AddProperty("ALPHASCINTILLATIONYIELD", pEne, aYields, true,true);

    fLXe_mt->AddProperty("IONSCINTILLATIONYIELD",  pEne, ionYields, true,true);

    fLXe->SetMaterialPropertiesTable(fLXe_mt);

    // Set the Birks Constant for the LXe scintillator (not needed?)
//    fLXe->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);

    std::vector<G4double> glass_Energy = { 2.004*eV, 3.191*eV };
    std::vector<G4double> glass_AbsLength = { 420. * cm, 420. * cm};
    G4MaterialPropertiesTable* glass_mt   = new G4MaterialPropertiesTable();
    glass_mt->AddProperty("ABSLENGTH", glass_Energy, glass_AbsLength);
    glass_mt->AddProperty("RINDEX", "Fused Silica");
    fGlass->SetMaterialPropertiesTable(glass_mt);

    G4MaterialPropertiesTable* vacuum_mt = new G4MaterialPropertiesTable();
    vacuum_mt->AddProperty("RINDEX", "Air");
    fVacuum->SetMaterialPropertiesTable(vacuum_mt);
    fAir->SetMaterialPropertiesTable(vacuum_mt);  // Give air the same rindex

    std::vector<G4double> wls_Energy = { 2.00 * eV, 2.87 * eV, 2.90 * eV,
                                         3.47 * eV };

    std::vector<G4double> rIndexPstyrene = { 1.5, 1.5, 1.5, 1.5 };
    std::vector<G4double> absorption1    = { 2. * cm, 2. * cm, 2. * cm, 2. * cm };
    std::vector<G4double> scintilFast    = { 0.0, 0.0, 1.0, 1.0 };
    fMPTPStyrene = new G4MaterialPropertiesTable();
    fMPTPStyrene->AddProperty("RINDEX", wls_Energy, rIndexPstyrene);
    fMPTPStyrene->AddProperty("ABSLENGTH", wls_Energy, absorption1);
    fMPTPStyrene->AddProperty("SCINTILLATIONCOMPONENT1", wls_Energy, scintilFast);
    fMPTPStyrene->AddConstProperty("SCINTILLATIONYIELD", 10. / keV);
    fMPTPStyrene->AddConstProperty("RESOLUTIONSCALE", 1.0);
    fMPTPStyrene->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 10. * ns);
    fPstyrene->SetMaterialPropertiesTable(fMPTPStyrene);

    // Set the Birks Constant for the Polystyrene scintillator
    fPstyrene->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);

    std::vector<G4double> AbsFiber    = { 9.0 * m, 9.0 * m, 0.1 * mm, 0.1 * mm };
    std::vector<G4double> EmissionFib = { 1.0, 1.0, 0.0, 0.0 };
    G4MaterialPropertiesTable* fiberProperty = new G4MaterialPropertiesTable();
    fiberProperty->AddProperty("RINDEX", "PMMA");
    fiberProperty->AddProperty("WLSABSLENGTH", wls_Energy, AbsFiber);
    fiberProperty->AddProperty("WLSCOMPONENT", wls_Energy, EmissionFib);
    fiberProperty->AddConstProperty("WLSTIMECONSTANT", 0.5 * ns);
    fPMMA->SetMaterialPropertiesTable(fiberProperty);

    std::vector<G4double> RefractiveIndexClad1 = { 1.49, 1.49, 1.49, 1.49 };
    G4MaterialPropertiesTable* clad1Property   = new G4MaterialPropertiesTable();
    clad1Property->AddProperty("RINDEX", wls_Energy, RefractiveIndexClad1);
    clad1Property->AddProperty("ABSLENGTH", wls_Energy, AbsFiber);
    fPethylene1->SetMaterialPropertiesTable(clad1Property);

    std::vector<G4double> RefractiveIndexClad2 = { 1.42, 1.42, 1.42, 1.42 };
    G4MaterialPropertiesTable* clad2Property   = new G4MaterialPropertiesTable();
    clad2Property->AddProperty("RINDEX", wls_Energy, RefractiveIndexClad2);
    clad2Property->AddProperty("ABSLENGTH", wls_Energy, AbsFiber);
    fPethylene2->SetMaterialPropertiesTable(clad2Property);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* LXeDetectorConstruction::Construct()
{
#ifdef NOGDML
    // The experimental hall walls are all 1m away from housing walls
    G4double world_x = fExpHall_x;// + 1. * cm;
    G4double world_y = fExpHall_y;// + 1. * cm;
    G4double world_z = fExpHall_z;// + 1. * cm;

    // Create experimental hall
    fWorld_box =
            new G4Box("expHall_box", world_x, world_y, world_z);
    fWorld_log =
            new G4LogicalVolume(fWorld_box, fVacuum, "expHall_log", 0, 0, 0);


    fWorld_phys = new G4PVPlacement(
                0, G4ThreeVector(), fWorld_log, "expHall", 0, false, 0);
    fCheckOverlaps = 0;

#else

#ifdef SUSFRAMEONLY
    fParser.Read("gdml/mother_simpl.gdml");
#else
    fParser.Read("gdml/mother_nohead.gdml");
#endif
    fWorld_phys = fParser.GetWorldVolume();
    fWorld_log = fWorld_phys->GetLogicalVolume();

    G4LogicalVolume* pmtheadcomp1 = fParser.GetVolume("pmthead1");
    G4LogicalVolume* pmtheadcomp2 = fParser.GetVolume("pmthead2");
    G4LogicalVolume* pmtheadcomp3 = fParser.GetVolume("pmthead3");

    G4VPhysicalVolume* pmtheadcomp1phys = fParser.GetPhysVolume("pmthead1_phys");
    G4VPhysicalVolume* pmtheadcomp2phys = fParser.GetPhysVolume("pmthead2_phys");
    G4VPhysicalVolume* pmtheadcomp3phys = fParser.GetPhysVolume("pmthead3_phys");

    G4UnionSolid* solidPmthead12 = new G4UnionSolid("pmthead12",pmtheadcomp1->GetSolid(),pmtheadcomp2->GetSolid(),0,G4ThreeVector(0.*mm, 0.*mm, 0.*mm));
    G4UnionSolid* solidPmthead123 = new G4UnionSolid("pmthead123",solidPmthead12,pmtheadcomp3->GetSolid(),0,G4ThreeVector(0.*mm, 0.*mm, 0.*mm));
    G4LogicalVolume* logicPmthead123 = new G4LogicalVolume(solidPmthead123, STL_NIST,"logicPmthead123");            //its name
#endif

    fWorld_log->SetVisAttributes(G4VisAttributes::GetInvisible());


    // Place the main volume
    if(fMainVolumeOn)
    {
        fMainVolume = new LXeMainVolume(0, G4ThreeVector(), fWorld_log,
                                        false, 0, this);
    }

    // Miscellaneous structure

#ifndef NOGDML

    //---- Place HOLDER in moderator----

    ifstream infile1("holdermapout.txt");
    G4int nholder = 0;
    G4double h_x[1000];
    G4double h_y[1000];
    G4double h_rot[1000];

    while (!infile1.eof()){
        G4int myid,tmp;
        G4double tmp2;
        infile1>>myid;
        infile1>>h_x[myid]>>h_y[myid]>>tmp2>>h_rot[myid]>>tmp;
        nholder++;
    }
    infile1.close();


    G4int holdercpy = 0;
    for (G4int i=0;i<nholder;i++)
    {
        G4RotationMatrix* rotmat = new G4RotationMatrix;
        rotmat->set(0,0,0);
        rotmat->rotateZ(90*degree-h_rot[i]*degree);
        if (holdercpy==0){
            pmtheadcomp1phys->SetTranslation(G4ThreeVector(h_x[i]*mm,h_y[i]*mm,390*mm));
            pmtheadcomp2phys->SetTranslation(G4ThreeVector(h_x[i]*mm,h_y[i]*mm,390*mm));
            pmtheadcomp3phys->SetTranslation(G4ThreeVector(h_x[i]*mm,h_y[i]*mm,390*mm));
            pmtheadcomp1phys->SetRotation(rotmat);
            pmtheadcomp2phys->SetRotation(rotmat);
            pmtheadcomp3phys->SetRotation(rotmat);
        }else{
            new G4PVPlacement(rotmat,
                              G4ThreeVector(h_x[i]*mm,h_y[i]*mm,390*mm),             //position
                              "PhysicsHead",             //its name
                              logicPmthead123,            //its logical volume
                              fMainVolume,             //its mother  volume
                              false,                 //no boolean operation
                              holdercpy,                 //copy number
                              fCheckOverlaps);       // checking overlaps
        }
        holdercpy++;
        G4RotationMatrix* rotmat2 = new G4RotationMatrix;
        rotmat2->set(0,0,0);
        rotmat2->rotateZ(90*degree-h_rot[i]*degree);
        rotmat2->rotateY(180*degree);
        new G4PVPlacement(rotmat2,
                          G4ThreeVector(h_x[i]*mm,h_y[i]*mm,-390*mm),             //position
                          "PhysicsHead",             //its name
                          logicPmthead123,            //its logical volume
                          fMainVolume,             //its mother  volume
                          false,                 //no boolean operation
                          holdercpy,                 //copy number
                          fCheckOverlaps);       // checking overlaps
        holdercpy++;
    }

    ifstream infile2("holdermapin.txt");
    nholder = 0;
    while (!infile2.eof()){
        G4int myid,tmp;
        G4double tmp2;
        infile2>>myid;
        infile2>>h_x[myid]>>h_y[myid]>>tmp2>>h_rot[myid]>>tmp;
        nholder++;
    }
    infile2.close();

    for (G4int i=0;i<nholder;i++)
    {
        G4RotationMatrix* rotmat = new G4RotationMatrix;
        rotmat->set(0,0,0);
        rotmat->rotateZ(90*degree-h_rot[i]*degree);
        rotmat->rotateY(180*degree);
        rotmat->rotateX(180*degree);
        new G4PVPlacement(rotmat,
                          G4ThreeVector(h_x[i]*mm,h_y[i]*mm,390*mm),             //position
                          "PhysicsHead",             //its name
                           logicPmthead123,            //its logical volume
                           fMainVolume,             //its mother  volume
                           false,                 //no boolean operation
                           holdercpy,                 //copy number
                           fCheckOverlaps);       // checking overlaps

        holdercpy++;
        G4RotationMatrix* rotmat2 = new G4RotationMatrix;
        rotmat2->set(0,0,0);
        rotmat2->rotateZ(90*degree-h_rot[i]*degree);
        rotmat2->rotateX(180*degree);
        new G4PVPlacement(rotmat2,
                          G4ThreeVector(h_x[i]*mm,h_y[i]*mm,-390*mm),             //position
                          "PhysicsHead",             //its name
                          logicPmthead123,            //its logical volume
                          fMainVolume,             //its mother  volume
                          false,                 //no boolean operation
                          holdercpy,                 //copy number
                          fCheckOverlaps);       // checking overlaps
        holdercpy++;
    }
#endif

#ifdef WGARI
    fgari_det=new GARI();
    fgari_det->SetCopyNo(500);
    //account for center of logical gari
//    G4double offsetYZ_GARI = 2*cm;
    G4double offsetYZ_GARI = fgari_det->GetOffsetSurfaceGlass();
    fgari_det->SetPosition(G4ThreeVector(0,offsetYZ_GARI,offsetYZ_GARI));

    G4RotationMatrix* rotmatGARI = new G4RotationMatrix;
//    rotmatGARI->set(0,0,0);
    rotmatGARI->rotateX((45)*degree);
    fgari_det->SetRotation(*rotmatGARI);
    fgari_det->Placement(500,fMainVolume);

    G4double gariMaxRadius = fgari_det->GetGlassMaxRadius()+2*mm;
    G4double beamPipeLength = 500*mm;
    G4double beamPipeThickness = 2*mm;
    G4double beamPipeOffset = beamPipeLength/2+fgari_det->GetSizeZ()/2/cos(45*pi/180.);
    G4CutTubs* solidBeamPipe = new G4CutTubs("solidBeamPipe",gariMaxRadius-beamPipeThickness,gariMaxRadius,beamPipeLength/2,0,360.*deg,
                                             G4ThreeVector(0,0.,-1),G4ThreeVector(cos(45*pi/180.),0,cos(45*pi/180.)));
    G4LogicalVolume* logicBeamPipe = new G4LogicalVolume(solidBeamPipe,          //its solid
                                    STL_NIST,         //its material
                                    "logicBeamPipe");        //its name
    logicBeamPipe->SetVisAttributes(ColorMap["invi"]);
//    logicBeamPipe->SetVisAttributes(ColorMap["Red"]);

    G4RotationMatrix* rotmatBeamPipe = new G4RotationMatrix;
//    rotmatGARI->set(0,0,0);
    rotmatBeamPipe->rotateZ(-90*degree);
    rotmatBeamPipe->rotateY(180*degree);
    new G4PVPlacement(rotmatBeamPipe,					// Rotation
                      G4ThreeVector(0,offsetYZ_GARI,beamPipeOffset+offsetYZ_GARI),	// Transformation (Rot&Transl)
                      "physiBeamPipe",			// its logical volume
                      logicBeamPipe,		// its name
                      fMainVolume,		// its physical mother volume
                      false,				// unknown "pMany"; def: false
                      0,			// copy number
                      fCheckOverlaps);		// checkOverlaps

    G4Tubs * solidEncap = new G4Tubs("solidEncap",0,gariMaxRadius,0.1*mm,0,360.*deg);
    G4LogicalVolume* logicEncap = new G4LogicalVolume(solidEncap,          //its solid
                                    G4Material::GetMaterial("TAPE"),         //its material
                                    "logicEncap");        //its name
    logicEncap->SetVisAttributes(ColorMap["invi"]);
//    logicEncap->SetVisAttributes(ColorMap["Red"]);
    G4RotationMatrix* rotmatEncap = new G4RotationMatrix;
    rotmatEncap->rotateZ(-90*degree);
    rotmatEncap->rotateY(180*degree);
    new G4PVPlacement(rotmatEncap,					// Rotation
                      G4ThreeVector(0,offsetYZ_GARI,beamPipeOffset+beamPipeLength/2+0.1*mm+offsetYZ_GARI),	// Transformation (Rot&Transl)
                      "physiBeamPipe",			// its logical volume
                      logicEncap,		// its name
                      fMainVolume,		// its physical mother volume
                      false,				// unknown "pMany"; def: false
                      0,			// copy number
                      fCheckOverlaps);		// checkOverlaps
#endif

#ifdef WCLOVER
    //
    //Place Clover in Moderator
    //
    Clover* clv=new Clover();
    //G4ThreeVector clvPos=G4ThreeVector((-11.6+5.8/2)*cm,0,-5.8/2*cm);
    G4ThreeVector clvPos=G4ThreeVector(40*cm,0.*cm,0.);
    G4RotationMatrix clvRot;
    clvRot.set(0,0,0);
    clvRot.rotateY(90.*deg);
    clvRot.rotateX((-90-90.)*deg);
    //clvRot.invert();
    clv->SetPosition(clvPos);
    clv->SetRotation(clvRot);
    //clv->CreateSolids();
    clv->Placement(0,fMainVolume);

    //2nd one
    Clover* clv2=new Clover();
//    G4ThreeVector clvPos2=G4ThreeVector((11.6/2+0.4)*cm-20*cm,0.*cm,0.);
    G4ThreeVector clvPos2=G4ThreeVector(-15*cm,0.*cm,0.);
    G4RotationMatrix clvRot2;
    clvRot2.set(0,0,0);
    clvRot2.rotateY(90.*deg);
    clvRot2.rotateX(0*deg);
    clv2->SetPosition(clvPos2);
    clv2->SetRotation(clvRot2);
    clv2->Placement(1,fMainVolume);
#endif

    // Place the WLS slab
    if(fWLSslab)
    {
        G4VPhysicalVolume* slab = new LXeWLSSlab(
                    0, G4ThreeVector(0., 0., -fScint_z / 2. - fSlab_z - 1. * cm),
                    fWorld_log, false, 0, this);

        // Surface properties for the WLS slab
        G4OpticalSurface* scintWrap = new G4OpticalSurface("ScintWrap");

        new G4LogicalBorderSurface("ScintWrap", slab, fWorld_phys,
                                   scintWrap);

        scintWrap->SetType(dielectric_metal);
        scintWrap->SetFinish(polished);
        scintWrap->SetModel(glisur);

        std::vector<G4double> pp           = { 2.0 * eV, 3.5 * eV };
        std::vector<G4double> reflectivity = { 1.0, 1.0 };
        std::vector<G4double> efficiency   = { 0.0, 0.0 };

        G4MaterialPropertiesTable* scintWrapProperty =
                new G4MaterialPropertiesTable();

        scintWrapProperty->AddProperty("REFLECTIVITY", pp, reflectivity);
        scintWrapProperty->AddProperty("EFFICIENCY", pp, efficiency);
        scintWrap->SetMaterialPropertiesTable(scintWrapProperty);
    }

    return fWorld_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::ConstructSDandField()
{
    if(!fMainVolume)
        return;

    // PMT SD

    // Add GARI
    LXePMTSD* pmt = fPmt_SD.Get();
    if(!pmt)
    {
        // Created here so it exists as pmts are being placed
        G4cout << "Construction /LXeDet/pmtSD" << G4endl;
        LXePMTSD* pmt_SD = new LXePMTSD("/LXeDet/pmtSD");
        fPmt_SD.Put(pmt_SD);

        pmt_SD->InitPMTs();
        pmt_SD->SetPmtPositions(fMainVolume->GetPmtPositions());
    }
    else
    {
        pmt->InitPMTs();
        pmt->SetPmtPositions(fMainVolume->GetPmtPositions());
    }
    G4SDManager::GetSDMpointer()->AddNewDetector(fPmt_SD.Get());
    // sensitive detector is not actually on the photocathode.
    // processHits gets done manually by the stepping action.
    // It is used to detect when photons hit and get absorbed & detected at the
    // boundary to the photocathode (which doesn't get done by attaching it to a
    // logical volume.
    // It does however need to be attached to something or else it doesn't get
    // reset at the begining of events

    SetSensitiveDetector(fMainVolume->GetLogPhotoCath(), fPmt_SD.Get());
#ifdef WGARI
    SetSensitiveDetector(fgari_det->GetLogPhotoCath(),fPmt_SD.Get());
#endif
    // Scint SD

    if(!fScint_SD.Get())
    {
        G4cout << "Construction /LXeDet/scintSD" << G4endl;
        LXeScintSD* scint_SD = new LXeScintSD("/LXeDet/scintSD");
        fScint_SD.Put(scint_SD);

        scint_SD->InitScints();
        scint_SD->SetScintPositions(fMainVolume->GetScintPositions());
    }
    G4SDManager::GetSDMpointer()->AddNewDetector(fScint_SD.Get());
    SetSensitiveDetector(fMainVolume->GetLogScint(), fScint_SD.Get());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetDimensions(G4ThreeVector dims)
{
    fScint_x = dims[0];
    fScint_y = dims[1];
    fScint_z = dims[2];
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetHousingThickness(G4double d_mtl)
{
    fD_mtl = d_mtl;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetNX(G4int nx)
{
    fNx = nx;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetNY(G4int ny)
{
    fNy = ny;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetNZ(G4int nz)
{
    fNz = nz;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetPMTRadius(G4double outerRadius_pmt)
{
    fOuterRadius_pmt = outerRadius_pmt;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetDefaults()
{
    // Resets to default values
    fD_mtl = 0.0635 * cm;

//    fD_mtl = 10 * cm;

    fExpHall_x = 300 * cm;
    fExpHall_y = 300 * cm;
    fExpHall_z = 300 * cm;


    fScint_x = 17.8 * cm;
    fScint_y = 17.8 * cm;
    fScint_z = 22.6 * cm;

    fNx = 2;
    fNy = 2;
    fNz = 3;

    fOuterRadius_pmt = 2.3 * cm;

    fSphereOn = false;
    fRefl     = 1.0;

    fNfibers      = 15;
    fWLSslab      = false;
    fMainVolumeOn = true;
    fMainVolume   = nullptr;
    fSlab_z       = 2.5 * mm;

    G4UImanager::GetUIpointer()->ApplyCommand(
                "/LXe/detector/scintYieldFactor 1.");

    if(fLXe_mt)
        fLXe_mt->AddConstProperty("SCINTILLATIONYIELD", 12000. / MeV);
    if(fMPTPStyrene)
        fMPTPStyrene->AddConstProperty("SCINTILLATIONYIELD", 10. / keV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetSphereOn(G4bool b)
{
    fSphereOn = b;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetHousingReflectivity(G4double r)
{
    fRefl = r;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetWLSSlabOn(G4bool b)
{
    fWLSslab = b;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetMainVolumeOn(G4bool b)
{
    fMainVolumeOn = b;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetNFibers(G4int n)
{
    fNfibers = n;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetMainScintYield(G4double y)
{
    fLXe_mt->AddConstProperty("SCINTILLATIONYIELD", y / MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetWLSScintYield(G4double y)
{
    fMPTPStyrene->AddConstProperty("SCINTILLATIONYIELD", y / MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetSaveThreshold(G4int save)
{
    // Sets the save threshold for the random number seed. If the number of
    // photons generated in an event is lower than this, then save the seed for
    // this event in a file called run###evt###.rndm

    fSaveThreshold = save;
    G4RunManager::GetRunManager()->SetRandomNumberStore(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
