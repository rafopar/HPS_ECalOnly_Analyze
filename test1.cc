#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TObjArray.h>

//================ LCIO ==================
#include "IO/LCReader.h"
#include "IOIMPL/LCFactory.h"
#include <IMPL/LCCollectionVec.h>
#include "EVENT/LCCollection.h"
#include "EVENT/Cluster.h"
#include <EVENT/CalorimeterHit.h>

#include <ECalCluster.hh>

#include <cmath>

#include <iostream>

#include <TApplication.h>

using namespace std;

bool fid_ECal(double, double);
EVENT::CalorimeterHit* getSeedHit(EVENT::Cluster *);
int get_ECal_Raw(double);

// ECalCluster::h_clust1 = new TH2D("h_clust1", "", 42, -283.762, 368.81772, 10, -96.62689, 96.2 );
// ECalCluster::c1 = new TCanvas("c1", "", 750, 750);
// ECalCluster::mark1 = new TMarker();


int main( int argc, char **argv )
{
  cout<<"Kuku"<<endl;
  IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance()->createLCReader();
  TApplication *app1 = new TApplication("app1", NULL, NULL);
  TCanvas *c1 = new TCanvas("c1", "", 2000, 600);
  c1->cd();

  TChain *ch1 = new TChain();

  typedef long long long64;
  const double phot_hole_nom_x = 42.52;

  Int_t n_max_events = 300000000;
  
  if( argc < 2 )
    {
      cout<<"Sorry, Please specify run number "<<endl;
      cout<<"The program is exiting "<<endl;
      exit(1);
    }

  string rnum = argv[1];
  if( argc == 3 )
    {
      n_max_events = atoi(argv[2]);
    }
  

  //ch1->Add(Form("/net/home/rafopar/pumpkin1/hps/MC_Data/recon/trident/2pt2_old_recon/tritrigv1-egsv3*.slcio", rnum.c_str()));
  //ch1->Add(Form("/WD/work/HPS/Data/recon/Comission/hps_00%s.*_recon_*.slcio", rnum.c_str()));
  //ch1->Add("/WD/work/HPS/MC_Data/recon/2pt2/tridents/tritrigv1-egsv3*.slcio");
  //ch1->Add("/home/rafopar/WORKDIR/HPS/Data/test/nullhitsfixed.slcio");

  // ch1->Add(Form("/net/home/rafopar/pumpkin1/hps/Data/engrun/pass0/recon/hps_00%s.*_recon_*.slcio", rnum.c_str()));

  //ch1->Add("/net/home/rafopar/pumpkin1/hps/Data/engrun/pass1/recon/runped+norunped.slcio");

  //=============== Engineering runs ==================================
  //=========pass 0 ============
  //ch1->Add(Form("/net/home/rafopar/pumpkin1/hps/Data/engrun/pass0/recon/hps_00%s.*_recon_*.slcio", rnum.c_str()));
  //=========pass 1 ============
  //ch1->Add(Form("/net/home/rafopar/pumpkin1/hps/Data/engrun/pass1/recon/hps_00%s.*.*.slcio", rnum.c_str()));

  
  //=============== 1p92 tridents simulations ==================================
  ch1->Add(Form("/net/home/rafopar/pumpkin1/hps/MC_Data/recon/trident/1p92/egsv3-triv2-g4v1_s2d6_HPS-ECalCommissioning-v3_20150413v2_pairs*.slcio", rnum.c_str()));  

  //=============== 1p92 Beam bgr simulations ==================================
  //ch1->Add(Form("/net/home/rafopar/pumpkin1/hps/MC_Data/recon/Background/1p92/egsv3-triv2-g4v1_s2d6_HPS-ECalCommissioning-v3_20150413v2_pairs*.slcio", rnum.c_str()));
  
  


  TObjArray *fname_list = ch1->GetListOfFiles();
  int n_files = fname_list->GetEntries();
  for( int i = 0; i < n_files; i++ )
    {
      cout<<(fname_list->At(i))->GetTitle()<<endl;
    }
  
  ECalCluster clust;
  //clust.InitHistCanvas();

  EVENT::LCEvent* ev = NULL;
  
  const double radian = TMath::RadToDeg();
  const int n_max_clust = 20; // assuming number of clusters will not exceed this number
  const int n_fid_raw = 6; // There are 6 ECal raws in Fiducial region
  const int n_ev_time_average = 5000; // Will average 
  const double nanosec = 1.e-9;
  
  double clust_x_[n_max_clust], clust_y_[n_max_clust], clust_z_[n_max_clust];
  double clust_E_[n_max_clust];
  double clust_impact_angle_[n_max_clust];
  double seed_hit_time_[n_max_clust];
  int    clust_size_[n_max_clust];
  double clust_tmin_[n_max_clust];
  double clust_dtmin_[n_max_clust];
  double clust_R_[n_max_clust];
  ECalCluster cluster_ECal_[n_max_clust];

  double clust2_E_[n_max_clust]; // An additional constrain is the energy should be less than a beam energy (let say 1.15 GeV)


  //============================= Here defined variables that there is no requirements for them to be in a fiducial region ======================
  double clust_x_all_[n_max_clust]; double  clust_y_all_[n_max_clust]; double  clust_z_all_[n_max_clust];
  double clust_E_all_[n_max_clust];
  double clust_impact_angle_all_[n_max_clust];
  double seed_hit_time_all_[n_max_clust];
  int    clust_size_all_[n_max_clust];
  double clust_tmin_all_[n_max_clust];
  double clust_dtmin_all_[n_max_clust];
  double clust_R_all_[n_max_clust];
  ECalCluster cluster_ECal_all_[n_max_clust];
  
  TFile *file_out = new TFile(Form("Ecal_only_%s.root", rnum.c_str()), "Recreate");

  TH1D *h_n_clusts1 = new TH1D("h_n_clusts1", "", 30, 0., 15);
  TH1D *h_n_clusts2 = new TH1D("h_n_clusts2", "", 30, 0., 15);
  TH1D *h_n_fid_EC_clust1 = new TH1D("h_n_fid_EC_clust1", "", 30, 0., 15);
  TH1D *h_n2_fid_EC_clust1 = new TH1D("h_n2_fid_EC_clust1", "", 30, 0., 15);
  TH2D *h_clust_yxc1 = new TH2D("h_clust_yxc1", "", 200, -300., 380, 200, -90., 90);
  TH2D *h_clust_yxc2 = new TH2D("h_clust_yxc2", "", 200, -300., 380, 200, -90., 90);
  TH2D *h_clust_yxc3 = new TH2D("h_clust_yxc3", "", 200, -300., 380, 200, -90., 90);
  TH2D *h_clust_yxc4 = new TH2D("h_clust_yxc4", "", 200, -300., 380, 200, -90., 90);
  TH2D *h_clust_yxc5 = new TH2D("h_clust_yxc5", "", 200, -300., 380, 200, -90., 90);
  TH2D *h_clust_yxc6 = new TH2D("h_clust_yxc6", "", 200, -300., 380, 200, -90., 90);
  TH2D *h_clust_yxc7 = new TH2D("h_clust_yxc7", "", 200, -300., 380, 200, -90., 90);
  TH2D *h_clust_yxc8 = new TH2D("h_clust_yxc8", "", 200, -300., 380, 200, -90., 90);
  TH2D *h_clust_yxc9 = new TH2D("h_clust_yxc9", "", 200, -300., 380, 200, -90., 90);
  TH2D *h_clust_yxc10 = new TH2D("h_clust_yxc10", "", 200, -300., 380, 200, -90., 90);
  TH2D *h_clust_yxc11 = new TH2D("h_clust_yxc11", "", 200, -300., 380, 200, -90., 90);
  TH2D *h_clust_yxc12 = new TH2D("h_clust_yxc12", "", 200, -300., 380, 200, -90., 90);
  TH2D *h_clust_yxc13 = new TH2D("h_clust_yxc13", "", 200, -300., 380, 200, -90., 90);
  TH2D *h_clust_yxc14 = new TH2D("h_clust_yxc14", "", 200, -300., 380, 200, -90., 90);
  TH2D *h_clust_yxc15 = new TH2D("h_clust_yxc15", "", 200, -300., 380, 200, -90., 90);
  TH1D *h_clust_size1 = new TH1D("h_clust_size1", "", 30, 0., 15);
  TH1D *h_clust_size2 = new TH1D("h_clust_size2", "", 30, 0., 15);
  TH1D *h_clust_size3 = new TH1D("h_clust_size3", "", 30, 0., 15);
  TH1D *h_clust_size4 = new TH1D("h_clust_size4", "", 30, 0., 15);
  TH1D *h_clust_E1 = new TH1D("h_clust_E1", "", 200, 0., 2.5);
  TH1D *h_clust_E2 = new TH1D("h_clust_E2", "", 200, 0., 2.5);
  TH1D *h_clust_E3 = new TH1D("h_clust_E3", "", 200, 0., 2.5);
  
  TH2D *h_clust_center_yxc1 = new TH2D("h_clust_center_yxc1", "", 200, -300., 380., 200, -90., 90.);
  TH2D *h_clust_center_yxc2 = new TH2D("h_clust_center_yxc2", "", 200, -300., 380., 200, -90., 90.);
  TH2D *h_clust_center_yxc3 = new TH2D("h_clust_center_yxc3", "", 200, -300., 380., 200, -90., 90.);
  TH2D *h_clust_center_yxc4 = new TH2D("h_clust_center_yxc4", "", 200, -300., 380., 200, -90., 90.);
  TH2D *h_clust_center_yxc5 = new TH2D("h_clust_center_yxc5", "", 200, -300., 380., 200, -90., 90.);
  
  TH1D *h_Etot1 = new TH1D("h_Etot1", "", 200, 0., 4.5);
  TH1D *h_Etot2 = new TH1D("h_Etot2", "", 200, 0., 4.5);
  TH1D *h_Etot3 = new TH1D("h_Etot3", "", 200, 0., 4.5);
  TH1D *h_Etot4 = new TH1D("h_Etot4", "", 200, 0., 4.5);
  TH1D *h_Etot5 = new TH1D("h_Etot5", "", 200, 0., 4.5);
  TH1D *h_Etot6 = new TH1D("h_Etot6", "", 200, 0., 4.5);
  TH1D *h_Etot7 = new TH1D("h_Etot7", "", 200, 0., 4.5);
  TH1D *h_Etot8 = new TH1D("h_Etot8", "", 200, 0., 4.5);
  TH1D *h_Etot9 = new TH1D("h_Etot9", "", 200, 0., 4.5);
  TH1D *h_Etot10 = new TH1D("h_Etot10", "", 200, 0., 4.5);
  TH1D *h_Etot11 = new TH1D("h_Etot11", "", 200, 0., 4.5);

  TH1D *h_E_diff1 = new TH1D("h_E_diff1", "", 200, 0., 3.5);
  TH1D *h_E_diff2 = new TH1D("h_E_diff2", "", 200, 0., 3.5);
  TH1D *h_E_diff3 = new TH1D("h_E_diff3", "", 200, 0., 3.5);
  TH1D *h_E_diff4 = new TH1D("h_E_diff4", "", 200, 0., 3.5);
 
  TH2D *h_Ecl_21_1 = new TH2D("h_Ecl_21_1", "", 200, 0., 2.5, 200, 0., 2.5);
  TH2D *h_Ecl_21_2 = new TH2D("h_Ecl_21_2", "", 200, 0., 2.5, 200, 0., 2.5);
  TH2D *h_Ecl_21_3 = new TH2D("h_Ecl_21_3", "", 200, 0., 2.5, 200, 0., 2.5);
  TH2D *h_Ecl_21_4 = new TH2D("h_Ecl_21_4", "", 200, 0., 2.5, 200, 0., 2.5);
  TH2D *h_Ecl_21_5 = new TH2D("h_Ecl_21_5", "", 200, 0., 2.5, 200, 0., 2.5);
  TH2D *h_Ecl_21_6 = new TH2D("h_Ecl_21_6", "", 200, 0., 2.5, 200, 0., 2.5);
  TH2D *h_Ecl_21_7 = new TH2D("h_Ecl_21_7", "", 200, 0., 2.5, 200, 0., 2.5);
  TH2D *h_Ecl_21_8 = new TH2D("h_Ecl_21_8", "", 200, 0., 2.5, 200, 0., 2.5);
  TH2D *h_Ecl_21_9 = new TH2D("h_Ecl_21_9", "", 200, 0., 2.5, 200, 0., 2.5);
  TH2D *h_Ecl_21_10 = new TH2D("h_Ecl_21_10", "", 200, 0., 2.5, 200, 0., 2.5);
  TH2D *h_Ecl_21_11 = new TH2D("h_Ecl_21_11", "", 200, 0., 2.5, 200, 0., 2.5);
  TH2D *h_Ecl_21_12 = new TH2D("h_Ecl_21_12", "", 200, 0., 2.5, 200, 0., 2.5);
  TH2D *h_Ecl_21_13 = new TH2D("h_Ecl_21_13", "", 200, 0., 2.5, 200, 0., 2.5);

  TH2D *h_clust_angle_21_1 = new TH2D("h_clust_angle_21_1", "", 200, 0., 360., 200, 0., 360.);
  TH2D *h_clust_angle_21_2 = new TH2D("h_clust_angle_21_2", "", 200, 0., 360., 200, 0., 360.);
  TH2D *h_clust_angle_21_3 = new TH2D("h_clust_angle_21_3", "", 200, 0., 360., 200, 0., 360.);
  TH2D *h_clust_angle_21_4 = new TH2D("h_clust_angle_21_4", "", 200, 0., 360., 200, 0., 360.);
  TH2D *h_clust_angle_21_5 = new TH2D("h_clust_angle_21_5", "", 200, 0., 360., 200, 0., 360.);
  TH2D *h_clust_angle_21_6 = new TH2D("h_clust_angle_21_6", "", 200, 0., 360., 200, 0., 360.);

  TH1D *h_complanarity1 = new TH1D("h_complanarity_1", "", 200, 0., 360.);
  TH1D *h_complanarity2 = new TH1D("h_complanarity_2", "", 200, 0., 360.);
  TH1D *h_complanarity3 = new TH1D("h_complanarity_3", "", 200, 0., 360.);
  TH1D *h_complanarity4 = new TH1D("h_complanarity_4", "", 200, 0., 360.);
  TH1D *h_complanarity5 = new TH1D("h_complanarity_5", "", 200, 0., 360.);

  TH1D *h_clust_dist1 = new TH1D("h_clust_dist1", "", 200, 0., 600.);
  TH1D *h_clust_dist2 = new TH1D("h_clust_dist2", "", 200, 0., 600.);
  TH1D *h_clust_dist3 = new TH1D("h_clust_dist3", "", 200, 0., 600.);
  TH1D *h_clust_dist4 = new TH1D("h_clust_dist4", "", 200, 0., 600.);
  
  TH2D *h_clust_dist_centre_x1 = new TH2D("h_clust_dist_centre_x1", "", 200, -300., 380., 200, 0., 600.);
  TH2D *h_clust_dist_centre_x2 = new TH2D("h_clust_dist_centre_x2", "", 200, -300., 380., 200, 0., 600.);
  TH2D *h_clust_dist_centre_x3 = new TH2D("h_clust_dist_centre_x3", "", 200, -300., 380., 200, 0., 600.);
  TH2D *h_clust_dist_centre_x4 = new TH2D("h_clust_dist_centre_x4", "", 200, -300., 380., 200, 0., 600.);
  TH2D *h_clust_dist_centre_x5 = new TH2D("h_clust_dist_centre_x5", "", 200, -300., 380., 200, 0., 600.);
  TH2D *h_clust_dist_centre_x6 = new TH2D("h_clust_dist_centre_x6", "", 200, -300., 380., 200, 0., 600.);

  TH2D *h_clust_dx_x1 = new TH2D("h_clust_dx_x1", "", 200, -300, 380, 200, 0., 600);
  TH2D *h_clust_dx_x2 = new TH2D("h_clust_dx_x2", "", 200, -300, 380, 200, 0., 600);
  TH2D *h_clust_dx_x3 = new TH2D("h_clust_dx_x3", "", 200, -300, 380, 200, 0., 600);
  TH2D *h_clust_dx_x4 = new TH2D("h_clust_dx_x4", "", 200, -300, 380, 200, 0., 600);

  TH1D *h_clust_t1 = new TH1D("h_clust_t1", "", 600, 0., 60.);
  
  TH2D *h_clust_t21_1 = new TH2D("h_clust_t21_1", "", 200, 0., 200, 200, 0., 200);
  TH2D *h_clust_t21_2 = new TH2D("h_clust_t21_2", "", 200, 0., 200, 200, 0., 200);
  TH2D *h_clust_t21_3 = new TH2D("h_clust_t21_3", "", 200, 0., 200, 200, 0., 200);
  TH2D *h_clust_t21_4 = new TH2D("h_clust_t21_4", "", 200, 0., 200, 200, 0., 200);
  TH2D *h_clust_t21_5 = new TH2D("h_clust_t21_5", "", 200, 0., 200, 200, 0., 200);
  TH2D *h_clust_t21_6 = new TH2D("h_clust_t21_6", "", 200, 0., 200, 200, 0., 200);

  TH2D *h_dt_min1 = new TH2D("h_dt_min1", "", 200, 0., 200., 200, 0., 30.);
  TH2D *h_dt_min2 = new TH2D("h_dt_min2", "", 200, 0., 200., 200, 0., 30.);
  TH2D *h_dt_min3 = new TH2D("h_dt_min3", "", 200, 0., 200., 200, 0., 30.);
  TH2D *h_dt_min4 = new TH2D("h_dt_min4", "", 200, 0., 200., 200, 0., 30.);
  TH2D *h_dt_min5 = new TH2D("h_dt_min5", "", 200, 0., 200., 200, 0., 30.);
  TH2D *h_dt_min6 = new TH2D("h_dt_min6", "", 200, 0., 200., 200, 0., 30.);

  TH2D *h_Ecl_xc1 = new TH2D("h_Ecl_xc1", "", 200, -300., 380., 200, 0., 2.5);
  TH2D *h_Ecl_xc2 = new TH2D("h_Ecl_xc2", "", 200, -300., 380., 200, 0., 2.5);
  TH2D *h_Ecl_xc3 = new TH2D("h_Ecl_xc3", "", 200, -300., 380., 200, 0., 2.5);
  TH2D *h_Ecl_xc4 = new TH2D("h_Ecl_xc4", "", 200, -300., 380., 200, 0., 2.5);

  TH1D *h_beamline_clust1 = new TH1D("h_beamline_clust1", "", 600, 0., 1.5);

  TH1D *h_x_coplanar1 = new TH1D("h_x_coplanar1", "", 200, -200, 280);
  TH1D *h_x_coplanar2 = new TH1D("h_x_coplanar2", "", 200, -200, 280);
  TH1D *h_x_coplanar3 = new TH1D("h_x_coplanar3", "", 200, -200, 280);
  TH2D *h_x_coplan_Etot1 = new TH2D("h_x_coplan_Etot1", "", 200, 0., 2.5, 200, -200., 280.);
  TH2D *h_x_coplan_Etot2 = new TH2D("h_x_coplan_Etot2", "", 200, 0., 2.5, 200, -200., 280.);

  TH2D *h_clust_E_R1 = new TH2D("h_clust_E_R1", "", 200, 0., 400., 200, 0., 2.5);
  TH2D *h_clust_E_R2 = new TH2D("h_clust_E_R2", "", 200, 0., 400., 200, 0., 2.5);
  TH2D *h_clust_E_R3 = new TH2D("h_clust_E_R3", "", 200, 0., 400., 200, 0., 2.5);
  TH2D *h_clust_E_R4 = new TH2D("h_clust_E_R4", "", 200, 0., 400., 200, 0., 2.5);

  TH1D *h_seed_E1 = new TH1D("h_seed_E1", "", 200, 0., 2.5);
  TH1D *h_seed_E2 = new TH1D("h_seed_E2", "", 200, 0., 2.5);

  TH2D *h_Etot_clust_av_cent_x1 = new TH2D("h_Etot_clust_av_cent_x1", "", 200, -100., 100, 200, 0.2, 2.2);
  TH2D *h_Etot_clust_av_cent_x2 = new TH2D("h_Etot_clust_av_cent_x2", "", 200, -100., 100, 200, 0.2, 2.2);
  TH2D *h_Etot_clust_av_cent_x3 = new TH2D("h_Etot_clust_av_cent_x3", "", 200, -100., 100, 200, 0.2, 2.2);

  TH2D *h_Etot_clust_Xcent_Xav_Average1 = new TH2D("h_Etot_clust_Xcent_Xav_Average1", "", 200, -100., 100., 200, 0.2, 2.2);
  TH2D *h_Etot_clust_Xcent_Xav_Average2 = new TH2D("h_Etot_clust_Xcent_Xav_Average2", "", 200, -100., 100., 200, 0.2, 2.2);

  TH2D *h_clust_coplanX_centre_clust_avX_centre1 = new TH2D("h_clust_coplanX_centre_clust_avX_centre1", "", 200, -100., 100., 200, -100., 100.);
  TH2D *h_clust_coplanX_centre_clust_avX_centre2 = new TH2D("h_clust_coplanX_centre_clust_avX_centre2", "", 200, -100., 100., 200, -100., 100.);
  TH2D *h_clust_coplanX_centre_clust_avX_centre3 = new TH2D("h_clust_coplanX_centre_clust_avX_centre3", "", 200, -100., 100., 200, -100., 100.);
  
  TH2D *h_Clust_mom_perp_paral1 = new TH2D("h_Clust_mom_perp_paral1", "", 200, 0., 25., 200, 0., 25);
  TH2D *h_Clust_mom_perp_paral2 = new TH2D("h_Clust_mom_perp_paral2", "", 200, 0., 25., 200, 0., 25);

  TH2D *h_dE_dt1 = new TH2D("h_dE_dt1", "", 200, 0., 1.3, 200, -30., 30.);
  TH2D *h_dE_dt2 = new TH2D("h_dE_dt2", "", 200, 0., 1.3, 200, -30., 30.);
  TH2D *h_dE_dt3 = new TH2D("h_dE_dt3", "", 200, 0., 1.3, 200, -30., 30.);
  TH2D *h_dE_dt4 = new TH2D("h_dE_dt4", "", 200, 0., 1.3, 200, -30., 30.);

  TH2D *h_Coplan_Etot1 = new TH2D("h_Coplan_Etot1", "", 200, 0., 2.5, 200, 0., 360);
  TH2D *h_Coplan_Etot2 = new TH2D("h_Coplan_Etot2", "", 200, 0., 2.5, 200, 0., 360);
  TH2D *h_Coplan_Etot3 = new TH2D("h_Coplan_Etot3", "", 200, 0., 2.5, 200, 0., 360);
  
  TH2D *h_Coplan_Esum1 = new TH2D("h_Coplan_Esum1", "", 200, 0., 2.5, 200, 0., 360);
  TH2D *h_Coplan_Esum2 = new TH2D("h_Coplan_Esum2", "", 200, 0., 2.5, 200, 0., 360);
  TH2D *h_Coplan_Esum3 = new TH2D("h_Coplan_Esum3", "", 200, 0., 2.5, 200, 0., 360);

  TH2D *h_coplan_Ediff1 = new TH2D("h_coplan_Ediff1", "", 200, 0., 3.5, 200, 0., 360);
  TH2D *h_coplan_Ediff2 = new TH2D("h_coplan_Ediff2", "", 200, 0., 3.5, 200, 0., 360);

  TH2D *h_time_ev_number = new TH2D("h_time_ev_number", "", 10000, 0., 2.e8, 10000, 0., 2.e13);
  TH1D *h_EventRates1 = new TH1D("h_EventRates1", "", 200, 0., 90000.);
  TH1D *h_EventRates2 = new TH1D("h_EventRates2", "", 500, 0., 90000.);
  TH2D *h_Ev_Rates_Ev_number1 = new TH2D("h_Ev_Rates_Ev_number1", "", 10000, 0., 2.e8, 500, 0., 100000.);

  TH2D *h_Ecl_xc1_[n_fid_raw];

  for( int i = 0; i < n_fid_raw; i++ )
    {
      h_Ecl_xc1_[i] = new TH2D(Form("h_Ecl_xc1_%d", i), "", 200, -300., 380., 200, 0., 2.5);
    }


  //=================================== Histograms, where no fiduccial cuts are applied ========================
  TH2D *h_Ecl_21_all_1 = new TH2D("h_Ecl_21_all_1", "", 200, 0., 2.5, 200, 0., 2.5);
  TH2D *h_Ecl_21_all_5 = new TH2D("h_Ecl_21_all_5", "", 200, 0., 2.5, 200, 0., 2.5);
  TH2D *h_Ecl_21_all_6 = new TH2D("h_Ecl_21_all_6", "", 200, 0., 2.5, 200, 0., 2.5);
  TH2D *h_Coplan_Etot_all1 = new TH2D("h_Coplan_Etot_all1", "", 200, 0., 2.5, 200, 0., 360);
  TH2D *h_Coplan_Etot_all2 = new TH2D("h_Coplan_Etot_all2", "", 200, 0., 2.5, 200, 0., 360);
  TH2D *h_Coplan_Etot_all3 = new TH2D("h_Coplan_Etot_all3", "", 200, 0., 2.5, 200, 0., 360);
  TH2D *h_Coplan_Etot_all4 = new TH2D("h_Coplan_Etot_all4", "", 200, 0., 2.5, 200, 0., 360);
  TH2D *h_Coplan_Etot_all5 = new TH2D("h_Coplan_Etot_all5", "", 200, 0., 2.5, 200, 0., 360);
  TH2D *h_Coplan_Etot_all6 = new TH2D("h_Coplan_Etot_all6", "", 200, 0., 2.5, 200, 0., 360);
  TH1D *h_x_coplanar_all1 = new TH1D("h_x_coplanar_all1", "", 200, -200, 280);
  TH1D *h_x_coplanar_all2 = new TH1D("h_x_coplanar_all2", "", 200, -200, 280);
  TH1D *h_x_coplanar_all3 = new TH1D("h_x_coplanar_all3", "", 200, -200, 280);
  TH2D *h_x_coplan_Etot_all1 = new TH2D("h_x_coplan_Etot_all1", "", 200, 0., 2.5, 200, -200., 280.);
  TH2D *h_x_coplan_Etot_all2 = new TH2D("h_x_coplan_Etot_all2", "", 200, 0., 2.5, 200, -200., 280.);
  TH2D *h_x_coplan_Etot_all3 = new TH2D("h_x_coplan_Etot_all3", "", 200, 0., 2.5, 200, -200., 280.);
  TH2D *h_x_coplan_Etot_all4 = new TH2D("h_x_coplan_Etot_all4", "", 200, 0., 2.5, 200, -200., 280.);
  TH2D *h_x_coplan_Etot_all5 = new TH2D("h_x_coplan_Etot_all5", "", 200, 0., 2.5, 200, -200., 280.);
  TH2D *h_x_coplan_Etot_all6 = new TH2D("h_x_coplan_Etot_all6", "", 200, 0., 2.5, 200, -200., 280.);
  TH2D *h_coplan_Ediff_all1 = new TH2D("h_coplan_Ediff_all1", "", 200, 0., 3.5, 200, 0., 360);
  TH2D *h_coplan_Ediff_all2 = new TH2D("h_coplan_Ediff_all2", "", 200, 0., 3.5, 200, 0., 360);
  TH2D *h_coplan_Ediff_all3 = new TH2D("h_coplan_Ediff_all3", "", 200, 0., 3.5, 200, 0., 360);
  TH2D *h_coplan_Ediff_all4 = new TH2D("h_coplan_Ediff_all4", "", 200, 0., 3.5, 200, 0., 360);
  TH2D *h_clust_yxc_all1 = new TH2D("h_clust_yxc_all1", "", 200, -300., 380, 200, -90., 90);
  TH2D *h_clust_yxc_all2 = new TH2D("h_clust_yxc_all2", "", 200, -300., 380, 200, -90., 90);
  TH2D *h_clust_t21_all1 = new TH2D("h_clust_t21_all1", "", 200, 0., 200, 200, 0., 200);
  TH2D *h_clust_t21_all2 = new TH2D("h_clust_t21_all2", "", 200, 0., 200, 200, 0., 200);
  TH2D *h_clust_t21_all3 = new TH2D("h_clust_t21_all3", "", 200, 0., 200, 200, 0., 200);
  TH1D *h_clust21_dt_all1 = new TH1D("h_clust21_dt_all1", "", 600, -150., 150);
  TH1D *h_clust21_dt_all2 = new TH1D("h_clust21_dt_all2", "", 600, -150., 150);
  TH1D *h_clust21_dt_all3 = new TH1D("h_clust21_dt_all3", "", 600, -150., 150);

  TH2D *h_temp1 = new TH2D("h_temp1", "", 200, -300., 380, 200, -90., 90);

  int ind_global = 0;
  long64 prev_time = 0.;
  int prev_hps_ev_num = 0.;
  int ev_average = 0;
  long64 time_begin = 0.;
  // ============== Loop over input files ==================
  for( int ifile = 0; ifile < n_files; ifile++ )
    {
      lcReader->open((fname_list->At(ifile))->GetTitle());
      //      int nev = lcReader->getNumberOfEvents();
  
      cout<<"file "<<(fname_list->At(ifile))->GetTitle()<<endl;
      
      int ev_number = 0;
      while( (ev = lcReader->readNextEvent()) != 0)
	{
	  if( ind_global%50000 == 0 )
	    {
	      cout.flush()<<"Processed "<<ind_global<<"\r";
	    }

	  long64 time_stamp = ev->getTimeStamp();
	  int hps_ev_number = ev->getEventNumber();

	  long64 dt = (time_stamp  - prev_time)*4.;
	  int d_ev_number = hps_ev_number - prev_hps_ev_num;
	  
	  if( ev_average >= n_ev_time_average)
	    {
	      double average_rate = double(ev_average)/(double(time_stamp - time_begin)*nanosec);
	      //cout<<"Average_rate = "<<average_rate<<endl;
	      h_EventRates2->Fill(average_rate);
	      h_Ev_Rates_Ev_number1->Fill(hps_ev_number, average_rate);
	      ev_average = 0;
	      time_begin = time_stamp;
	    }
	  ev_average = ev_average + 1;
	  

	  //cout<<"hps_ev_numer = "<<hps_ev_number<<endl;
	  //cout<<"d_ev_number = "<<d_ev_number<<endl;

	  double Ev_rate = 1.e9*d_ev_number/dt; // 1e9 for converting ns^{-1} into sec^{-1}
	  h_EventRates1->Fill(Ev_rate);
	  prev_time = time_stamp;
	  prev_hps_ev_num = hps_ev_number;
	  //cout<<"time_stamp = "<<time_stamp<<"    hps_ev_number = "<<hps_ev_number<<endl;
	  
	  h_time_ev_number->Fill(hps_ev_number, time_stamp);

	  ind_global = ind_global + 1;
	  const vector<string> *col_names = ev->getCollectionNames();
	  int col_size = col_names->size();
	  
	  //	  if( ev_number == 5 ) {ev_number = ev_number + 1; continue;}
	  
	  //cout<<"Ev number = "<<ev_number<<endl;

	  for( int icol = 0; icol < col_size; icol++ )
	    {
	      const string cur_col_name = col_names->at(icol);
	      //cout<<"Collection is "<<cur_col_name<<endl;
	      //========================= Calorimeter information ============================
	      if( strcmp("EcalClusters", cur_col_name.c_str()) == 0 )
		//if( strcmp("EcalClustersIC", cur_col_name.c_str()) == 0 )
		// if( strcmp("EcalClustersGTP", cur_col_name.c_str()) == 0 )
		  {
		  IMPL::LCCollectionVec *clusters = (IMPL::LCCollectionVec*)ev->getCollection(cur_col_name);
		  int n_clusters = clusters->getNumberOfElements();
		  h_n_clusts1->Fill(n_clusters);
		  
		  int n_clust_in_fid = 0; 
		  double E_tot_in_fid = 0;
		  double E_tot_all = 0;
		  double E_diff = -100.;
		  double complanarity = -100.; // differenc between angles of two clusters
		  double coplanarity_all = -100.; // differenc between angles of two clusters: Note no fiducial cuts are applied

		  int n2_clust_in_fid = 0; 
		  double E2_tot_in_fid = 0;

		  for( int iclust = 0; iclust < n_clusters; iclust++ )
		    {
		      EVENT::Cluster *cur_cluster = (EVENT::Cluster*)clusters->getElementAt(iclust);
		      const float *clust_pos = cur_cluster->getPosition();
		      double clust_x = clust_pos[0];
		      double clust_y = clust_pos[1];
		      double clust_z = clust_pos[2];
		      
		      double clust_E = cur_cluster->getEnergy();
		      
		      double clust_t = (getSeedHit(cur_cluster))->getTime();
		      double seed_E = (getSeedHit(cur_cluster))->getEnergy();
		      
		      E_tot_all = E_tot_all + clust_E;
		      //cout<<"Time = "<<clust_t<<endl;
		      h_clust_yxc1->Fill(clust_x, clust_y);
		      
		      double clust_impact_angle = atan2(clust_y, clust_x - phot_hole_nom_x)*radian;
		      if( clust_impact_angle < 0 ) {clust_impact_angle = clust_impact_angle + 360.;}

		      EVENT::CalorimeterHitVec hits_vec = (EVENT::CalorimeterHitVec)cur_cluster->getCalorimeterHits();
		      int clust_size = hits_vec.size();
		      h_clust_size1->Fill(clust_size);
		      
		      h_clust_E1->Fill(clust_E);
		      h_seed_E1->Fill(seed_E);
		      
		      bool in_fid = fid_ECal(clust_x, clust_y);
		      
		      clust.SetLCIOCluster(cur_cluster);

		      double t_min = clust.Tmin();
		      double t_max = clust.Tmax();
		      double dt = t_max - t_min;

		      h_dt_min1->Fill(t_min, dt);
		      
		      double clust_r = sqrt((clust_x-phot_hole_nom_x)*(clust_x - phot_hole_nom_x) + clust_y*clust_y);

		      clust.Fill_dE_dt(h_dE_dt1);
		      
		      clust_x_all_[iclust] = clust_x;
		      clust_y_all_[iclust] = clust_y;
		      clust_z_all_[iclust] = clust_z;
		      clust_E_all_[iclust] = clust_E;
		      clust_impact_angle_all_[iclust] = clust_impact_angle;
		      seed_hit_time_all_[iclust] = clust_t;

		      if( in_fid && clust_E > 0.1 ) // > 0.1 because with new clustering algorithm the cluster threshold is 100 MeV (Probably :-) )
			{
			  h_clust_yxc2->Fill(clust_x, clust_y);
			  h_clust_size2->Fill(clust_size);
			  h_clust_t1->Fill(clust_t);
			  h_clust_E2->Fill(clust_E);
			  h_seed_E2->Fill(seed_E);					      

			  clust_x_[n_clust_in_fid] = clust_x;
			  clust_y_[n_clust_in_fid] = clust_y;
			  clust_z_[n_clust_in_fid] = clust_z;
			  clust_E_[n_clust_in_fid] = clust_E;
			  seed_hit_time_[n_clust_in_fid] = clust_t;
			  clust_size_[n_clust_in_fid] = clust_size;
			  clust_tmin_[n_clust_in_fid] = t_min;
			  clust_dtmin_[n_clust_in_fid] = dt;
			  cluster_ECal_[n_clust_in_fid] = clust;
			  clust_R_[n_clust_in_fid] = clust_r;

			  E_tot_in_fid = E_tot_in_fid + clust_E;
			  
			  clust_impact_angle_[n_clust_in_fid] = clust_impact_angle;
			  
			  h_dt_min2->Fill(t_min, dt);

			  h_Ecl_xc1->Fill(clust_x, clust_E);
			  
			  int bl_raw = get_ECal_Raw(clust_y);
			  //cout<<"yc = "<<clust_y<<"   bl_raw = "<<bl_raw<<endl;
			  
			  if( seed_E > 0.6*clust_E  )
			    {
			      const float *seed_hit_pos = (getSeedHit(cur_cluster))->getPosition();
			      //(getSeedHit(cur_cluster))->getPosition;
			      double seed_hit_x = seed_hit_pos[0];
			      h_Ecl_xc1_[bl_raw]->Fill(seed_hit_x, clust_E);
			    }
			  
			  clust.Fill_dE_dt(h_dE_dt2);
			  if( clust_E > 1 )
			    {
			      h_clust_yxc3->Fill(clust_x, clust_y);

			      if( clust_E > 1.5 )
				{
				  h_clust_yxc4->Fill(clust_x, clust_y);
				}
			    }
			  
			  h_clust_E_R1->Fill(clust_r, clust_E);

			  n_clust_in_fid = n_clust_in_fid + 1;

			  if( clust_E < 1.15 )
			    {
			      clust2_E_[n2_clust_in_fid] = clust_E;
			      E2_tot_in_fid = E2_tot_in_fid + clust_E;
			      
			      n2_clust_in_fid = n2_clust_in_fid + 1;
			    }
			}
		    }
		  if( n_clusters == 2 )
		    {
		      h_Ecl_21_all_1->Fill(clust_E_all_[0], clust_E_all_[1]);
		      
		      double clust_dt_all = seed_hit_time_all_[0] - seed_hit_time_all_[1];

		      h_clust_t21_all1->Fill(seed_hit_time_all_[0], seed_hit_time_all_[1]);
		      
		      if( seed_hit_time_all_[0] > 15 && seed_hit_time_all_[1] > 15 )
			{
			  h_clust21_dt_all1->Fill(clust_dt_all);
			}
		      
		      if( seed_hit_time_all_[0] > 15 && seed_hit_time_all_[1] > 15 && clust_dt_all > -3. && clust_dt_all < 3. )
			{
			  coplanarity_all = clust_impact_angle_all_[0] - clust_impact_angle_all_[1];
			  if( coplanarity_all < 0 ) {coplanarity_all = coplanarity_all + 360.;}
			  
			  h_Coplan_Etot_all1->Fill(E_tot_all, coplanarity_all);
			  
			  double x_coplan_centr = (clust_x_all_[0] + clust_x_all_[1]*TMath::Abs(clust_y_all_[0]/clust_y_all_[1]))/(1 + TMath::Abs(clust_y_all_[0]/clust_y_all_[1]));
			  h_x_coplanar_all1->Fill(x_coplan_centr);
			  h_x_coplan_Etot_all1->Fill(E_tot_all, x_coplan_centr);
			  double E_diff_all = TMath::Abs(clust_E_all_[0] - clust_E_all_[1]);
			  h_coplan_Ediff_all1->Fill(E_diff_all, coplanarity_all);
			  
			  if( E_diff_all < 0.5 )
			    {
			      h_Coplan_Etot_all3->Fill(E_tot_all, coplanarity_all);
			      h_x_coplan_Etot_all3->Fill(E_tot_all, x_coplan_centr);
			    }
			}
		    }
		  else if( n_clusters > 2 )
		    {
		      for( int ii_clust = 0; ii_clust < n_clusters - 1; ii_clust++ )
			{
			  for( int jj_clust = ii_clust + 1; jj_clust < n_clusters; jj_clust++ )
			    {
			      double E_sum = clust_E_all_[ii_clust] + clust_E_all_[jj_clust];
			      double clust_dt_all = seed_hit_time_all_[ii_clust] - seed_hit_time_all_[jj_clust];
			      double x_coplan_centr = (clust_x_all_[ii_clust] + clust_x_all_[jj_clust]*
						       TMath::Abs(clust_y_all_[ii_clust]/clust_y_all_[jj_clust]))/(1 + TMath::Abs(clust_y_all_[ii_clust]/clust_y_all_[jj_clust]));
			      
			      coplanarity_all = clust_impact_angle_all_[ii_clust] - clust_impact_angle_all_[jj_clust];
			      if( coplanarity_all < 0 ) { coplanarity_all = coplanarity_all + 360.;}
			      
			      h_Coplan_Etot_all2->Fill(E_sum, coplanarity_all);
			      
			      h_Ecl_21_all_5->Fill(clust_E_all_[ii_clust], clust_E_all_[jj_clust]);
			      h_x_coplan_Etot_all2->Fill(E_sum, x_coplan_centr);
			      h_clust21_dt_all2->Fill(clust_dt_all);
			      h_clust_t21_all2->Fill(seed_hit_time_all_[ii_clust], seed_hit_time_all_[jj_clust]);
			      if( seed_hit_time_all_[ii_clust] > 15 && seed_hit_time_all_[jj_clust] > 15 && clust_dt_all > -3. && clust_dt_all < 3. )
				{
				  h_Ecl_21_all_6->Fill(clust_E_all_[ii_clust], clust_E_all_[jj_clust]);
				  
				  h_Coplan_Etot_all5->Fill(E_sum, coplanarity_all);
				  
				  h_x_coplanar_all2->Fill(x_coplan_centr);
				  h_x_coplan_Etot_all5->Fill(E_sum, x_coplan_centr);
				  double E_diff_all = TMath::Abs(clust_E_all_[ii_clust] - clust_E_all_[jj_clust]);
				  h_coplan_Ediff_all2->Fill(E_diff_all, coplanarity_all);
				  
				  //  There are many cluster pairs having sum energy more than the beam energy
				  if( clust_E_all_[ii_clust] + clust_E_all_[jj_clust] > 1.95 )
				    {
				      h_clust21_dt_all3->Fill(clust_dt_all);
				      
				      h_clust_yxc_all2->Fill(clust_x_all_[ii_clust], clust_y_all_[ii_clust]);
				      h_clust_yxc_all2->Fill(clust_x_all_[jj_clust], clust_y_all_[jj_clust]);
				      h_clust_t21_all3->Fill(seed_hit_time_all_[ii_clust], seed_hit_time_all_[jj_clust]);
				    }
				  if( E_diff_all < 0.5 )
				    {
				      h_Coplan_Etot_all4->Fill(E_sum, coplanarity_all);
				      h_x_coplan_Etot_all4->Fill(E_sum, x_coplan_centr);
				    }
				}

			      // There are strange events peaked at coplanarity around 0 and 360
// 			      if(  coplanarity_all < 10 ||  coplanarity_all > 350 )
// 				{
// 				  h_clust_yxc_all1->Fill(clust_x_all_[ii_clust], clust_y_all_[ii_clust]);
// 				  h_clust_yxc_all1->Fill(clust_x_all_[jj_clust], clust_y_all_[jj_clust]);
// 				  h_temp1->Fill(clust_x_all_[ii_clust], clust_y_all_[ii_clust]);
// 				  h_temp1->Fill(clust_x_all_[jj_clust], clust_y_all_[jj_clust]);
// 				  h_temp1->Draw("colz");
// 				  c1->Modified();
// 				  c1->Update();
// 				  cout<<"clust_angle ["<<ii_clust<<"] = "<<clust_impact_angle_all_[ii_clust]<<endl;
// 				  cout<<"clust_angle ["<<jj_clust<<"] = "<<clust_impact_angle_all_[jj_clust]<<endl;
// 				  cin.ignore();
// 				  h_temp1->Reset();
// 				}
			      
			    }
			}
		      
		    }
		  
		  if( n_clust_in_fid > 0 )
		    {
		      h_n_fid_EC_clust1->Fill(n_clust_in_fid);
		      h_Etot1->Fill(E_tot_in_fid);
		      
		      if( n_clust_in_fid == 2 )
			{
			  E_diff = TMath::Abs(clust_E_[0] - clust_E_[1]);
			  h_Etot2->Fill(E_tot_in_fid);
			  h_Ecl_21_1->Fill(clust_E_[0], clust_E_[1]);
			  h_E_diff1->Fill(E_diff);
			  
			  h_clust_t21_1->Fill(seed_hit_time_[0], seed_hit_time_[1]);
			  
			  cluster_ECal_[0].Fill_dE_dt(h_dE_dt3);
			  cluster_ECal_[1].Fill_dE_dt(h_dE_dt3);

			  double clust_dist = sqrt((clust_x_[1] - clust_x_[0])*(clust_x_[1] - clust_x_[0]) + (clust_y_[1] - clust_y_[0])*(clust_y_[1] - clust_y_[0]));
			  double clust_average_cent_x = (clust_x_[0]*clust_E_[0] + clust_x_[1]*clust_E_[1])/(clust_E_[0] + clust_E_[1]);
			  double clust_average_cent_y = (clust_y_[0]*clust_E_[0] + clust_y_[1]*clust_E_[1])/(clust_E_[0] + clust_E_[1]);
			  // double clust_average_cent_x = (clust_x_[0] + clust_x_[1])/2.;
			  // double clust_average_cent_y = (clust_y_[0] + clust_y_[1])/2.;
			  
			  double clust_dX = TMath::Abs(clust_x_[1] - clust_x_[0]);
			  
			  // This is the complanarity center, it is calculated assuming two clusters are exactly coplanar w.r.t. that point
			  double x_complan_centr = (clust_x_[0] + clust_x_[1]*TMath::Abs(clust_y_[0]/clust_y_[1]))/(1 + TMath::Abs(clust_y_[0]/clust_y_[1]));
			  h_x_coplanar1->Fill(x_complan_centr);
			  
			  double x_av_clustCentX_avX = (clust_average_cent_x + x_complan_centr)/2.;
			  h_Etot_clust_Xcent_Xav_Average1->Fill(x_av_clustCentX_avX, E_tot_in_fid);

			  h_clust_angle_21_1->Fill(clust_impact_angle_[0], clust_impact_angle_[1]);
			  
			  h_clust_dist1->Fill(clust_dist);
			  h_clust_center_yxc1->Fill(clust_average_cent_x, clust_average_cent_y);
			  
			  h_clust_dist_centre_x1->Fill(clust_average_cent_x, clust_dist);
			  h_clust_dx_x1->Fill(clust_average_cent_x, clust_dX);
			  
			  complanarity = clust_impact_angle_[0] - clust_impact_angle_[1];
			  if( complanarity < 0 ) {complanarity = complanarity + 360.;}
			  
			  h_complanarity1->Fill(complanarity);
			  
			  h_x_coplan_Etot1->Fill(E_tot_in_fid, x_complan_centr);

			  h_clust_coplanX_centre_clust_avX_centre1->Fill(clust_average_cent_x, x_complan_centr);
			  //======== Calculate cluster moments ===========
			  // The following parameter aa and bb define y = aa + bb*x line, that connects two cluster centers
			  double bb = (clust_y_[1] - clust_y_[0])/(clust_x_[1] - clust_x_[0]);
			  double aa = clust_y_[1] - bb*clust_x_[1];
			  double moment_parallel_0 = cluster_ECal_[0].GetShowerMoment(aa, bb);

			  // Now the line that is perpendicular to the above line (with aa and bb)
			  // Note while 1st line is the same for two clusters, the perpendicular line
			  // is different from cluster to cluster, so two different lines should be calculated
			  double bb_perp_0 = -1./bb;
			  double aa_perp_0 = clust_y_[0] - bb_perp_0*clust_x_[0];
			  double moment_perp_0 = cluster_ECal_[0].GetShowerMoment(aa_perp_0, bb_perp_0);
			  
			  h_Clust_mom_perp_paral1->Fill(moment_parallel_0, moment_perp_0);
			  
			  double bb_perp_1 = bb_perp_0;
			  double moment_parallel_1 = cluster_ECal_[1].GetShowerMoment(aa, bb);
			  double aa_perp_1 = clust_y_[1] - bb_perp_0*clust_x_[1];
			  double moment_perp_1 = cluster_ECal_[1].GetShowerMoment(aa_perp_1, bb_perp_1);
			  h_Clust_mom_perp_paral1->Fill(moment_parallel_1, moment_perp_1);
			  
			  h_Etot_clust_av_cent_x1->Fill(clust_average_cent_x, E_tot_in_fid);
			  
			  h_clust_size3->Fill(clust_size_[0]);
			  h_clust_size3->Fill(clust_size_[1]);
			  
			  h_dt_min3->Fill(clust_tmin_[0], clust_dtmin_[0]);
			  h_dt_min3->Fill(clust_tmin_[1], clust_dtmin_[1]);

			  h_Ecl_xc2->Fill(clust_x_[0], clust_E_[0]);
			  h_Ecl_xc2->Fill(clust_x_[1], clust_E_[1]);
			  
			  h_coplan_Ediff1->Fill(E_diff, complanarity);
			  h_Coplan_Etot1->Fill(E_tot_in_fid, complanarity);
			  
			  if( clust_E_[0] > clust_E_[1] )
			    {
			      h_clust_E_R2->Fill(clust_R_[1], clust_E_[1]);
			    }
			  else { h_clust_E_R2->Fill(clust_R_[0], clust_E_[0]); }
			  
			  
			  if( complanarity > 165. && complanarity < 195. )
			    {
			      h_clust_yxc15->Fill(clust_x_[0], clust_y_[0]);
			      h_clust_yxc15->Fill(clust_x_[1], clust_y_[1]);

			      if( clust_E_[0] > clust_E_[1] )
				{
				  h_clust_E_R3->Fill(clust_R_[1], clust_E_[1]);
				}
			      else { h_clust_E_R3->Fill(clust_R_[0], clust_E_[0]); }
			    }


			  if( E_diff < 0.5 )
			    {
			      h_Coplan_Etot2->Fill(E_tot_in_fid, complanarity);
			      if( E_diff < 0.3 )
				{
				  h_Coplan_Etot3->Fill(E_tot_in_fid, complanarity);
				}
			    }
			  if( E_tot_in_fid < 1.1 )
			    {
			      h_clust_yxc13->Fill(clust_x_[0], clust_y_[0]);
			      h_clust_yxc13->Fill(clust_x_[1], clust_y_[1]);
			      
			      if( complanarity > 160. && complanarity < 200. )
				{
				  h_clust_yxc14->Fill(clust_x_[0], clust_y_[0]);
				  h_clust_yxc14->Fill(clust_x_[1], clust_y_[1]);
				  if( clust_E_[0] > clust_E_[1] )
				    {
				      h_clust_E_R4->Fill(clust_R_[1], clust_E_[1]);
				    }
				  else { h_clust_E_R4->Fill(clust_R_[0], clust_E_[0]); }
				}

			    }
			  
			  // Want to check whether high/low energy component is symmetric in e- and e+ side of Calo or, it is
			  // favorable to one side
			  if( clust_x_[0] > clust_x_[1] ) 
			    {
			      h_Ecl_21_10->Fill(clust_E_[0], clust_E_[1]);
			    }
			  else
			    {
			      h_Ecl_21_10->Fill(clust_E_[1], clust_E_[0]);
			    }
			  
			  if( E_tot_in_fid < 1.1 )
			    {
			      h_complanarity5->Fill(complanarity);
			    }

                          if(  clust_E_[0] < 1.1 && clust_E_[1] < 1.1 )
                            {
                              h_Etot8->Fill(E_tot_in_fid);
			      
			      if( E_tot_in_fid < 1.6 && E_tot_in_fid > 1.3 )
				{
				  h_clust_size4->Fill(clust_size_[0]);
				  h_clust_size4->Fill(clust_size_[1]);
				  h_clust_angle_21_5->Fill(clust_impact_angle_[0], clust_impact_angle_[1]);
				  h_complanarity4->Fill(complanarity);
				  
				  h_x_coplanar2->Fill(x_complan_centr);
				  h_Clust_mom_perp_paral2->Fill(moment_parallel_0, moment_perp_0);
				  h_Clust_mom_perp_paral2->Fill(moment_parallel_1, moment_perp_1);
				  h_clust_dist_centre_x6->Fill(clust_average_cent_x, clust_dist);
				  h_clust_center_yxc5->Fill(clust_average_cent_x, clust_average_cent_y);
				  
				  h_clust_yxc11->Fill(clust_x_[0], clust_y_[0]);
				  h_clust_yxc11->Fill(clust_x_[1], clust_y_[1]);
				  h_clust_t21_5->Fill(seed_hit_time_[0], seed_hit_time_[1]);


				  if( clust_y_[0]*clust_y_[1] < 0 )
				    {
				      h_x_coplanar3->Fill(x_complan_centr);
				      h_Etot_clust_av_cent_x3->Fill(clust_average_cent_x, E_tot_in_fid);
				    }
				  
				}
			      
                            }
                          else
                            {
			      h_Etot9->Fill(clust_E_[0]);
			      h_Etot9->Fill(clust_E_[1]);
                            }
			  
			  
			  if( clust_y_[0] > 0 )
			    {
			      h_clust_yxc5->Fill(clust_x_[1], clust_y_[1]);
			    }
			  else if( clust_y_[1] > 0 )
			    {
			      h_clust_yxc5->Fill(clust_x_[0], clust_y_[0]);
			    }
			  
			  if( clust_E_[0] > 1.2 && clust_E_[1] > 1.2 )
			    {
			      h_clust_t21_6->Fill(seed_hit_time_[0], seed_hit_time_[1]);
			    }
			  
			  if( complanarity > 165. && complanarity < 195. )
			    {
			      h_Ecl_21_9->Fill(clust_E_[0], clust_E_[1]);
			      h_clust_yxc9->Fill(clust_x_[0], clust_y_[0]);
			      h_clust_yxc9->Fill(clust_x_[1], clust_y_[1]);
			    }

			  
			  if( clust_y_[0]*clust_y_[1] < 0 )
			    {
			      h_Ecl_21_2->Fill(clust_E_[0], clust_E_[1]);
			      h_Etot3->Fill(E_tot_in_fid);
			      h_E_diff2->Fill(E_diff);
			      h_clust_t21_2->Fill(seed_hit_time_[0], seed_hit_time_[1]);

			      h_clust_angle_21_2->Fill(clust_impact_angle_[0], clust_impact_angle_[1]);
			      h_complanarity2->Fill(complanarity);
			      
			      h_clust_dist2->Fill(clust_dist);
			      h_clust_center_yxc2->Fill(clust_average_cent_x, clust_average_cent_y);
			      h_clust_dist_centre_x2->Fill(clust_average_cent_x, clust_dist);
			      h_clust_dx_x2->Fill(clust_average_cent_x, clust_dX);
			      
			      h_Ecl_xc3->Fill(clust_x_[0], clust_E_[0]);
			      h_Ecl_xc3->Fill(clust_x_[1], clust_E_[1]);
			      
			      h_Etot_clust_av_cent_x2->Fill(clust_average_cent_x, E_tot_in_fid);
			      h_clust_coplanX_centre_clust_avX_centre2->Fill(clust_average_cent_x, x_complan_centr);
			      h_Etot_clust_Xcent_Xav_Average2->Fill(x_av_clustCentX_avX, E_tot_in_fid);
			      
			      /*
			      if( clust_average_cent_x > -30. && clust_average_cent_x < -20. && E_tot_in_fid < 1.55 && E_tot_in_fid > 1.3 )
				{
				  cluster_ECal_[1].Reset_clust_Hist();
				  cluster_ECal_[0].DrawCluster(cluster_ECal_[1]);
				  //cluster_ECal_[1].DrawCluster();
				  cin.ignore();
				  c1->Modified();
				  c1->Update();
				  c1->Print("Cluster_view.eps");
				}
			      */
			      if( (clust_x_[0] - phot_hole_nom_x)*(clust_x_[1] - phot_hole_nom_x) < 0 )
				{
				  h_Ecl_21_3->Fill(clust_E_[0], clust_E_[1]);
				  h_clust_t21_3->Fill(seed_hit_time_[0], seed_hit_time_[1]);
				  h_Etot4->Fill(E_tot_in_fid);
				  h_clust_yxc6->Fill(clust_x_[0], clust_y_[0]);
				  h_clust_yxc6->Fill(clust_x_[1], clust_y_[1]);
				  h_E_diff3->Fill(E_diff);
				  h_clust_angle_21_3->Fill(clust_impact_angle_[0], clust_impact_angle_[1]);
				  h_complanarity3->Fill(complanarity);
				  
				  h_clust_dist3->Fill(clust_dist);
				  h_clust_center_yxc3->Fill(clust_average_cent_x, clust_average_cent_y);
				  h_clust_dist_centre_x3->Fill(clust_average_cent_x, clust_dist);
				  h_clust_dx_x3->Fill(clust_average_cent_x, clust_dX);

				  if( complanarity > 150. && complanarity < 210. )
				    {
				      h_Ecl_21_5->Fill(clust_E_[0], clust_E_[1]);
				      h_Etot5->Fill(E_tot_in_fid);
				      h_clust_angle_21_4->Fill(clust_impact_angle_[0], clust_impact_angle_[1]);
				      h_clust_dist4->Fill(clust_dist);
				      h_clust_center_yxc4->Fill(clust_average_cent_x, clust_average_cent_y);
				      h_clust_dist_centre_x4->Fill(clust_average_cent_x, clust_dist);
				      h_clust_dx_x4->Fill(clust_average_cent_x, clust_dX);
				      h_clust_yxc7->Fill(clust_x_[0], clust_y_[0]);
				      h_clust_yxc7->Fill(clust_x_[1], clust_y_[1]);
				      h_clust_t21_4->Fill(seed_hit_time_[0], seed_hit_time_[1]);				      

				      if( E_tot_in_fid < 1.7 && E_tot_in_fid > 1.25 )
					{
					  h_clust_yxc8->Fill(clust_x_[0], clust_y_[0]);
					  h_clust_yxc8->Fill(clust_x_[1], clust_y_[1]);
					  h_Ecl_xc4->Fill(clust_x_[0], clust_E_[0]);
					  h_Ecl_xc4->Fill(clust_x_[1], clust_E_[1]);
					}
				      
				      if( E_tot_in_fid < 1.1 )
					{
					  h_clust_dist_centre_x5->Fill(clust_average_cent_x, clust_dist);
					}

				      if( clust_average_cent_x > 0 )
					{
					  h_Ecl_21_6->Fill(clust_E_[0], clust_E_[1]);
					  h_Etot6->Fill(E_tot_in_fid);
					}
				      else
					{
					  h_Ecl_21_7->Fill(clust_E_[0], clust_E_[1]);
					  h_Etot7->Fill(E_tot_in_fid);
					  if( clust_average_cent_x < -120. )
					    {
					      h_Ecl_21_8->Fill(clust_E_[0], clust_E_[1]);
					    }

					}
				    }
				}
			      else
				{
				  h_Ecl_21_4->Fill(clust_E_[0], clust_E_[1]);
				}
			    }
			  else
			    {
			      h_n_clusts2->Fill(n_clusters);
			    }
			}
		      else
			{
			  if( n_clust_in_fid == 3 )
			    {
			      h_Etot11->Fill(E_tot_in_fid);
			      
			      for( int i_clust = 0; i_clust < 3; i_clust++ )
				{
				  h_clust_E3->Fill(clust_E_[i_clust]);
				  h_clust_yxc12->Fill(clust_x_[i_clust], clust_y_[i_clust]);
				}
			    }
			  
			  if( n_clust_in_fid > 2 )
			    {
			      for( int ii_clust = 0; ii_clust < n_clust_in_fid - 1; ii_clust++ )
				{
				  for( int jj_clust = ii_clust+1; jj_clust < n_clust_in_fid; jj_clust++ )
				    {
				      h_Ecl_21_13->Fill(clust_E_[ii_clust], clust_E_[jj_clust]);
				      
				      double E_sum = clust_E_[ii_clust] + clust_E_[jj_clust];
				      
				      E_diff = TMath::Abs(clust_E_[ii_clust] - clust_E_[jj_clust]);

				      complanarity = clust_impact_angle_[ii_clust] - clust_impact_angle_[jj_clust];
				      if( complanarity < 0 ) {complanarity = complanarity + 360.;}
				      
				      h_Coplan_Esum1->Fill(E_sum, complanarity);
				      
				      h_coplan_Ediff2->Fill(E_diff, complanarity);
				      
				      double x_complan_centr = (clust_x_[ii_clust] + clust_x_[jj_clust]*
								TMath::Abs(clust_y_[ii_clust]/clust_y_[jj_clust]))/(1 + TMath::Abs(clust_y_[ii_clust]/clust_y_[jj_clust]));

				      h_x_coplan_Etot2->Fill(E_sum, x_complan_centr);


				      if( E_diff < 0.5 )
					{
					  h_Coplan_Esum2->Fill(E_sum, complanarity);
					  
					  if( E_diff < 0.3 )
					    {
					      h_Coplan_Esum3->Fill(E_sum, complanarity);
					    }
					}
				    }
				}
			    }
			}
		      h_n2_fid_EC_clust1->Fill(n2_clust_in_fid);
		      if( n2_clust_in_fid == 2 )
			{
			  h_Ecl_21_11->Fill(clust2_E_[0], clust2_E_[1]);
			  
			  if( n_clust_in_fid > 2 )
			    {
			      h_Ecl_21_12->Fill(clust2_E_[0], clust2_E_[1]);
			    }
			}
		    }
		}
	      
// 	      else if(strcmp("EcalCalHits", cur_col_name.c_str()) == 0)
// 		{
// 		  IMPL::LCCollectionVec *hits = (IMPL::LCCollectionVec*)ev->getCollection(cur_col_name);
// 		  int n_hits = hits->getNumberOfElements();
		  
// 		  for( int ihit = 0; ihit < n_hits; ihit++ )
// 		    {
// 		      EVENT::CalorimeterHit *cur_hit = (EVENT::CalorimeterHit*)hits->getElementAt(ihit);
		      
// 		      //cout<<"id0 = "<<cur_hit->getCellID0()<<"    id2 = "<<cur_hit->getCellID1()<<endl;
// 		      const float *pos = cur_hit->getPosition();
// 		      //cout<<"hit x = "<<pos[0]<<"   hit y = "<<pos[1]<<endl;
// 		      double bl_x = pos[0];
// 		      double bl_y = pos[1];
		      
// 		      h_clust_yxc10->Fill(bl_x, bl_y);
		      
// 		      double bl_E = cur_hit->getEnergy();
		      
// 		      if( bl_x > -126. && bl_x < 53. && bl_y > -50. && bl_y < 50. )
// 			{
// 			  h_beamline_clust1->Fill(bl_E);
// 			}
		      
// 		    }
		  
// 		}
	      
	    }
	  
	  ev_number = ev_number + 1;
	  if( ind_global > n_max_events ) {break;}
	}
      if( ind_global > n_max_events ) {break;}
    }
  h_n_clusts1->Write();
  h_n_clusts2->Write();
  h_clust_yxc1->Write();
  h_clust_yxc2->Write();
  h_clust_yxc3->Write();
  h_clust_yxc4->Write();
  h_clust_yxc5->Write();
  h_clust_yxc6->Write();
  h_clust_yxc7->Write();
  h_clust_yxc8->Write();
  h_clust_yxc9->Write();
  h_clust_yxc10->Write();
  h_clust_yxc11->Write();
  h_clust_yxc12->Write();
  h_clust_yxc13->Write();
  h_clust_yxc14->Write();
  h_clust_yxc15->Write();
  h_clust_size1->Write();
  h_clust_size2->Write();
  h_clust_size3->Write();
  h_clust_size4->Write();
  h_clust_E1->Write();
  h_clust_E2->Write();
  h_clust_E3->Write();
  h_n_fid_EC_clust1->Write();
  h_Etot1->Write();
  h_Etot2->Write();
  h_Etot3->Write();
  h_Etot4->Write();
  h_Etot5->Write();
  h_Etot6->Write();
  h_Etot7->Write();
  h_Etot8->Write();
  h_Etot9->Write();
  h_Etot11->Write();
  h_Ecl_21_1->Write();
  h_Ecl_21_2->Write();
  h_Ecl_21_3->Write();
  h_Ecl_21_4->Write();
  h_Ecl_21_5->Write();
  h_Ecl_21_6->Write();
  h_Ecl_21_7->Write();
  h_Ecl_21_8->Write();
  h_Ecl_21_9->Write();
  h_E_diff1->Write();
  h_E_diff2->Write();
  h_E_diff3->Write();
  h_clust_angle_21_1->Write();
  h_clust_angle_21_2->Write();
  h_clust_angle_21_3->Write();
  h_clust_angle_21_4->Write();
  h_clust_angle_21_5->Write();
  h_complanarity1->Write();
  h_complanarity2->Write();
  h_complanarity3->Write();
  h_complanarity4->Write();
  h_complanarity5->Write();
  h_clust_dist1->Write();
  h_clust_dist2->Write();
  h_clust_dist3->Write();
  h_clust_dist4->Write();
  h_clust_center_yxc1->Write();
  h_clust_center_yxc2->Write();
  h_clust_center_yxc3->Write();
  h_clust_center_yxc4->Write();
  h_clust_center_yxc5->Write();
  h_clust_dist_centre_x1->Write();
  h_clust_dist_centre_x2->Write();
  h_clust_dist_centre_x3->Write();
  h_clust_dist_centre_x4->Write();
  h_clust_dist_centre_x5->Write();
  h_clust_dist_centre_x6->Write();
  h_clust_dx_x1->Write();
  h_clust_dx_x2->Write();
  h_clust_dx_x3->Write();
  h_clust_dx_x4->Write();
  h_clust_t1->Write();
  h_clust_t21_1->Write();
  h_clust_t21_2->Write();
  h_clust_t21_3->Write();
  h_clust_t21_4->Write();
  h_clust_t21_5->Write();
  h_clust_t21_6->Write();
  h_dt_min1->Write();
  h_dt_min2->Write();
  h_dt_min3->Write();
  h_Ecl_xc1->Write();
  h_Ecl_xc2->Write();
  h_Ecl_xc3->Write();
  h_Ecl_xc4->Write();
  h_beamline_clust1->Write();
  h_x_coplanar1->Write();
  h_x_coplanar2->Write();
  h_x_coplanar3->Write();
  h_clust_E_R1->Write();
  h_clust_E_R2->Write();
  h_clust_E_R3->Write();
  h_clust_E_R4->Write();
  h_seed_E1->Write();
  h_seed_E2->Write();
  h_Clust_mom_perp_paral1->Write();
  h_Clust_mom_perp_paral2->Write();
  h_Etot_clust_av_cent_x1->Write();
  h_Etot_clust_av_cent_x2->Write();
  h_Etot_clust_av_cent_x3->Write();
  h_clust_coplanX_centre_clust_avX_centre1->Write();
  h_clust_coplanX_centre_clust_avX_centre2->Write();
  h_Etot_clust_Xcent_Xav_Average1->Write();
  h_Etot_clust_Xcent_Xav_Average2->Write();
  h_dE_dt1->Write();
  h_dE_dt2->Write();
  h_dE_dt3->Write();
  
  h_Ecl_21_10->Write();
  h_Ecl_21_11->Write();
  h_Ecl_21_12->Write();
  h_Ecl_21_13->Write();

  h_Coplan_Etot1->Write();
  h_Coplan_Etot2->Write();
  h_Coplan_Etot3->Write();
  h_coplan_Ediff1->Write();
  h_coplan_Ediff2->Write();
  h_time_ev_number->Write();
  h_EventRates1->Write();
  h_EventRates2->Write();
  h_Ev_Rates_Ev_number1->Write();
  for( int i = 0; i < n_fid_raw; i++ )
    {
      h_Ecl_xc1_[i]->Write();
    }
  
  h_n2_fid_EC_clust1->Write();
  h_Coplan_Esum1->Write();
  h_Coplan_Esum2->Write();
  h_Coplan_Esum3->Write();
  h_x_coplan_Etot1->Write();
  h_x_coplan_Etot2->Write();

  h_Ecl_21_all_1->Write();
  h_Coplan_Etot_all1->Write();
  h_Ecl_21_all_5->Write();
  h_Ecl_21_all_6->Write();
  h_Coplan_Etot_all2->Write();
  h_Coplan_Etot_all3->Write();
  h_Coplan_Etot_all4->Write();
  h_Coplan_Etot_all5->Write();
  h_x_coplanar_all1->Write();
  h_x_coplanar_all2->Write();
  h_x_coplan_Etot_all1->Write();
  h_x_coplan_Etot_all2->Write();
  h_x_coplan_Etot_all3->Write();
  h_x_coplan_Etot_all4->Write();
  h_x_coplan_Etot_all5->Write();
  h_coplan_Ediff_all1->Write();
  h_coplan_Ediff_all2->Write();
  h_coplan_Ediff_all3->Write();
  h_coplan_Ediff_all4->Write();
  h_clust_yxc_all1->Write();
  h_clust_yxc_all2->Write();
  h_clust_t21_all1->Write();
  h_clust_t21_all2->Write();
  h_clust_t21_all3->Write();
  h_clust21_dt_all1->Write();
  h_clust21_dt_all2->Write();
  h_clust21_dt_all3->Write();
  
  return 0;
}

bool fid_ECal(double x, double y)
{
  bool in_fid = false;
  const double x_edge_low = -262.74;
  const double x_edge_high = 347.7;
  const double y_edge_low = 33.54;
  const double y_edge_high = 75.18;

  const double x_gap_low = -106.66;
  const double x_gap_high = 42.17;
  const double y_gap_high = 47.18;

  y = TMath::Abs(y);

  if( x > x_edge_low && x < x_edge_high && y > y_edge_low && y < y_edge_high )
    {
      if( !(x > x_gap_low && x < x_gap_high && y > y_edge_low && y < y_gap_high) )
        {
          in_fid = true;
        }
    }

  return in_fid;
}

EVENT::CalorimeterHit* getSeedHit(EVENT::Cluster *cluster)
{
  EVENT::CalorimeterHitVec hits_vec = (EVENT::CalorimeterHitVec)cluster->getCalorimeterHits();
  int clust_size = hits_vec.size();
  
  double seedEnergy = 0;
  int seed_index = -1;
  
  for( int ind = 0; ind < clust_size; ind++ )
    {
      //cout<<"hits_vec = "<<hits_vec[ind]<<endl;
      double cur_hit_energy = hits_vec[ind]->getEnergy();
      
      if( cur_hit_energy > seedEnergy )
	{
	  seedEnergy = cur_hit_energy;
	  seed_index = ind;
	}
    }
  
  return hits_vec[seed_index];
}


int get_ECal_Raw(double yc)
{
  const double y_edge_low = 33.54;
  const double y_edge_high = 75.18;
  const double bl_width = (y_edge_high - y_edge_low)/3.;
  
  int bl_raw = 3 + int((yc/TMath::Abs(yc))*(TMath::Abs(yc) - y_edge_low)/bl_width);
  if( yc < 0 ) {bl_raw = bl_raw - 1;}
  return bl_raw;
}
