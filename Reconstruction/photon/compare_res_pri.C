#include <stdio.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <string.h>
#include "TLatex.h"
#include <TTree.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <TGraphErrors.h>
#include <TFile.h>
#include "TLatex.h"

using namespace std;

double GausM(double *x, double *par)
{
  return par[0] * exp(-0.5 * TMath::Power(((x[0] - par[1]) / par[2]), 2)) + par[3];
}


Double_t Eresolution_fit(double *x, double *par)
{
  return sqrt(TMath::Power(par[0] / sqrt(x[0]), 2) + TMath::Power(par[1], 2));
}


Double_t E_resolu_fit_quardratic_sum(double *x, double *par)
{
  return sqrt( par[0] * par[0] + (par[1] * par[1]) / x[0] + (par[2] * par[2]) / (x[0] * x[0]) );
}


Double_t E_resolu_fit_high(double *x, double *par)
{
  return ( par[0] + par[1] / sqrt(x[0]) + par[2] / x[0]);
}


Double_t E_resolu_fit(double *x, double *par)
{
  return ( par[0] + par[1] / sqrt(x[0]) );
}


void compare_res_pri()
{
  
  std::string path = "../data/";
  //  std::string path = "../data/glass/";
  std::string fileName = path + "outCluster.root";
    
  cout << "************ " << fileName << " ************" << endl;

  
  TFile *f1 = new TFile(fileName.c_str(), "Read");
  TTree *outTree = (TTree *) f1->Get("outTree");

  TTreeReader fReader("outTree", f1);

  TTreeReaderValue<int>        N_cluster = {fReader, "N_cluster"};
  TTreeReaderArray<double>     Cl_seed_energy = {fReader, "Cl_seed_energy"};
  TTreeReaderArray<double>     Cl_seed_x = {fReader, "Cl_seed_x"};
  TTreeReaderArray<double>     Cl_seed_y = {fReader, "Cl_seed_y"};
  TTreeReaderArray<double>     Cl_seed_z = {fReader, "Cl_seed_z"};
  TTreeReaderArray<double>     Cl_x = {fReader, "Cl_x"};
  TTreeReaderArray<double>     Cl_y = {fReader, "Cl_y"};
  TTreeReaderArray<double>     Cl_radius = {fReader, "Cl_radius"};
  TTreeReaderArray<double>     Cl_theta = {fReader, "Cl_theta"};
  TTreeReaderArray<double>     Cl_phi = {fReader, "Cl_phi"};
  TTreeReaderArray<double>     Cl_Energy_tot_simul = {fReader, "Cl_Energy_tot_simul"};
  TTreeReaderArray<double>     Cl_Energy_pe = {fReader, "Cl_Energy_pe"};
  TTreeReaderArray<double>     Cl_size_simul = {fReader, "Cl_size_simul"};
  TTreeReaderArray<double>     Cl_par_a = {fReader, "Cl_par_a"};
  TTreeReaderArray<double>     Cl_x_corr = {fReader, "Cl_x_corr"};
  TTreeReaderArray<double>     Cl_y_corr = {fReader, "Cl_y_corr"};

  TTreeReaderArray<double>     g_pjt_x = {fReader, "g_pjt_x"};
  TTreeReaderArray<double>     g_pjt_y = {fReader, "g_pjt_y"};
  TTreeReaderArray<double>     g_E_all = {fReader, "g_E_all"};

  /*
  TTreeReaderValue<double>     e_hit_emcal_x = {fReader, "e_hit_emcal_x"};
  TTreeReaderValue<double>     e_hit_emcal_y = {fReader, "e_hit_emcal_y"};
  TTreeReaderValue<double>     e_hit_emcal_z = {fReader, "e_hit_emcal_z"};
  TTreeReaderValue<double>     g_hit_emcal_x = {fReader, "g_hit_emcal_x"};
  TTreeReaderValue<double>     g_hit_emcal_y = {fReader, "g_hit_emcal_y"};
  TTreeReaderValue<double>     g_hit_emcal_z = {fReader, "g_hit_emcal_z"};
  TTreeReaderValue<double>     e_pjt_emcal_x = {fReader, "e_pjt_emcal_x"};
  TTreeReaderValue<double>     e_pjt_emcal_y = {fReader, "e_pjt_emcal_y"};
  TTreeReaderValue<double>     g_pjt_emcal_x = {fReader, "g_pjt_emcal_x"};
  TTreeReaderValue<double>     g_pjt_emcal_y = {fReader, "g_pjt_emcal_y"};
  TTreeReaderValue<double>     e_px = {fReader, "e_px"};
  TTreeReaderValue<double>     e_py = {fReader, "e_py"};
  TTreeReaderValue<double>     e_pz = {fReader, "e_pz"};
  TTreeReaderValue<double>     e_E = {fReader, "e_E"};
  TTreeReaderValue<double>     g_px = {fReader, "g_px"};
  TTreeReaderValue<double>     g_py = {fReader, "g_py"};
  TTreeReaderValue<double>     g_pz = {fReader, "g_pz"};
  TTreeReaderValue<double>     g_E = {fReader, "g_E"};
  TTreeReaderValue<int>        e_flag_emcal = {fReader, "e_flag_emcal"};
  TTreeReaderValue<int>        g_flag_emcal = {fReader, "g_flag_emcal"};
  TTreeReaderValue<int>        N_hit_emcal = {fReader, "N_hit_emcal"};
  */

  
  // TH1F
  //  auto diffE_e_res_pri = new TH1F("diffE_e_res_pri", "diffE_e_res_pri", 100, -5., 5.);
  auto diffE_g_res_pri = new TH1F("diffE_g_res_pri", "diffE_g_res_pri", 50, -2., 2.);
  auto diff_posx_g_res_pri = new TH1F("diff_posx_g_res_pri", "diff_posx_g_res_pri", 60, -300., 300.);
  auto diff_posy_g_res_pri = new TH1F("diff_posy_g_res_pri", "diff_posy_g_res_pri", 60, -300., 300.);
  auto diff_posx_cor_g_res_pri = new TH1F("diff_posx_cor_g_res_pri", "diff_posx_cor_g_res_pri", 50, -50., 50.);
  auto diff_posy_cor_g_res_pri = new TH1F("diff_posy_cor_g_res_pri", "diff_posy_cor_g_res_pri", 50, -50., 50.);
  auto hit_cl_size = new TH1F("hit_cl_size", "hit_cl_size", 20, 0., 20.);
  auto Ratio_E_seed_E_recons = new TH1F("Ratio_E_seed_E_recons", "Ratio_E_seed_E_recons", 50, 0., 1.);
  auto E_reso_all_particles = new TH1F("E_reso_all_particles", "E_reso_all_particles", 10, 0., 10.);
  auto E_correction_all_particles = new TH1F("E_correction_all_particles", "E_correction_all_particles", 10, 0., 10.);
  auto difference_E_simul_E_adc = new TH1F("difference_E_simul_E_adc", "difference_E_simul_E_adc", 20, -100., 100.);
 
  
  // TH2F
  auto diffx_g_E = new TH2F("diffx_g_E", "diffx_g_E", 40, -200., 200., 10, 0., 10.);
  auto diffx_g_x = new TH2F("diffx_g_x", "diffx_g_x", 40, 0., 200., 50, -1500., 0.);
  auto g_a_E = new TH2F("g_a_E", "g_a_E", 20, 60., 100., 20, 0., 10.);
  auto g_diffx_x = new TH2F("g_diffx_x", "g_diffx_x", 40, -100., 100., 50, -1500., 0.);
  auto g_diffy_y = new TH2F("g_diffy_y", "g_diffy_y", 40, -100., 100., 50, -500., 500.);
  auto g_diffx_cor_x = new TH2F("g_diffx_cor_x", "g_diffx_cor_x", 40, -100., 100., 50, -1500., 0.);
  auto g_diffy_cor_y = new TH2F("g_diffy_cor_y", "g_diffy_cor_y", 40, -100., 100., 50, -500., 500.);
  auto Ratio_E_seed_E_recons_xpos = new TH2F("Ratio_E_seed_E_recons_xpos", "Ratio_E_seed_E_recons_xpos", 50, 0., 1., 150, -1500., 0.);
  auto Ratio_E_seed_recons_vs_E_recons = new TH2F("Ratio_E_seed_recons_vs_E_recons", "Ratio_E_seed_recons_vs_E_recons", 50, 0., 1., 100, 0., 10.);
  auto delta_E_pos = new TH2F("delta_E_pos", "delta_E_pos", 140, -1400., 1400., 140, -1400., 1400.);
  auto DE_simul_E_adc_vs_E_primary = new TH2F("DE_simul_E_adc_vs_E_primary", "DE_simul_E_adc_vs_E_primary", 20, -100., 100., 100, 0., 10000.);
    

  TLegend *legend[4];

  TH1F* E_diff_per_bin[10];
  for(int i = 0 ; i < 10 ; i++)
    E_diff_per_bin[i] = new TH1F(Form("E_diff_per_bin_%d", i), Form("E_diff_per_bin_%d", i), 35, -0.5, 3.);


  double g_eD, g_poxD, g_poyD, g_poxD_cor, g_poyD_cor;
  double e_eD, e_poxD, e_poyD, e_poxD_cor, e_poyD_cor;
  double Eseed_vs_Et;

  
  //========================================
  // Read the events and fill the plots
  //========================================


  int good_events = 0;
  size_t events_numer = 0;  
  
  while (fReader.Next())
    {
      if(++events_numer > 20000)
	break;
      auto count = static_cast<size_t>(*N_cluster.Get());
	      
      //      double ge = *g_E.Get(), g_pjx = *g_pjt_emcal_x.Get(), g_pjy = *g_pjt_emcal_y.Get();

      if( (count < 4) && (count > 0) )
	{
	  /*
	      cout << "Without correction  "; 
	      for(int i = 0 ; i < count ; i++)
		cout << "[" << Cl_x[i] << ", " << Cl_y[i] << "] || ";
	      cout << endl;
	  
	      cout << "Correction  "; 
	      for(int i = 0 ; i < count ; i++)
		cout << "[" << Cl_x_corr[i] << ", " << Cl_y_corr[i] << "] || ";
	      cout << endl;
	  
	      for(int i = 0 ; i < count ; i++)
		//	    cout << g_E_all[i] << " ";
		cout << "[" << g_pjt_x[i] << ", " << g_pjt_y[i] << "] || ";
	      cout << endl << endl;
	  */
	  
	  for(int i = 0 ; i < count ; i++)
	    {	      
	      //	  cout <<  Cl_seed_energy[i] << " " << Cl_x[i] << " " << Cl_y[i] << " ";
	      double DE_simul_adc = Cl_Energy_tot_simul[i] - Cl_Energy_pe[i];
	      difference_E_simul_E_adc->Fill(DE_simul_adc);
	      
	      for(int j = 0 ; j < count ; j++)
		{
		  Eseed_vs_Et = Cl_seed_energy[i] / Cl_Energy_tot_simul[i];
	      
		  //		  g_eD = g_E_all[j] - (Cl_Energy_tot_simul[i] / 1000.);
		  g_eD = g_E_all[j] - (Cl_Energy_pe[i] / 1000.);
		  g_poxD = g_pjt_x[j] - Cl_x[i];
		  g_poyD = g_pjt_y[j] - Cl_y[i];
		  g_poxD_cor = g_pjt_x[j] - Cl_x_corr[i];
		  g_poyD_cor = g_pjt_y[j] - Cl_y_corr[i];
		  double criteria = TMath::Sqrt(g_poxD_cor * g_poxD_cor + g_poyD_cor * g_poyD_cor);
		  

		  if( criteria < 50. )
		    {

		      // if( ((g_poxD < 30.) && (g_poxD > -30.)) || ((g_poyD < 30.) && (g_poyD > -30.)) )
		      // 	cout << "[" << Cl_x[i] << ", " << Cl_y[i] << "] || [" << g_pjt_x[j] << ", " << g_pjt_y[j] << "]" << endl;

		      DE_simul_E_adc_vs_E_primary->Fill(DE_simul_adc, g_E_all[j]*1000.);
		      diffE_g_res_pri->Fill(g_eD);

		      double rect_bord = 700.;
		      //		      if( ((Cl_x_corr[i] < rect_bord) && (Cl_x_corr[i] > -rect_bord)) && ( (Cl_y_corr[i] < rect_bord) && (Cl_y_corr[i] > -rect_bord) ) )
		      //		      if( ((Cl_x_corr[i] < -200.) && (Cl_x_corr[i] > -400.)) && ( (Cl_y_corr[i] < 400.) && (Cl_y_corr[i] > 200.) ) )

		      if( ((Cl_x_corr[i] < -850.) || (Cl_x_corr[i] > 850.)) || ( (Cl_y_corr[i] < -850.) || (Cl_y_corr[i] > 850.) ) )
		      //		      if( Cl_x_corr[i] < -850. )
			{
			  diff_posx_g_res_pri->Fill(g_poxD);
			  diff_posy_g_res_pri->Fill(g_poyD);
			  diff_posx_cor_g_res_pri->Fill(g_poxD_cor);
			  diff_posy_cor_g_res_pri->Fill(g_poyD_cor);
			}

		      if( (Cl_Energy_pe[i] / 1000.) < 10. )
			{
			  //			  int bin = (Cl_Energy_tot_simul[i] / 1000.);
			  int bin = (Cl_Energy_pe[i] / 1000.);
			  E_diff_per_bin[bin]->Fill(g_eD);
			  good_events++;
			}
		      
		      
		      if( g_eD > 1.5 )
			{
			  //			  if( g_pjt_x[j] < 200. && g_pjt_x[j] > -200. && g_pjt_y[j] < 200. && g_pjt_y[j] > -200. )
			  //			    {
			      cout << "Primary photon energy: " << g_E_all[j] << endl;
			      cout << "Primary photon position: [" << g_pjt_x[j] << ", " << g_pjt_y[j] << "]" << endl;
			      cout << "Reconstructed energy: " << (Cl_Energy_pe[i] / 1000.) << endl;
			      cout << "Reconstructed position: [" << Cl_x_corr[i] << ", " << Cl_y_corr[i] << "]" << endl;
			      cout << "Energy difference: " << g_eD << "GeV || Seed E: " << Cl_seed_energy[i] << "GeV" << endl;
			      cout << "Position difference: [" << g_poxD_cor << ", " << g_poyD_cor << "]" << endl << endl;
			      delta_E_pos->Fill(g_pjt_x[j], g_pjt_y[j]);
			      //			    }
			}
		      //		      break;
		    }
		}


	      /*
	      if( (Cl_x[i] < 0.) && (*g_flag_emcal.Get()) == 1 )
		{
		  g_eD = *g_E.Get() - Cl_Energy_tot_simul[i] / 1000.;
		  g_poxD = *g_pjt_emcal_x.Get() - Cl_x[i];
		  g_poyD = *g_pjt_emcal_y.Get() - Cl_y[i];
		  g_poxD_cor = *g_pjt_emcal_x.Get() - Cl_x_corr[i];
		  g_poyD_cor = *g_pjt_emcal_y.Get() - Cl_y_corr[i];
		  if( (Cl_Energy_tot_simul[i] / 1000.) < 10. )
		    {
		      int bin = (Cl_Energy_tot_simul[i] / 1000.);
		      E_diff_per_bin[bin]->Fill(g_eD);
		    }
		}

	      if( (Cl_x[i] > 0.) && ((*e_flag_emcal.Get()) == 1) )
		{
		  e_eD = *e_E.Get() - Cl_Energy_tot_simul[i] / 1000.;
		  e_poxD = *e_pjt_emcal_x.Get() - Cl_x[i];
		  e_poyD = *e_pjt_emcal_y.Get() - Cl_y[i];
		  e_poxD_cor = *e_pjt_emcal_x.Get() - Cl_x_corr[i];
		  e_poyD_cor = *e_pjt_emcal_y.Get() - Cl_y_corr[i];
		  if( (Cl_Energy_tot_simul[i] / 1000.) < 10. )
		    {
		      int bin = (Cl_Energy_tot_simul[i] / 1000.);
		      E_diff_per_bin[bin]->Fill(e_eD);
		    }
		}
	    	  
	      if( (Cl_x[i] < 0.) && ((*g_flag_emcal.Get()) == 1) )
		{	      
		  Ratio_E_seed_E_recons->Fill(Eseed_vs_Et);
		  diffE_g_res_pri->Fill(g_eD);
		  diff_posx_g_res_pri->Fill(g_poxD);
		  diff_posy_g_res_pri->Fill(g_poyD);
		  diff_posx_cor_g_res_pri->Fill(g_poxD_cor);
		  diff_posy_cor_g_res_pri->Fill(g_poyD_cor);
		      
		  //	      diffx_g_E->Fill(g_poxD, g_eD);
		  //	      diffx_g_x->Fill(g_poyD, g_pjy);
		  g_a_E->Fill(Cl_par_a[i], (Cl_Energy_tot_simul[i] / 1000.));
		  g_diffx_x->Fill(g_poxD, g_pjx);
		  g_diffy_y->Fill(g_poyD, g_pjy);
		  g_diffx_cor_x->Fill(g_poxD_cor, g_pjx);
		  g_diffy_cor_y->Fill(g_poyD_cor, g_pjy);
		  Ratio_E_seed_E_recons_xpos->Fill(Eseed_vs_Et, g_pjx);
		  Ratio_E_seed_recons_vs_E_recons->Fill(Eseed_vs_Et, (Cl_Energy_tot_simul[i] / 1000.));
		}// photon hit only 
	      */
	    }// loop the N cluster
	}// N cluster < 3
    }


  cout << "Good reconstructed events: " << good_events << endl;
  

  //============================
  // Draw and Fit the results
  //============================
  

  /*
  auto c1 = new TCanvas("c1", "c1", 1000, 1000);
  c1->Divide(2,2);
  c1->cd(1);
  diff_posx_g_res_pri->SetStats(0);
  diff_posx_g_res_pri->GetXaxis()->SetTitle("[mm]");
  diff_posx_g_res_pri->Draw();
  c1->cd(2);
  diff_posy_g_res_pri->SetStats(0);
  diff_posy_g_res_pri->GetXaxis()->SetTitle("[mm]");
  diff_posy_g_res_pri->Draw();
  c1->cd(3);
  g_diffx_x->SetStats(0);
  g_diffx_x->GetXaxis()->SetTitle("[mm]");
  g_diffx_x->Draw("colorz");
  c1->cd(4);
  g_diffy_y->SetStats(0);
  g_diffy_y->GetXaxis()->SetTitle("[mm]");
  g_diffy_y->Draw("colorz");
  
  */

  
  auto *fun_e_res = new TF1("fun_e_res", GausM, -0.5, 2., 4);
  auto c10 = new TCanvas("c10", "c10", 1200, 900);
  c10->Divide(4,3);
  for(int i = 0 ; i < 10 ; i++)
    {
      c10->cd(i+1);
      //      cout << E_diff_per_bin[i]->GetEntries() << endl;
      double up_lim = E_diff_per_bin[i]->GetEntries();
      double plot_mean = E_diff_per_bin[i]->GetMean();
      double plot_std = E_diff_per_bin[i]->GetStdDev();
      fun_e_res->SetParLimits(0, 0.1, up_lim);
      fun_e_res->SetParLimits(1, -0.1, 1.5);
      fun_e_res->SetParLimits(2, 0.01, 3.);
      fun_e_res->SetParLimits(3, 0.1, 50.);
      if( i <= 2 )
	E_diff_per_bin[i]->Fit("fun_e_res", "", "R", ( plot_mean - 3. * plot_std), (plot_mean + plot_std) );
      else
	E_diff_per_bin[i]->Fit("fun_e_res", "", "R", ( plot_mean - 3. * plot_std), plot_mean);
      //      E_diff_per_bin[i]->SetStats(0);
      E_diff_per_bin[i]->GetXaxis()->SetTitle("Delta E [GeV]");
      E_diff_per_bin[i]->SetTitle(Form("Energy %d ~ %d GeV", i, i+1));
      E_diff_per_bin[i]->Draw();
      cout << fun_e_res->GetChisquare() / fun_e_res->GetNDF() << endl;
      
      double e_resolu = (fun_e_res->GetParameter(2)) / (i + 0.5) * 100.;  // x100% 

      if( (i >= 1) )
	E_reso_all_particles->SetBinContent( i+1, e_resolu );
      //      E_reso_all_particles->SetBinError( i+1, (fun_e_res->GetParError(2)) );

      double e_correction = fun_e_res->GetParameter(1);
      E_correction_all_particles->SetBinContent( i+1, e_correction);

    }

  //  auto *fun_recons_e_res = new TF1("fun_recons_e_res", E_resolu_fit_high, 0., 10., 3);
  auto *fun_recons_e_res = new TF1("E_resolu_fit_quardratic_sum", E_resolu_fit, 0., 10., 3);
  fun_recons_e_res->SetParLimits(0, 0.1, 10.);
  fun_recons_e_res->SetParLimits(1, 0.1, 10.);
  fun_recons_e_res->SetParLimits(2, 0.1, 10.);
  auto c9 = new TCanvas("c9", "c9", 1200, 600);
  c9->Divide(2,1);
  c9->cd(1);
  E_reso_all_particles->SetStats(0);
  E_reso_all_particles->SetMarkerStyle(2);
  E_reso_all_particles->SetMarkerSize(1);
  E_reso_all_particles->GetXaxis()->SetTitle("E [GeV]");
  E_reso_all_particles->GetYaxis()->SetTitle("E_sig / E [%]");
  E_reso_all_particles->Fit("E_resolu_fit_quardratic_sum", "", "R", 0., 10.);
  E_reso_all_particles->Draw("p");
  cout << fun_recons_e_res->GetChisquare() / fun_recons_e_res->GetNDF() << endl;
  legend[0] = new TLegend(0.5, 0.6, 0.75, 0.8);
  legend[0]->SetBorderSize(0);
  legend[0]->AddEntry((TObject*)0, Form("Par A:  %.3f", (fun_recons_e_res->GetParameter(0)) ), "");
  legend[0]->AddEntry((TObject*)0, Form("Par B:  %.3f", (fun_recons_e_res->GetParameter(1)) ), "");
  legend[0]->AddEntry((TObject*)0, Form("Par C:  %.3f", (fun_recons_e_res->GetParameter(2)) ), "");
  legend[0]->Draw("same");
        
  c9->cd(2);
  E_correction_all_particles->SetStats(0);
  E_correction_all_particles->SetMarkerStyle(2);
  E_correction_all_particles->SetMarkerSize(1);
  E_correction_all_particles->GetXaxis()->SetTitle("E [GeV]");
  E_correction_all_particles->GetYaxis()->SetTitle("delta E [GeV]");
  E_correction_all_particles->Draw("p");
  

  
  /*  
  auto c3 = new TCanvas("c3", "c3", 1000, 1000);
  c3->Divide(2,2);
  c3->cd(1);
  diffE_g_res_pri->Draw();
  auto *fun_e = new TF1("fun_e", GausM, 0., 0.4, 4);
  fun_e->SetParLimits(0, 5., 2000.);
  fun_e->SetParLimits(1, 0.05, 0.3);
  fun_e->SetParLimits(2, 0.001, 1.);
  fun_e->SetParLimits(3, 0.1, 10.);
  diffE_g_res_pri->Fit("fun_e", "R");
  diffE_g_res_pri->SetStats(0);
  diffE_g_res_pri->GetXaxis()->SetTitle("difference [GeV]");
  //  fun_e->Draw("n");
  legend[0] = new TLegend(0.12, 0.6, 0.48, 0.8);
  legend[0]->SetBorderSize(0);
  legend[0]->AddEntry((TObject*)0, Form("Resolution: %.3fGeV", (fun_e->GetParameter(2)) ), "");
  legend[0]->AddEntry((TObject*)0, Form("Difference: %.3fGeV", (fun_e->GetParameter(1)) ), "");
  legend[0]->AddEntry((TObject*)0, Form("Chi / NDF: %.3f", (fun_e->GetChisquare() / fun_e->GetNDF()) ), "");
  legend[0]->Draw("same");

  c3->cd(2);
  Ratio_E_seed_E_recons->SetStats(0);
  Ratio_E_seed_E_recons->GetXaxis()->SetTitle("E_seed / E_recons");
  Ratio_E_seed_E_recons->Draw();

  c3->cd(3);
  Ratio_E_seed_E_recons_xpos->SetStats(0);
  Ratio_E_seed_E_recons_xpos->GetXaxis()->SetTitle("E_seed / E_recons");
  Ratio_E_seed_E_recons_xpos->Draw("colorz");

  c3->cd(4);
  Ratio_E_seed_recons_vs_E_recons->SetStats(0);
  Ratio_E_seed_recons_vs_E_recons->GetXaxis()->SetTitle("E_seed / E_recons");
  Ratio_E_seed_recons_vs_E_recons->Draw("colorz");
  */
  
  
  auto c4 = new TCanvas("c4", "c4", 1000, 1000);
  c4->Divide(2,2);

  c4->cd(1);
  auto *fun_x = new TF1("fun_x", GausM, -20., 30., 4);
  fun_x->SetParLimits(0, 5., 20000.);
  fun_x->SetParLimits(1, 5., 20.);
  fun_x->SetParLimits(2, 1., 10.);
  fun_x->SetParLimits(3, 0.1, 100.);
  //  diff_posx_g_res_pri->Fit("fun_x");
  //  cout << "Before correction: " << fun_x->GetChisquare() / fun_x->GetNDF() << endl;
  diff_posx_g_res_pri->SetStats(0);
  diff_posx_g_res_pri->GetXaxis()->SetTitle("difference [mm]");
  diff_posx_g_res_pri->Draw();
  
  // legend[0] = new TLegend(0.12, 0.6, 0.45, 0.8);
  // legend[0]->SetBorderSize(0);
  // legend[0]->AddEntry((TObject*)0, Form("Resolution: %.3fmm", (fun_x->GetParameter(2)) ), "");
  // legend[0]->AddEntry((TObject*)0, Form("Difference: %.3fmm", (fun_x->GetParameter(1)) ), "");
  // legend[0]->AddEntry((TObject*)0, Form("Chi / NDF: %.3f", (fun_x->GetChisquare() / fun_x->GetNDF()) ), "");
  // legend[0]->Draw("same");
  
  
  c4->cd(2);
  auto *fun_y = new TF1("fun_y", GausM, -40., 40., 4);
  fun_y->SetParLimits(0, 10., 20000.);
  fun_y->SetParLimits(1, -3., 3.);
  fun_y->SetParLimits(2, 2., 15.);
  fun_y->SetParLimits(3, 0.1, 100.);
  //  diff_posy_g_res_pri->Fit("fun_y");
  //  cout << "Before correction: " << fun_y->GetChisquare() / fun_y->GetNDF() << endl;
  diff_posy_g_res_pri->SetStats(0);
  diff_posy_g_res_pri->GetXaxis()->SetTitle("difference [mm]");
  diff_posy_g_res_pri->Draw();
  
  // legend[0] = new TLegend(0.12, 0.6, 0.45, 0.8);
  // legend[0]->SetBorderSize(0);
  // legend[0]->AddEntry((TObject*)0, Form("Resolution: %.3fmm", (fun_y->GetParameter(2)) ), "");
  // legend[0]->AddEntry((TObject*)0, Form("Difference: %.3fmm", (fun_y->GetParameter(1)) ), "");
  // legend[0]->AddEntry((TObject*)0, Form("Chi / NDF: %.3f", (fun_y->GetChisquare() / fun_y->GetNDF()) ), "");
  // legend[0]->Draw("same");
  
  
  c4->cd(3);
  auto *fun_x_cor = new TF1("fun_x_cor", GausM, -40., 40., 4);
  fun_x_cor->SetParLimits(0, 5., 20000.);
  fun_x_cor->SetParLimits(1, -3., 3.);
  fun_x_cor->SetParLimits(2, 2., 15.);
  fun_x_cor->SetParLimits(3, 0.1, 100.);
  diff_posx_cor_g_res_pri->Fit("fun_x_cor", "R", "", -40., 40.);
  //  cout << "After correction: " << fun_x_cor->GetChisquare() / fun_x_cor->GetNDF() << endl;
  diff_posx_cor_g_res_pri->SetStats(0);
  diff_posx_cor_g_res_pri->GetXaxis()->SetTitle("difference [mm]");
  diff_posx_cor_g_res_pri->Draw();
  legend[0] = new TLegend(0.12, 0.6, 0.4, 0.8);
  legend[0]->SetBorderSize(0);
  legend[0]->AddEntry((TObject*)0, Form("Resolution: %.3fmm", (fun_x_cor->GetParameter(2)) ), "");
  legend[0]->AddEntry((TObject*)0, Form("Difference: %.3fmm", (fun_x_cor->GetParameter(1)) ), "");
  legend[0]->AddEntry((TObject*)0, Form("Chi / NDF: %.3f", (fun_x_cor->GetChisquare() / fun_x_cor->GetNDF()) ), "");
  legend[0]->Draw("same");

  
  c4->cd(4);
  diff_posy_cor_g_res_pri->Fit("fun_y", "R", "", -40., 40.);
  cout << "After correction: " << fun_y->GetChisquare() / fun_y->GetNDF() << endl;
  diff_posy_cor_g_res_pri->SetStats(0);
  diff_posy_cor_g_res_pri->GetXaxis()->SetTitle("difference [mm]");
  diff_posy_cor_g_res_pri->Draw();
  legend[0] = new TLegend(0.12, 0.6, 0.4, 0.8);
  legend[0]->SetBorderSize(0);
  legend[0]->AddEntry((TObject*)0, Form("Resolution: %.3fmm", (fun_y->GetParameter(2)) ), "");
  legend[0]->AddEntry((TObject*)0, Form("Difference: %.3fmm", (fun_y->GetParameter(1)) ), "");
  legend[0]->AddEntry((TObject*)0, Form("Chi / NDF: %.3f", (fun_y->GetChisquare() / fun_y->GetNDF()) ), "");
  legend[0]->Draw("same");


  auto c5 = new TCanvas("c5", "c5", 800, 400);
  c5->Divide(2,1);
  c5->cd(1);
  difference_E_simul_E_adc->SetStats(0);
  difference_E_simul_E_adc->GetXaxis()->SetTitle("[MeV]");
  difference_E_simul_E_adc->Draw();
  c5->cd(2);
  DE_simul_E_adc_vs_E_primary->SetStats(0);
  DE_simul_E_adc_vs_E_primary->GetXaxis()->SetTitle("[MeV]");
  DE_simul_E_adc_vs_E_primary->GetYaxis()->SetTitle("[MeV]");
  DE_simul_E_adc_vs_E_primary->Draw("colorz");
  // delta_E_pos->Draw("colorz");
  

}
