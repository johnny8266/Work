void pi0_decay()
{
  int N_event = 0;
  const int N_bin = 8;
  const double M_pi0 = 0.134977, zpos = 190., R0 = 29.;
  Double_t Momentum[N_bin] = {2., 4., 7., 10., 14., 18., 22., 30.};
  //  Double_t Momentum[N_bin] = {2.};

  Double_t xpos, ypos, L, R;
  Double_t R_Ctheta, R_Stheta, pi_Sphi, pi_Cphi, phi_plane;

  Double_t E;
  Double_t cosalpha, phi;
  Double_t hit_x, hit_y, good_pos1, good_pos2;
  TVector3 p1_3, p2_3;

  TH1F *angle[N_bin];
  TH1F *E_dis[N_bin];
  TH2F *hit_pos[N_bin];


  TCanvas *c1 = new TCanvas("c1", "c1", 1600, 800);
  c1->Divide(4,2);

  for(int i = 0 ; i < N_bin ; i++)
    {

      angle[i] = new TH1F(Form("angle_%d", i), Form("angle_%d", i), 200, 0., 10.);      
      E_dis[i] = new TH1F(Form("E_dis_%d", i), Form("E_dis_%d", i), (Momentum[i] + 1.)*10, 0., (Momentum[i] + 1.));
      hit_pos[i] = new TH2F(Form("hit_pos_%d", i), Form("hit_pos_%d", i), 60, -60., 60., 60, -60., 60.);

      E = std::sqrt(M_pi0 * M_pi0 + Momentum[i] * Momentum[i]);

      cout << E << endl;
      
      /*
	ofstream myfile;
	myfile.open ("/vol0/pwang-l/Singularity/my_det/sub_crystal/pi0_2g.txt");

	myfile << "MILOU32     Wang, Pu-Kai        IJCLab" << endl;
	myfile << "MILOU EVENT FILE" << endl;
	myfile << "============================================" << endl;
	myfile << "I, ievent, linesnum, weight, genprocess, radcorr,        truex, trueQ2, truey, truet, treuphi, phibelgen, phibelres,       phibelrec" << endl;
	myfile << "============================================" << endl;
	myfile << "I, K(I,1)  K(I,2)  K(I,3)  K(I,4)  K(I,5)             P(I,1)  P(I,2)  P(I,3)  P(I,4)  P(I,5) V(I,1)  V(I,2)  V(I,3)" << endl;
	myfile << "============================================" << endl;
      */

  
      while(1)
	{

	  R = R0 + 5. * gRandom->Rndm();
	  phi_plane = TMath::Pi() * ((2. * gRandom->Rndm()) - 1.);
	  pi_Cphi = TMath::Cos(phi_plane);
	  pi_Sphi = TMath::Sin(phi_plane);
	  L = std::sqrt(R * R + zpos * zpos);
	  R_Stheta = R / L;
	  R_Ctheta = std::sqrt(1. - R_Stheta * R_Stheta);

	  TLorentzVector pi0( (Momentum[i] * R_Stheta * pi_Sphi), (Momentum[i] * R_Stheta * pi_Cphi), (Momentum[i] * R_Ctheta), E);
	  TVector3 boost_V = pi0.BoostVector();
	  
	  cosalpha = (2. * gRandom->Rndm()) - 1.; // Cos alpha uniformly between -1 and 1
	  phi = 2. * TMath::Pi() * (gRandom->Rndm()); // phi between 0 and 2*pi
  
	  TLorentzVector g1(0.5*M_pi0*TMath::Sin(TMath::ACos(cosalpha))*TMath::Cos(phi),
			    0.5*M_pi0*TMath::Sin(TMath::ACos(cosalpha))*TMath::Sin(phi),
			    0.5*M_pi0*cosalpha,
			    0.5*M_pi0);
  
	  TLorentzVector g2(-0.5*M_pi0*TMath::Sin(TMath::ACos(cosalpha))*TMath::Cos(phi),
			    -0.5*M_pi0*TMath::Sin(TMath::ACos(cosalpha))*TMath::Sin(phi),
			    -0.5*M_pi0*cosalpha,
			    0.5*M_pi0);

	  g1.Boost(boost_V); // boost to lab
	  g2.Boost(boost_V); // boost to lab

	  // g1.Print();
	  // g2.Print();
	  // cout << endl << endl;	  

	  hit_x = g1.Px() / g1.Pz() * -1. * L;
	  hit_y = g1.Py() / g1.Pz() * -1. * L;
	  good_pos1 = std::sqrt(hit_x * hit_x + hit_y * hit_y);
	  hit_x = g2.Px() / g2.Pz() * -1. * L;
	  hit_y = g2.Py() / g2.Pz() * -1. * L;
	  good_pos2 = std::sqrt(hit_x * hit_x + hit_y * hit_y);
	  
	  if( ((good_pos1 < 48.) && (good_pos1 > 15.)) && ((good_pos2 < 48.) && (good_pos2 > 15.)) )
	    {
	      N_event++;
	      
	      p1_3.SetXYZ(g1.Px(), g1.Py(), g1.Pz());
	      p2_3.SetXYZ(g2.Px(), g2.Py(), g2.Pz());

	      angle[i]->Fill( (p1_3.Angle(p2_3) / TMath::Pi() * 180.) );
	      E_dis[i]->Fill(g1.E());

	      hit_pos[i]->Fill(hit_x, hit_y);
	      hit_x = g1.Px() / g1.Pz() * -1. * L;
	      hit_y = g1.Py() / g1.Pz() * -1. * L;
	      // hit_x = pi0.Px() / pi0.Pz() * -1. * L;
	      // hit_y = pi0.Py() / pi0.Pz() * -1. * L;
	      hit_pos[i]->Fill(hit_x, hit_y);
	      /*
	      if( g1.E() < g2.E() )
		E_dis->Fill(g1.E());
	      else
		E_dis->Fill(g2.E());
	      */
	      
	      /*
		myfile << "0 ";
		myfile << N_event;
		myfile << " 5 1.0 2 0 0.158214897 5.30490398 4.09338698E-02 1.32313299 0. 0.182215601 1.80325413 0.178929791" << endl;
		myfile << "============================================" << endl;
		myfile << "1 21 11   0 0 0 0.       0.     5.    5.       5.09999983E-04  0.0 0.0 0.0" << endl;
		myfile << "2 21 2212 0 0 0 0.       0.     41.   41.0107  0.9383          0.0 0.0 0.0" << endl;
		myfile << "3 1  11   1 2 0 -2.2035  0.4819 4.53  5.061    5.09999983E-04  0.0 0.0 0.0" << endl;
		myfile << "4 1  22   0 0 0 " << g1.Px() << " " << g1.Py() << " " << g1.Pz() << " " << g1.E() << " 0.0 0.0 0.0 0.0" << endl;
		myfile << "5 1  22   0 0 0 " << g2.Px() << " " << g2.Py() << " " << g2.Pz() << " " << g2.E() << " 0.0 0.0 0.0 0.0" << endl;
		myfile << "6 1  2212 1 2 0 -0.4883  -0.849 31.27 31.3     0.9383          0.0 0.0 0.0" << endl;
		myfile << "=============== Event finished ===============" << endl;
	      */
	    }
	  else
	    continue;

	  if(N_event == 20000)
	    break;

	}

      N_event = 0;
      cout << "Finish " << Momentum[i] << " ......" << endl;

      
      c1->cd(i+1);
      angle[i]->SetStats(0);
      angle[i]->GetXaxis()->SetTitle("#theta^{#circ}");
      angle[i]->SetTitle(Form("#pi0 %.1f GeV, angle between 2 photons", Momentum[i]));
      E_dis[i]->SetStats(0);
      E_dis[i]->GetXaxis()->SetTitle("E [GeV]");
      E_dis[i]->SetTitle(Form("#pi0 %.1f GeV, single photon E", Momentum[i]));
      hit_pos[i]->SetStats(0);
      hit_pos[i]->GetXaxis()->SetTitle("cm");
      hit_pos[i]->GetYaxis()->SetTitle("cm");
      hit_pos[i]->SetTitle(Form("#pi0 %.1f GeV, hit pos of photons", Momentum[i]));
      //     angle[i]->Draw();
      //      E_dis[i]->Draw();
      hit_pos[i]->Draw("colorz");
      
    }



  //  myfile.close();

  
}
