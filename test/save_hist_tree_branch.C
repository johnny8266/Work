void save_hist_tree_branch()
{
  gBenchmark->Start("hsimple");
  TFile f("ht.root","recreate");
  auto T     = new TTree("T","test");
  auto hpx   = new TH1F("hpx","This is the px distribution",100,-4,4);
  auto hpxpy = new TH2F("hpxpy","py vs px",40,-4,4,40,-4,4);
  auto hprof  = new TProfile("hprof","Profile of pz versus px",100,-4,4,0,20);
  T->Branch("hpx","TH1F",&hpx,32000,0);
  T->Branch("hpxpy","TH2F",&hpxpy,32000,0);
  T->Branch("hprof","TProfile",&hprof,32000,0);
  Float_t px, py, pz;
  for (Int_t i = 0; i < 25000; i++) {
    if (i%1000 == 0) printf("at entry: %d\n",i);
    gRandom->Rannor(px,py);
    pz = px*px + py*py;
    hpx->Fill(px);
    hpxpy->Fill(px,py);
    hprof->Fill(px,pz);

  }
  T->Fill();
  T->Print();
  f.Write();
  gBenchmark->Show("hsimple");
}
