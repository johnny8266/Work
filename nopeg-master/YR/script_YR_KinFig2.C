void script(const char* name, int part, int config)
{
    //Please label the plot with indicating the particle ID. Normalize your counts to an integrated luminosity of 10 fb-1
    TString title;//(#intL=10 fb^{-1} )"); // Title of your plot
    switch(part){
        case 1: title="x_{B} distribution"; break;//(#intL=10 fb^{-1} )"); // Title of your plot
        case 2: title="-t distribution";   break;//(#intL=10 fb^{-1} )"); // Title of your plot
    }

    switch(config){
        case 1: title+=" ( 5x 41)"; break;
        case 2: title+=" (10x110)"; break;
        case 3: title+=" (18x110)"; break;
    }

    TFile *f=new TFile(name);
    TTree *tree=(TTree*)f->Get("TOPEG");

    //Please keep the following binning so that we can compare across different WG
    Double_t limit=3.;
    switch(part){
        case 1: limit=0.02; break; // Title of your plot
        case 2: limit=0.02; break; // Title of your plot
    }

    TH1D *h=new TH1D("h",title.Data(),50,0.,limit);

    //Fill your 2D histogram with Theta (in deg) and Momenta (in GeV) for a given particle ID
    switch(part){
        case 1: tree->Draw("Xbj","Xbj<.03","");         break;
        case 2: tree->Draw("t","t<.4","");         break;
    }

    gStyle->SetOptStat(0);
    gStyle->SetTitleX(2.5);
    gStyle->SetTitleY(2.87);

    gPad->SetRightMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetLogy();

    h->SetTitle("test");
    h->DrawCopy("POL SAME");

    // ---- Output -------

    TFile *fout=new TFile("output_plot.root","recreate");
    h->Write(title.Data());
    gPad->Write(title.Data());
    title+=".gif";
    gPad->Print(title.Data()); 

}
