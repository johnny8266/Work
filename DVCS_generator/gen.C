{
  Double_t M=0.938, t=-0.2, Q2=2.3, phi=180*TMath::DegToRad();

  //Hadronic initial state
  TLorentzVector p(0,0,0,M);
  TLorentzVector q(-0.947635,0,3.60565,3.40559);

  TLorentzVector cms = q + p;
  q.Boost(-cms.BoostVector() );
  p.Boost(-cms.BoostVector() );
  Double_t s=cms.M2();
  cout << "S_var: " << s << endl;
  q.Print(); p.Print();
  
  //We align p and q along the Z-axis
  TVector3 oz(0.,0.,1.);
  TVector3 perpvec=(q.Vect()).Cross(oz);
  Double_t angle=(q.Vect()).Angle(oz);
  cout << "Align angle: " << angle << endl;
  q.Rotate(angle,perpvec.Unit());
  p.Rotate(angle,perpvec.Unit());
  cout << "After rotate: " << endl;
  q.Print(); p.Print();

  //Center-of-mass variables
  Double_t E1cm=(s-Q2-M*M)/(2.*TMath::Sqrt(s));
  Double_t E3cm=(s-M*M)/(2.*TMath::Sqrt(s));
  Double_t p1cm=TMath::Sqrt(E1cm*E1cm+Q2);
  Double_t p3cm=TMath::Sqrt(E3cm*E3cm);
  cout << "Energy and momentum: " << E1cm << " " << E3cm << " " << p1cm << " " << p3cm << endl;
  
  Double_t t0=TMath::Power((-Q2-M*M+M*M)/(2.*TMath::Sqrt(s)),2)-TMath::Power(p1cm-p3cm,2);
  Double_t thetacm=TMath::ASin(2*TMath::Sqrt(-(t-t0)/(4*p1cm*p3cm)));
  cout << "t0: " << t0 << " " << "theta: " << thetacm << " Cos: " << TMath::Cos(thetacm) << endl;
  
  TLorentzVector g(E3cm*TMath::Sin(thetacm),0.,E3cm*TMath::Cos(thetacm),E3cm);
  TLorentzVector pf(q.Px()+p.Px()-g.Px(),q.Py()+p.Py()-g.Py(),q.Pz()+p.Pz()-g.Pz(),q.E()+p.E()-g.E());
  
  //We rotate back the particles
  pf.Rotate(-angle,perpvec.Unit());
  g.Rotate(-angle,perpvec.Unit());
  q.Rotate(-angle,perpvec.Unit());
  
  //We boost them back to the lab
  pf.Boost(cms.BoostVector());
  g.Boost(cms.BoostVector());
  q.Boost(cms.BoostVector());
  g.Print(); pf.Print();

  
  //We rotate around the virtual photon
  //First, we align vectors to the Z-axis
  angle=(q.Vect()).Angle(oz);
  pf.RotateY(angle);
  g.RotateY(angle);
  //Then rotate an angle phi around Z
  pf.RotateZ(phi);
  g.RotateZ(phi);
  //Then rotate back
  pf.RotateY(-angle);
  g.RotateY(-angle);

  //Print out the results
  g.Print();
  pf.Print();
  
}
