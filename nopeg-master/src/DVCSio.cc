#define _USE_MATH_DEFINES
#include "DVCSio.h"
#include "DVCSmodel.h"
#include "TRandom3.h"
#include <iostream>
#include <stdlib.h>

DVCSio::DVCSio(char* name){

    file = new TFile(name,"recreate");
    tree = new TTree("TOPEG","Tree with simulated events from TOPEG");
    this->branching();

}

DVCSio::~DVCSio(){
    file->Write();
    delete tree;
    delete file;
}

// Function to create branches of the TFile and connect them to the proper values
void DVCSio::branching(){
    tree->Branch("ievent"    , &ievent           , "ievent/I"             );

    tree->Branch("ElBeam"    , &ElBeam           , "ElBeam/F"             );
    tree->Branch("PIDlBeam"  , &PIDlBeam         , "PIDlBeam/I"           );
    tree->Branch("EhBeam"    , &EhBeam           , "EhBeam/F"             );
    tree->Branch("PIDhBeam"  , &PIDhBeam         , "PIDhBeam/I"           );

    tree->Branch("Q2"        , &Q2               , "Q2/F"                 );
    tree->Branch("W"         , &W                , "W/F"                  );
    tree->Branch("Gamnu"     , &Nu               , "Gamnu/F"              );
    tree->Branch("Xbj"       , &XBj              , "Xbj/F"                );
    tree->Branch("y"         , &y                , "y/F"                  );
    tree->Branch("t"         , &t                , "t/F"                  );
    tree->Branch("phih"      , &phih             , "phih/F"               );

    tree->Branch("Nb_part"   , &Nb               , "Nb_part/I"            );
    tree->Branch("part_id"   , PID               , "part_id[Nb_part]/I"   );
    tree->Branch("part_px"   , px                , "part_px[Nb_part]/F"   );
    tree->Branch("part_py"   , py                , "part_py[Nb_part]/F"   );
    tree->Branch("part_pz"   , pz                , "part_pz[Nb_part]/F"   );
    tree->Branch("part_e"    , E                 , "part_e[Nb_part]/F"    );

    return; 
}

void DVCSio::Fill(double var0,double var1, double var2, double var3){

    DVCSio::InitVars();
    ievent++;
    ElBeam   = DVCSmodel::GetEBene();
    PIDlBeam = 11;
    double PhBeam = DVCSmodel::GetHBmom(); //Beam momentum
    EhBeam   = DVCSmodel::GetHBene();
    PIDhBeam = 1000020040;
    double mas = DVCSmodel::Getmass();
  

    Q2   = var1;
    y    = var0;
    t    = var2;
    phih = var3;

    // Calculate the electron observables
    double ScaEne = (1-y)*ElBeam + PhBeam*Q2/(2*ElBeam*(EhBeam+PhBeam));
    double ThetaE = 2*std::asin(std::sqrt(Q2/4/ElBeam/ScaEne));

    double s = DVCSmodel::GetS();
    XBj = Q2 / y / (s-mas*mas); // This is xA
    double csi = XBj/(2.-XBj);
    double tmin = 4.*mas*mas*csi*csi/(1.-csi*csi);
    XBj *= mas/mN; // This is xB
    Nu  = y/2/mas*(s-mas*mas);
    W   = std::sqrt(mas*mas-Q2+y*(s-mas*mas));

    // Calculate the electron kinematic
    // TODO define this in class
    TRandom *Random = new TRandom3();
    Random->SetSeed((int) ThetaE*100000);

    double PhiEle = Random->Rndm()*2*M_PI;

    double Ek = ScaEne;
    double kx = ScaEne*std::sin(ThetaE)*std::cos(PhiEle);
    double ky = ScaEne*std::sin(ThetaE)*std::sin(PhiEle);
    double kz =-ScaEne*std::cos(ThetaE);

    // Calculate virtual photon
    double Eq = ElBeam - Ek;
    double qx = - kx;
    double qy = - ky;
    double qz = - ElBeam - kz;
    double Pq = std::sqrt(qx*qx+qy*qy+qz*qz);

    // Rotate into COM around x

    double thx = std::atan(qy/(qz+PhBeam));

    double kyp = ky*std::cos(thx)-kz*std::sin(thx);
    double kzp = kz*std::cos(thx)+ky*std::sin(thx);

    double qyp = qy*std::cos(thx)-qz*std::sin(thx);
    double qzp = qz*std::cos(thx)+qy*std::sin(thx);

    double pyp =-PhBeam*std::sin(thx);
    double pzp = PhBeam*std::cos(thx);

    // Rotate into COM around y

    double thy = std::atan(-qx/(qzp+pzp));

    double kxpp = kx*std::cos(thy) + kzp*std::sin(thy);
    double kzpp = kzp*std::cos(thy) - kx*std::sin(thy);

    double qxpp = qx*std::cos(thy) + qzp*std::sin(thy);
    double qzpp = qzp*std::cos(thy) - qx*std::sin(thy);

    double pxpp = pzp*std::sin(thy);
    double pzpp = pzp*std::cos(thy);

    // Boost into COM in z
    double bet = (qzpp+pzpp)/(Eq+EhBeam);
    double gam = 1/std::sqrt(1-bet*bet);

    double kzb = gam*(kzpp-bet*Ek);
    double Ekb = gam*(Ek-bet*kzpp);
    double Pkb = std::sqrt(kxpp*kxpp + kyp*kyp + kzb*kzb);

    double qzb = gam*(qzpp-bet*Eq);
    double Eqb = gam*(Eq-bet*qzpp);
    double Pqb = std::sqrt(qxpp*qxpp + qyp*qyp + qzb*qzb);

    double pzb = gam*(pzpp-bet*EhBeam);
    double Epb = gam*(EhBeam-bet*pzpp);
    double Ppb = std::sqrt(pxpp*pxpp + pyp*pyp + pzb*pzb);

    // Calculate photon kinematic base
    double Eqpb = ((Eqb+Epb)*(Eqb+Epb)-mas*mas)/(Eqb+Epb)/2;
    double Thqb = std::acos(Eqb/Pqb + (Q2-t)/(2*Eqpb*Pqb));

    double tverif0 = Q2 + 2*(Eqb*Eqpb - Pqb*Eqpb*std::cos(Thqb));

    // Longitudinal unit vector
    double lx = qxpp/Pqb;
    double ly = qyp /Pqb;
    double lz = qzb /Pqb;

    // Calculate unit vector for phi_h = 0
    double ux = qyp *(kxpp*qyp -kyp *qxpp)-qzb *(kzb *qxpp-kxpp*qzb );
    double uy = qzb *(kyp *qzb -kzb *qyp )-qxpp*(kxpp*qyp -kyp *qxpp);
    double uz = qxpp*(kzb *qxpp-kxpp*qzb )-qyp *(kyp *qzb -kzb *qyp );
    double Pu = std::sqrt(ux*ux+uy*uy+uz*uz);
    ux /= Pu; uy /= Pu; uz /= Pu;

    // Rotate unit vector by phi_h around q boosted
    double cph = std::cos(phih/180*M_PI);
    double sph = std::sin(phih/180*M_PI);

    double tx = ux*(cph+lx*lx*(1-cph))    + uy*(lx*ly*(1-cph)-lz*sph) + uz*(lx*lz*(1-cph)+ly*sph);
    double ty = ux*(lx*ly*(1-cph)+lz*sph) + uy*(cph+ly*ly*(1-cph))    + uz*(ly*lz*(1-cph)-lx*sph);
    double tz = ux*(lx*lz*(1-cph)-ly*sph) + uy*(ly*lz*(1-cph)+lx*sph) + uz*(cph+lz*lz*(1-cph))   ;

    // Boosted photon kine
    double cth = std::cos(Thqb);
    double sth = std::sin(Thqb);

    double bjkx = (cth*lx + sth*tx)*Eqpb;
    double bjky = (cth*ly + sth*ty)*Eqpb;
    double bjkz = (cth*lz + sth*tz)*Eqpb;

    double tverif1 = Q2 + 2*(Eqb*Eqpb - qxpp*bjkx - qyp *bjky - qzb *bjkz);

    // Now boost back in z
    double jzpp = gam*(bjkz+bet*Eqpb);
    double Ejpp = gam*(Eqpb+bet*bjkz);

    // Rotate back around y
    double jxp = bjkx*std::cos(thy) - jzpp*std::sin(thy);
    double jzp = jzpp*std::cos(thy) + bjkx*std::sin(thy);

    // Rotate back around x
    double jx = jxp;
    double jy = bjky*std::cos(thx)+jzp *std::sin(thx);
    double jz = jzp *std::cos(thx)-bjky*std::sin(thx);

    double tverif2 = - (Eq-Ejpp)*(Eq-Ejpp) + (qx-jx)*(qx-jx) + (qy-jy)*(qy-jy) + (qz-jz)*(qz-jz);

    // Calculate Helium
    double Ppx = qx-jx ;
    double Ppy = qy-jy ;
    double Ppz = qz+PhBeam-jz;
    double Ppp = std::sqrt(Ppx*Ppx + Ppy*Ppy + Ppz*Ppz);
    double Epp = std::sqrt(Ppx*Ppx + Ppy*Ppy + Ppz*Ppz + mas*mas );

    // Set particles kinematic for recording
    Nb = 3;
    PID[0] = 11;
    px [0] = kx;
    py [0] = ky;
    pz [0] = kz;
    E  [0] = Ek;

    PID[1] = 22;
    px [1] = jx;
    py [1] = jy;
    pz [1] = jz;
    E  [1] = Ejpp;

    PID[2] = 1000020040;
    px [2] = Ppx;
    py [2] = Ppy;
    pz [2] = Ppz;
    E  [2] = Epp;

    //All the following is for cross check purposes
    if (0) {
    double delta = Eq+EhBeam-Ejpp-Epp;

    double ax = qyp *kzb -qzb *kyp ;
    double ay = qzb *kxpp-qxpp*kzb ;
    double az = qxpp*kyp -qyp *kxpp;
    double Pa = std::sqrt(ax*ax+ay*ay+az*az);
    ax/=Pa; ay/=Pa; az/=Pa;

    double bx = qyp *bjkz-qzb *bjky;;
    double by = qzb *bjkx-qxpp*bjkz;;
    double bz = qxpp*bjky-qyp *bjkx;;
    double Pb = std::sqrt(bx*bx+by*by+bz*bz);
    bx/=Pb; by/=Pb; bz/=Pb;

    double cx = ay*bz-az*by;
    double cy = az*bx-ax*bz;
    double cz = ax*by-ay*bx;
    double qc = qxpp*cx+qyp *cy+qzb *cz;

    double Sol1 = std::acos(ax*bx+ay*by+az*bz)*180/M_PI;
    if (qc<0) Sol1 = 360-Sol1;


    double dx = qy*kz-qz*ky;
    double dy = qz*kx-qx*kz;
    double dz = qx*ky-qy*kx;
    double Pd = std::sqrt(dx*dx+dy*dy+dz*dz);
    dx/=Pd; dy/=Pd; dz/=Pd;

    double ex = qy*jz-qz*jy;
    double ey = qz*jx-qx*jz;
    double ez = qx*jy-qy*jx;
    double Pe = std::sqrt(ex*ex+ey*ey+ez*ez);
    ex/=Pe; ey/=Pe; ez/=Pe;

    double fx = dy*ez-dz*ey;
    double fy = dz*ex-dx*ez;
    double fz = dx*ey-dy*ex;
    double qf = qxpp*fx+qyp *fy+qzb *fz;

    double Sol2 = std::acos(dx*ex+dy*ey+dz*ez)*180/M_PI;
    if (qf<0) Sol2 = 360-Sol2;


        std::cout << ievent << " event !" << std::endl;
        std::cout << " beam energies " << ElBeam << " " << PhBeam << " " << EhBeam << " " << s <<  std::endl;
        std::cout << " input variables " << Q2 << " " << y << " " << t << " " << phih << std::endl;
        std::cout << " E and angle " << ScaEne << " " << ThetaE << " " << PhiEle <<  " " << t-tmin << std::endl;
        std::cout << " ele obs " << XBj << " " << Nu << " " << W << std::endl;
        std::cout << " ele kine " << px[0] << " " << py[0] << " " << pz[0] << " " << E[0] << std::endl;
        std::cout << " V phot kine " << qx << " " << qy << " " << qz << " " << Eq << " " << Pq << std::endl;
        std::cout << " R phot kine " << px[1] << " " << py[1] << " " << pz[1] << " " << E[1] << std::endl;
        std::cout << " He kine " << px[2] << " " << py[2] << " " << pz[2] << " " << E[2] << std::endl;
        std::cout << " t verif " << t << " " << tverif0 << " " << tverif1 << " " << tverif2 << std::endl;
        std::cout << " Energy and phih verif " << delta << " " << phih << " " << Sol1 << " " << Sol2 << std::endl;

    }

    tree->Fill();
}

void DVCSio::InitVars(){
    PIDlBeam=0;
    PIDhBeam=0;
    ElBeam=0.;
    EhBeam=0.;
    Q2=0.;
    W=0.;
    Nu=0.;
    XBj=0.;
    y=0.;
    t=0.;
    phih=0.;

    Nb=0;
    for (int i=0; i<10; i++){
        PID[10]=0; 
        E[10]=0.;
        px[10]=0.;
        py[10]=0.;
        pz[10]=0.; 
    } 
}
