#ifndef TauHelper_h
#define TauHelper_h


#include <vector>
#include "TLorentzVector.h"
#include "TComplex.h"
#include <string.h>
#include <iostream>
#include <string.h>
using namespace std;


class TauHelper {

 public:
  TauHelper();
  //  TauHelper();
  ~TauHelper();
  void Configure(TLorentzVector OSPion, TLorentzVector SSPion1, TLorentzVector SSPion2,TLorentzVector TauA1, TLorentzVector TauMu );

  void Initialize(TLorentzVector t, TLorentzVector mu);
  bool OmegaIsValid(){return isValid_;}
 


  void SetParametersReco(TLorentzVector tau, TLorentzVector mu );
  void SetBoost(TLorentzVector TauA1, TLorentzVector TauMu );
  TLorentzVector Boost(TLorentzVector pB, TLorentzVector frame);

  double costheta();
  double costheta1();
  double costheta2();
  float CosBeta();
  double CosBeta1();
  std::vector<float> Sin2Cos2Gamma(TLorentzVector p1,TLorentzVector p2, TLorentzVector p3);
  float CosPsi();
  float CosPsi1();


  double I_squared_cospsi(double limit);
  double Integral_squared_cospsi();
  double I_cospsi(double limit);
  double Integral_cospsi();
  double I_squared_cospsi_costheta(double limit);
  double Integral_squared_cospsi_costheta();
  double I_cospsi_costheta(double limit);
  double Integral_cospsi_costheta();
  double I_sinpsi_sintheta(double limit);
  double Integral_sinpsi_sintheta();
  double I_sin2psi_sintheta(double limit);
  double Integral_sin2psi_sintheta();

  double getA1omegaIntegratedOverTheta();



  double lambda(double x, double y, double z);
  float Scalar(TLorentzVector p1, TLorentzVector p2);


  //--------------------------- Hadronic current ---------------------
  /* float WA(TLorentzVector s1,TLorentzVector s2,TLorentzVector s3,float QQ); */
  /* float WC(TLorentzVector s1,TLorentzVector s2,TLorentzVector s3,float QQ); */
  /* float WD(TLorentzVector s1,TLorentzVector s2,TLorentzVector s3,float QQ); */
  /* float WE(TLorentzVector s1,TLorentzVector s2,TLorentzVector s3,float QQ); */
  TComplex  F(float s1, float s2,float QQ);
  float VV1(float SS1 ,float SS2, float SS3, float QQ);
  float VV2(float SS1 ,float SS2, float SS3, float QQ);
  float V1V2(float SS1 ,float SS2, float SS3, float QQ);
  float h0(float SS1 ,float SS2, float SS3, float QQ);
  float h(float SS1 ,float SS2, float SS3, float QQ);
  TComplex  BWa1(float QQ);
  TComplex  BWrho(float QQ);
  TComplex  BWrhoPrime(float QQ);
  float GammaA1(float QQ);
  float gForGammaA1(float QQ);
  float GammaRho(float QQ);
  float  GammaRhoPrime(float QQ);





  double GetOmegaA1();

  TLorentzVector Get1TauSolution(){return TauMu1_;}
  TLorentzVector Get2TauSolution(){return TauMu2_;}
  bool BothSidesAreOK(){return Flag_;}

  double VV1();
  double VV2();
  double V1V2();
  double h0();
  double h();
  double WA();
  double WC();
  double WD();
  double WE();
  TComplex  BreitWigner(double Q, string type="rho");
  TComplex  BRho(double Q);
  TComplex F1();
  TComplex F2();
  TComplex F4();
  TComplex Conjugate(TComplex a);
  TVector3 Rotate(TVector3 LVec, TVector3 Rot);
  double  Widths(double Q, string type="rho");
  double ppi(double QQ);
  double ga1(double  Q);
  double Mass(string type="rho");


  TComplex BWIGML(double S, double M,  double G, double m1, double m2, int L);
  TComplex FPIKM(double W, double XM1, double XM2);
  TComplex F3PI(double IFORM,double QQ,double SA,double SB);
  TComplex FA1A1P(double XMSQ);
  double WGA1(double QQ);
  double WGA1C(double S);
  double WGA1N(double S);
  //


  vector<double> StructureFunction();
  vector<double> StructureFunctionKuhn();

//---------------------------------------  hadronic structure functions (Kuhn Model) ---------------------------


 private:
 
  TLorentzVector KFitTau_;
  TLorentzVector RecoMuon_;
  TLorentzVector Tau1_;
  TLorentzVector Tau2_;

  TLorentzVector TauMu1_;
  TLorentzVector TauMu2_;
  double _Q;
  double _s1;
  double _s2;
  double _s3;
  bool Flag_;





  TLorentzVector Z_;
  TLorentzVector TauA1_;
  TLorentzVector TauMu_;


  TLorentzVector A1ZFrame_;
  TLorentzVector OSPionZFrame_;
  TLorentzVector SSPion1ZFrame_;
  TLorentzVector SSPion2ZFrame_;



  TLorentzVector _s12;
  TLorentzVector _s13;
  TLorentzVector _s23;

  bool isValid_;

  double mpi;
  double mpi0;
  double mtau;
  double coscab;
  double mrho;
  double mrhoprime;
  double ma1;
  double mpiprime;
  double Gamma0rho; 
  double Gamma0rhoprime; 
  double Gamma0a1;
  double Gamma0piprime;
  double fpi;
  double fpiprime;
  double gpiprimerhopi;
  double grhopipi;
  double beta;
  double COEF1;
  double COEF2;
  double COEF3;
  int SIGN;



};
#endif
