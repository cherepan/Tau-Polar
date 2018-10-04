#include "TauHelper.h"
#include <string.h>


TauHelper::TauHelper(){}
void 
TauHelper::Configure(TLorentzVector OSPion, TLorentzVector SSPion1, TLorentzVector SSPion2,TLorentzVector TauA1, TLorentzVector TauMu ){


  Z_ =TauA1 + TauMu;
  TLorentzVector A1LabFrame = OSPion + SSPion1 + SSPion2;

   // OSPionZFrame_ = Boost(OSPion,Z_);
   // SSPion1ZFrame_= Boost(SSPion1,Z_);
   // SSPion2ZFrame_= Boost(SSPion2,Z_);
   // A1ZFrame_     = Boost(A1LabFrame,Z_);




   // std::cout<<" tau1: ";TauA1.Print();
   // std::cout<<" tau2: ";TauMu.Print();
   // std::cout<<" Z_: ";Z_.Print();

   // std::cout<<" a1 per boost: ";A1LabFrame.Print();
   // std::cout<<" a1 post boost: ";A1ZFrame_.Print();
  OSPionZFrame_ = OSPion;
  SSPion1ZFrame_= SSPion1;
  SSPion2ZFrame_= SSPion2;
  A1ZFrame_     = A1LabFrame;
  TauA1_ = TauA1;
  _Q = A1ZFrame_.M();


  _s12 = SSPion1ZFrame_  + SSPion2ZFrame_;
  _s13 = SSPion1ZFrame_  + OSPionZFrame_;
  _s23 = SSPion2ZFrame_  + OSPionZFrame_;
  _s1  =  _s23.M2(); 
  _s2  =  _s13.M2();
  _s3  =  _s12.M2();


   mpi  = 0.13957018; // GeV 
   mpi0 = 0.1349766;   // GeV
   mtau = 1.776; // GeV
   coscab = 0.975; 
   mrho = 0.773; // GeV
   mrhoprime = 1.370; // GeV
   ma1 = 1.251; // GeV
   mpiprime = 1.300; // GeV
   Gamma0rho  =0.145; // GeV
   Gamma0rhoprime = 0.510; // GeV
   Gamma0a1 = 0.599; // GeV
   Gamma0piprime = 0.3; // GeV
   fpi= 0.093; // GeV
   fpiprime = 0.08; //GeV
   gpiprimerhopi = 5.8; //GeV
   grhopipi = 6.08;  //GeV
   beta = -0.145;
   COEF1 =2.0*sqrt(2.)/3.0;
   COEF2 =-2.0*sqrt(2.)/3.0;
   COEF3 = 2.0*sqrt(2.)/3.0; //C AJW 2/98: Add in the D-wave and I=0 3pi substructure:
   SIGN = -1;



  // std::cout<<"A1 Lab Frame  "<< A1LabFrame.Px() << "  " <<A1LabFrame.Py() << " "<< A1LabFrame.Pz() << "  " <<A1LabFrame.E()<<std::endl;
  // std::cout<<"A1 Z frame    "<< A1ZFrame_.Px() << "  "<< A1ZFrame_.Py() << " "<< A1ZFrame_.Pz() << "  " <<A1ZFrame_.E()<<std::endl;
  // std::cout<<"Z   "<< Z_.Px() << "  "<< Z_.Py() << " "<< Z_.Pz() << "  " <<Z_.E()<<std::endl;

}




void 
TauHelper::SetParametersReco(TLorentzVector tau, TLorentzVector mu ){
 Initialize(tau,mu);
}
void 
TauHelper::SetBoost(TLorentzVector TauA1, TLorentzVector TauMu ){
  Tau1_ = TauA1;
  Tau2_ = TauMu;
}



TauHelper::~TauHelper(){
}



void 
TauHelper::Initialize(TLorentzVector t, TLorentzVector mu){
  RecoMuon_=mu;
  KFitTau_=t;
}

double 
TauHelper::costheta(){
  
  TLorentzVector p1 = SSPion1ZFrame_;
  TLorentzVector p2 = SSPion2ZFrame_;
  TLorentzVector p3 = OSPionZFrame_;
  TLorentzVector Z=Z_;

  
  double zmass = Z.M();
  double mtau = 1.777;
  TLorentzVector a1 = p1+p2+p3;


  float QQ = a1.E()*a1.E() - a1.P()*a1.P();

  double x = 2*a1.E()/zmass;
  double ctheta = (2*x*mtau*mtau - mtau*mtau - QQ)/((mtau*mtau - QQ)*sqrt(1 - 4*mtau*mtau/zmass/zmass));
  // std::cout<<"p1 "<< p1.Px() << "  " <<p1.Py() << " "<< p1.Pz() << "  " <<p1.M()<<std::endl;
  // std::cout<<"p2 "<< p2.Px() << "  "<< p2.Py() << " "<< p2.Pz() << "  " <<p2.M()<<std::endl;
  // std::cout<<"p3 "<< p3.Px() << "  "<< p3.Py() << " "<< p3.Pz() << "  " <<p3.M()<<std::endl;
  // std::cout<<"a1 "<< a1.Px() << "  "<< a1.Py() << " "<< a1.Pz() << "  " <<a1.M()<<"  " <<a1.M()*a1.M() <<std::endl;
  // std::cout<<"QQ "<< QQ << " x "<< x <<" zmass  " <<zmass <<std::endl;

  //  std::cout<<"2*x*mtau*mtau/mtau*mtau - QQ "<< 2*x*mtau*mtau/(mtau*mtau - QQ)<< "  "<< x <<std::endl;
  return ctheta;
}
double 
TauHelper::costheta2(){
  
  TLorentzVector p1 = SSPion1ZFrame_;
  TLorentzVector p2 = SSPion2ZFrame_;
  TLorentzVector p3 = OSPionZFrame_;
  TLorentzVector Z=Z_;

  // double x = LFa1LV.E()/LFtauLV.E();
  // double s = 4*LFtauLV.E()*LFtauLV.E();

  double zmass = Z.M();
  double mt = 1.777;
  TLorentzVector a1 = p1+p2+p3;
  double ma  = a1.M();
  double diffmass = mt*mt - ma*ma;
  float QQ = a1.E()*a1.E() - a1.P()*a1.P();

  double x = a1.E()/TauA1_.E();
  double ctheta = 2*x*mt*mt/diffmass   - (mt*mt + ma*ma)/diffmass;
  return ctheta;

}


double 
TauHelper::costheta1(){
  
  TLorentzVector p1 = SSPion1ZFrame_;
  TLorentzVector p2 = SSPion2ZFrame_;
  TLorentzVector p3 = OSPionZFrame_;
  TLorentzVector Z=Z_;

  double zmass = Z.M();
  double mt = 1.777;
  TLorentzVector a1 = p1+p2+p3;
  double ma  = a1.M();
  double diffmass = mt*mt - ma*ma;
  float QQ = a1.E()*a1.E() - a1.P()*a1.P();

  double x = 2*a1.E()/zmass;
  double ctheta = 4*mt*mt*a1.E()/zmass/diffmass   - (mt*mt + ma*ma)/diffmass;
  return ctheta;

}


float 
TauHelper::CosBeta(){
  TLorentzVector p1 = SSPion1ZFrame_;
  TLorentzVector p2 = SSPion2ZFrame_;
  TLorentzVector p3 = OSPionZFrame_;

  float mpi  = 0.139;
//   float E = p1.E() +  p2.E() +  p3.E(); 
//   float P = p1.P() +  p2.P() +  p3.P(); 
//   float QQ = E*E - P*P;


//   std::cout<<"  cosbeta --------------------- "<<std::endl;

  TLorentzVector a1 = p1+p2+p3;
  float QQ = a1.E()*a1.E() - a1.P()*a1.P();

  float B1 = (pow(p1.E()*a1.E()   - Scalar(p1,a1),2 ) - QQ*mpi*mpi)/QQ;
  float B2 = (pow(p2.E()*a1.E()   - Scalar(p2,a1),2 ) - QQ*mpi*mpi)/QQ;
  float B3 = (pow(p3.E()*a1.E()   - Scalar(p3,a1),2 ) - QQ*mpi*mpi)/QQ;

  float T = 0.5*sqrt(-lambda(B1,B2,B3));

  TLorentzVector p1Timesp2(p1.Py()*p2.Pz() - p1.Pz()*p2.Py(),p1.Pz()*p2.Px() - p1.Px()*p2.Pz(),p1.Px()*p2.Py() - p1.Py()*p2.Px(),1);
    
    
  float cbeta = Scalar(p3,p1Timesp2)/a1.P()/T;

//   std::cout<<"  T  "<< T<<std::endl;
//   std::cout<<"  B1  "<< B1<<std::endl;
//   std::cout<<"  B2  "<< B2<<std::endl;
//   std::cout<<"  B3  "<< B3<<std::endl;
//   std::cout<<"  QQ  "<< QQ<<std::endl;
//   std::cout<<"  Scalar(p3,p1Timesp2)  "<< Scalar(p3,p1Timesp2)<<std::endl;
//   std::cout<<"  cbeta  "<< cbeta  <<std::endl;


  return cbeta;

}


double 
TauHelper::CosBeta1(){
  TLorentzVector p1 = SSPion1ZFrame_;
  TLorentzVector p2 = SSPion2ZFrame_;
  TLorentzVector p3 = OSPionZFrame_;


  //std::cout<<"  cosbeta1  ================================= "<<std::endl;
  double mpi  = 0.139;
//   double E = p1.E() +  p2.E() +  p3.E(); 
//   double P = p1.P() +  p2.P() +  p3.P(); 
//   double QQ = E*E - P*P;

  TLorentzVector a1 = p1+p2+p3;
  //  float P = 
  TLorentzVector s12 = p1+p2;
  TLorentzVector s13 = p1+p3;
  TLorentzVector s23 = p2+p3;

  double QQ = a1.E()*a1.E() - a1.P()*a1.P();


  TLorentzVector p1Timesp2(p1.Py()*p2.Pz() - p1.Pz()*p2.Py(),p1.Pz()*p2.Px() - p1.Px()*p2.Pz(),p1.Px()*p2.Py() - p1.Py()*p2.Px(),1);
  float mm=a1.M()*a1.M();
  float mm12=s12.M()*s12.M();
  float mm13=s13.M()*s13.M();
  float mm23=s23.M()*s23.M();
  float mmpi=mpi*mpi;

  float l1  = lambda( mm, mm12 , mmpi);
  float l2  = lambda( mm, mm13 , mmpi);
  float l3  = lambda( mm, mm23 , mmpi);

 
  double cbeta = /*8*a1.M()*a1.M()**/Scalar(p3,p1Timesp2)*a1.P()/sqrt(-lambda(l1,l2,l3));

//   std::cout<<"  QQ  "<< QQ<<std::endl;
//   std::cout<<"  mm12  "<< mm12<<std::endl;
//   std::cout<<"  mm13  "<< mm13<<std::endl;
//   std::cout<<"  mm23  "<< mm23<<std::endl;
//   std::cout<<"  mm  "<< mm<<std::endl;

//   std::cout<<"  lambda1 = "<<l1<<std::endl;
//   std::cout<<"  lambda2 = "<<l2<<std::endl;
//   std::cout<<"  lambda3 = "<<l3<<std::endl;
//   std::cout<<"  lambda  = "<<-lambda(l1,l2,l3)<<std::endl;

  
//   std::cout<<"  Scalar(p3,p1Timesp2)*a1.P()  "<< Scalar(p3,p1Timesp2)<<std::endl;
//   std::cout<<"  a1.P()  "<< a1.P() <<std::endl;

//   std::cout<<"  a1.P()  "<< a1.P() << "   *  "<</*8*a1.M()*a1.M()**/Scalar(p3,p1Timesp2)/sqrt(-lambda(l1,l2,l3)) <<std::endl;



//   std::cout<<"  cbeta  "<< cbeta  <<std::endl;



  return cbeta;

}




double 
TauHelper::lambda(double x, double y, double z){
  
  return x*x +y*y +z*z - 2*x*y - 2*x*z - 2*z*y;

}




TLorentzVector 
TauHelper::Boost(TLorentzVector pB, TLorentzVector frame){
  
  


  float zmass = frame.M();
  float g  = sqrt(zmass*zmass + pow(frame.Pz(),2  ) +pow(frame.Px(),2  ) +pow(frame.Py(),2  )  )/zmass;
   // float g  = sqrt(zmass*zmass + pow(frame.Pz(),2  ) /*+pow(Frame1_.Px() + Frame2_.Px(),2  ) +pow(Frame1_.Py() + Frame2_.Py(),2  ) */ )/zmass;
  float b  = sqrt(1 - 1/g/g);


  int signx;
  int signy;
  int signz;


  float bx = frame.Px()/frame.E();
  float by = frame.Py()/frame.E();
  float bz = frame.Pz()/frame.E();


  float E = g*pB.E() - g*bx*pB.Px() - g*by*pB.Py() - g*bz*pB.Pz();
  float Px= -g*bx*pB.E() + (1+ (g-1)*bx*bx/b/b)*pB.Px() + ((g-1)*bx*by/b/b)*pB.Py() + ((g-1)*bx*bz/b/b)*pB.Pz();
  float Py= -g*by*pB.E() + ((g-1)*by*bx/b/b)*pB.Px() + (1 + (g-1)*by*by/b/b)*pB.Py() + ((g-1)*by*bz/b/b)*pB.Pz();
  float Pz= -g*bz*pB.E() + ((g-1)*bz*bx/b/b)*pB.Px() + ((g-1)*bz*by/b/b)*pB.Py()    + (1 + (g-1)*bz*bz/b/b)*pB.Pz();

//   float Pz = -sign*b*g*particleToBoost.E() + g*particleToBoost.Pz();
//   float E  = g*particleToBoost.E() - sign*b*g*particleToBoost.Pz();
//   //  std::cout<< "gamma:  "<< g << " beta:  "<<b<<" zmass " <<tau.M() <<std::endl;
  return TLorentzVector(Px,Py,Pz,E);
}


 


float 
TauHelper::Scalar(TLorentzVector p1, TLorentzVector p2){
  
  return p1.Px()*p2.Px() +  p1.Py()*p2.Py() +  p1.Pz()*p2.Pz();
}

std::vector<float> 
TauHelper::Sin2Cos2Gamma(TLorentzVector p1,TLorentzVector p2, TLorentzVector p3){

  std::vector<float> sin2cos2;
  float mpi  = 0.139;
  TLorentzVector a1 = p1+p2+p3;
  float QQ = a1.E()*a1.E() - a1.P()*a1.P();

  float B1 = (pow(p1.E()*a1.E()   - Scalar(p1,a1),2 ) - QQ*mpi*mpi)/QQ;
  float B2 = (pow(p2.E()*a1.E()   - Scalar(p2,a1),2 ) - QQ*mpi*mpi)/QQ;
  float B3 = (pow(p3.E()*a1.E()   - Scalar(p3,a1),2 ) - QQ*mpi*mpi)/QQ;

  float T = 0.5*sqrt(-lambda(B1,B2,B3));

  float A1=(a1.E()*Scalar(a1,p1) - p1.E()*a1.P()*a1.P())/QQ;
  float A2=(a1.E()*Scalar(a1,p2) - p2.E()*a1.P()*a1.P())/QQ;
  float A3=(a1.E()*Scalar(a1,p3) - p3.E()*a1.P()*a1.P())/QQ;


  float cosgamma = A3/a1.P()/sqrt(B3)/sqrt(1 - CosBeta()*CosBeta());
  float singamma = -cosgamma*(B3*A1/A3 - 0.5*(B2 - B1 - B3))/T;

  sin2cos2.push_back(2*singamma*cosgamma);
  sin2cos2.push_back(2*cosgamma*cosgamma - 1);
  return sin2cos2;

}

float 
TauHelper::CosPsi1(){

  TLorentzVector p1 = SSPion1ZFrame_;
  TLorentzVector p2 = SSPion2ZFrame_;
  TLorentzVector p3 = OSPionZFrame_;
  TLorentzVector Z=Z_;

  //  float mtau =1.776;
  TLorentzVector a1 = p1+p2+p3;
  float QQ = a1.E()*a1.E() - a1.P()*a1.P();
  //  float cos = (costheta()*(mtau*mtau  + QQ)   + (mtau*mtau  - QQ))/(costheta()*(mtau*mtau  - QQ)   + (mtau*mtau  + QQ));


  double s = 4*TauA1_.E()*TauA1_.E();
  double x = a1.E()/TauA1_.E();
  if(x*x  - 4*QQ/s <= 0 ){std::cout<<"Warning! In a1Helper::cospsi root square <=0! return 0"<<std::endl; return 0;}
  return    ( x*(mtau*mtau + QQ)  - 2*QQ  )   /   ( mtau*mtau  - QQ   ) / sqrt(x*x  - 4*QQ/s); 
}

 
float 
TauHelper::CosPsi(){

  TLorentzVector p1 = SSPion1ZFrame_;
  TLorentzVector p2 = SSPion2ZFrame_;
  TLorentzVector p3 = OSPionZFrame_;
  TLorentzVector Z=Z_;

  float mtau =1.777;
  TLorentzVector a1 = p1+p2+p3;
  float QQ = a1.E()*a1.E() - a1.P()*a1.P();
  float cos = (costheta()*(mtau*mtau  + QQ)   + (mtau*mtau  - QQ))/(costheta()*(mtau*mtau  - QQ)   + (mtau*mtau  + QQ));
  return cos;





}

 



// //---------------------------------------  hadronic current ---------------------------
// float TauHelper::WA(TLorentzVector s1,TLorentzVector s2,TLorentzVector s3,float QQ){

//   float SS1 = (s2+s3).M()*(s2+s3).M();
//   float SS2 = (s1+s3).M()*(s1+s3).M();
//   float SS3 = (s1+s2).M()*(s1+s2).M();

//  return  VV1(SS1,SS2,SS3,QQ)*F(SS1,SS2,QQ).Rho2() + VV2(SS1,SS2,SS3,QQ)*F(SS2,SS1,QQ).Rho2()  + 2*V1V2(SS1,SS2,SS3,QQ)*( F(SS2,SS1,QQ)*Conjugate(F(SS1,SS2,QQ)) ).Re();
// }

// float TauHelper::WC(TLorentzVector s1,TLorentzVector s2,TLorentzVector s3,float QQ){

//   float SS1 = (s2+s3).M()*(s2+s3).M();
//   float SS2 = (s1+s3).M()*(s1+s3).M();
//   float SS3 = (s1+s2).M()*(s1+s2).M();

//   return  (VV1(SS1,SS2,SS3,QQ) - 2*h(SS1,SS2,SS3,QQ))*F(SS1,SS2,QQ).Rho2() + (VV2(SS1,SS2,SS3,QQ) - 2*h(SS1,SS2,SS3,QQ))*F(SS2,SS1,QQ).Rho2() +   (2*V1V2(SS1,SS2,SS3,QQ) + 4*h(SS1,SS2,SS3,QQ))*( F(SS2,SS1,QQ)*Conjugate(F(SS1,SS2,QQ)) ).Re();
// } 

// float
// TauHelper::WD(TLorentzVector s1,TLorentzVector s2,TLorentzVector s3,float QQ){
//   float SS1 = (s2+s3).M()*(s2+s3).M();
//   float SS2 = (s1+s3).M()*(s1+s3).M();
//   float SS3 = (s1+s2).M()*(s1+s2).M();
//   float mpi = 0.139;
//   float undersqrt1 = VV1(SS1,SS2,SS3,QQ)  -h(SS1,SS2,SS3,QQ);
//   float undersqrt2 = VV2(SS1,SS2,SS3,QQ)  -h(SS1,SS2,SS3,QQ);

//   if(undersqrt1 < 0) undersqrt1 =0;
//   if(undersqrt2 < 0) undersqrt2 =0;

//   //std::cout<<" debug Wd "<<h(SS1,SS2,SS3,QQ)<<" VV1(SS1,SS2,SS3,QQ)   "<<VV1(SS1,SS2,SS3,QQ) <<"  "<<VV2(SS1,SS2,SS3,QQ)<<"   "<<-h0(SS1,SS2,SS3,QQ) <<"  "<< ( F(SS2,SS1,QQ)*Conjugate(F(SS1,SS2,QQ)) ).Re()<<"   " <<VV1(SS1,SS2,SS3,QQ)  -h(SS1,SS2,SS3,QQ) << "    " << VV2(SS1,SS2,SS3,QQ)  -h(SS1,SS2,SS3,QQ)<<std::endl;
//   //  return  -sqrt(h(SS1,SS2,SS3,QQ))*(2*sqrt(VV1(SS1,SS2,SS3,QQ)  -h(SS1,SS2,SS3,QQ)) *F(SS1,SS2,QQ).Rho2() - 2*sqrt(VV2(SS1,SS2,SS3,QQ)  -h(SS1,SS2,SS3,QQ))*F(SS2,SS1,QQ).Rho2()  
//   //				    + (QQ- mpi*mpi + SS3)*(SS1 - SS2 )*( F(SS2,SS1,QQ)*Conjugate(F(SS1,SS2,QQ)) ).Re()/QQ/sqrt(-h0(SS1,SS2,SS3,QQ)));

//   return  -sqrt(h(SS1,SS2,SS3,QQ))*(2*sqrt(undersqrt1) *F(SS1,SS2,QQ).Rho2() - 2*sqrt(undersqrt2)*F(SS2,SS1,QQ).Rho2()  
// 				    + (QQ- mpi*mpi + SS3)*(SS1 - SS2 )*( F(SS2,SS1,QQ)*Conjugate(F(SS1,SS2,QQ)) ).Re()/QQ/sqrt(-h0(SS1,SS2,SS3,QQ)));


// }

// float
// TauHelper::WE(TLorentzVector s1,TLorentzVector s2,TLorentzVector s3,float QQ){
//   float SS1 = (s2+s3).M()*(s2+s3).M();
//   float SS2 = (s1+s3).M()*(s1+s3).M();
//   float SS3 = (s1+s2).M()*(s1+s2).M();
//   return  -3*sqrt(-h(SS1,SS2,SS3,QQ)*h0(SS1,SS2,SS3,QQ))*( F(SS2,SS1,QQ)*Conjugate(F(SS1,SS2,QQ)) ).Im();

// }


TComplex 
TauHelper::F(float s1, float s2,float QQ){


  float beta= -0.145;
  float fpi= 0.093;


  TComplex factor(0,-2*sqrt(2)/3/fpi/(1+beta));
  TComplex BreighWignerRhoPrime(beta*BWrhoPrime(s2).Re(), beta*BWrhoPrime(s2).Im());
  TComplex out = factor*BWa1(QQ)*(BWrho(s2) + BreighWignerRhoPrime);
  return out;
}


float
TauHelper::VV1(float SS1 ,float SS2, float SS3, float QQ){
  float mpi = 0.139;
  return  SS2 - 4*mpi*mpi + pow(SS3 - SS1,2)/4/QQ;
}

float
TauHelper::VV2(float SS1 ,float SS2, float SS3, float QQ){
  float mpi = 0.139;
  return  SS1 - 4*mpi*mpi + pow(SS3 - SS2,2)/4/QQ;
}

float
TauHelper::V1V2(float SS1 ,float SS2, float SS3, float QQ){
  float mpi = 0.139;
  return  (QQ/2 - SS3 - mpi*mpi/2) + (SS3 - SS1)*(SS3 - SS2)/4/QQ;
}


float
TauHelper::h0(float SS1 ,float SS2, float SS3, float QQ){
  float mpi = 0.139;
  return 4*mpi*mpi - pow(2*mpi*mpi - SS1 - SS2,2)/QQ;
}

float
TauHelper::h(float SS1 ,float SS2, float SS3, float QQ){
  float mpi = 0.139;
  return -(SS1*SS2*SS3 - mpi*mpi*pow(QQ - mpi*mpi,2))/h0(SS1,SS2,SS3,QQ)/QQ;
}



TComplex 
TauHelper::BWa1(float QQ){
  float m =  1.251;
  TComplex re,im;
  re = (m*m*pow(m*m - QQ,2))/(pow(m*m - QQ,2) - m*m*GammaA1(QQ)*GammaA1(QQ));
  im = m*m*m*GammaA1(QQ)/(pow(m*m - QQ,2) - m*m*GammaA1(QQ)*GammaA1(QQ));
  TComplex out(re,im);
  return out;
}



TComplex 
TauHelper::BWrho(float QQ){
  float m =  0.773;
  TComplex re,im;
  re = (m*m*pow(m*m - QQ,2))/(pow(m*m - QQ,2) - m*m*GammaRho(QQ)*GammaRho(QQ));
  im = m*m*m*GammaRho(QQ)/(pow(m*m - QQ,2) - m*m*GammaRho(QQ)*GammaRho(QQ));
  TComplex out(re,im);
  return out;
}




TComplex 
TauHelper::BWrhoPrime(float QQ){
  float m =  1.251;
  TComplex re,im;
  re = (m*m*pow(m*m - QQ,2))/(pow(m*m - QQ,2) - m*m*GammaRhoPrime(QQ)*GammaRhoPrime(QQ));
  im = m*m*m*GammaRhoPrime(QQ)/(pow(m*m - QQ,2) - m*m*GammaRhoPrime(QQ)*GammaRhoPrime(QQ));
  TComplex out(re,im);
  return out;
}


float
TauHelper::GammaA1(float QQ){
  float ma1 = 1.251;
  float ga1 = 0.599;
  float out = ga1*gForGammaA1(QQ)/gForGammaA1(ma1*ma1);
  return out;
}

float
TauHelper::gForGammaA1(float QQ){
  float mpi  = 0.139;
  float mrho = 0.773;
  float out;
  if(QQ > pow((mrho + mpi),2)){ out = QQ*(1.623 + 10.38/QQ - 9.34/QQ/QQ + 0.65/QQ/QQ/QQ);}
  else out = 4.1*pow((QQ - 9*mpi*mpi),3)*(1- 3.3*(QQ - 9*mpi*mpi) + 5.8*pow(QQ - 9*mpi*mpi,2));
  return out;
}


float
TauHelper::GammaRho(float QQ){
  float mpi  = 0.139;
  float mrho = 0.773;
  float grho = 0.599;
  float out  =grho*mrho*pow( sqrt(QQ - 4*mpi*mpi)/sqrt(mrho*mrho - mpi*mpi)   ,3)/sqrt(QQ);
  return out;
}


float
TauHelper::GammaRhoPrime(float QQ){

  float mrhoPrime = 1.370;
  float grhoPrime = 0.510;
  float out  =grhoPrime*QQ/mrhoPrime/mrhoPrime;
  return out;
}


double
TauHelper::GetOmegaA1(){
        isValid_ = false;
        double omega(-999.);
	TLorentzVector pi1 = SSPion1ZFrame_;
	TLorentzVector pi2 = SSPion2ZFrame_;
	TLorentzVector pi3 = OSPionZFrame_;
	TLorentzVector a = A1ZFrame_;
	float mtau = 1.777;
	float cospsi  = CosPsi1();
	float sinpsi  = sqrt(1 - cospsi*cospsi);
	float sin2psi = 2*sinpsi*cospsi;
	
	float sin2gamma = Sin2Cos2Gamma(pi1,pi2,pi3).at(0);
	float cos2gamma = Sin2Cos2Gamma(pi1,pi2,pi3).at(1);
	
	
	float cosbeta=CosBeta();
	float sinbeta = sqrt(1  - cosbeta*cosbeta);
	
	  
	float costheta=costheta2();
	float sintheta = sqrt(1 - costheta*costheta);

	double RR  = mtau*mtau/a.M()/a.M();
	float U = 0.5*(3*cospsi*cospsi - 1)*(1 - RR);
	float V = 0.5*(3*cospsi*cospsi - 1)*(1 + RR)*costheta + 0.5*3*sin2psi*sintheta*sqrt(RR);

	// float Wa =WA();
	// float Wc =WC();
	// float Wd =WD();
	// float We =WE();


	// std::cout<<" Wa "<< WA() << " cleo wa  "<< StructureFunctionKuhn().at(0) <<std::endl;
	// std::cout<<" Wc "<< WC() << " cleo wc  "<< StructureFunctionKuhn().at(1) <<std::endl;
	// std::cout<<" Wd "<< WD() << " cleo wd  "<< StructureFunctionKuhn().at(2) <<std::endl;
	// std::cout<<" We "<< WE() << " cleo we  "<< StructureFunctionKuhn().at(3) <<std::endl;


	   float Wa =StructureFunctionKuhn().at(0);
	   float Wc =StructureFunctionKuhn().at(1);
	   float Wd =StructureFunctionKuhn().at(2);
	   float We =StructureFunctionKuhn().at(3);


	//   float Wa =StructureFunction().at(0);
	//   float Wc =StructureFunction().at(1);
	//   float Wd =StructureFunction().at(2);
	//   float We =StructureFunction().at(3);

	// std::cout<<"  "<< StructureFunction().at(0)*StructureFunction().at(0) << " =  "<< StructureFunction().at(1)*StructureFunction().at(1)+StructureFunction().at(2)*StructureFunction().at(2)+St	  ructureFunction().at(3)*StructureFunction().at(3) <<std::endl;

	float fa1 = (2  + RR + 0.5*(3*cosbeta*cosbeta - 1)*U)*Wa/3 - 0.5*sinbeta*sinbeta*cos2gamma*U*Wc + 0.5*sinbeta*sinbeta*sin2gamma*U*Wd + cospsi*cosbeta*We;
	float ga1 = (costheta*(RR -2) - 0.5*(3*cosbeta*cosbeta - 1)*V)*Wa/3 + 0.5*sinbeta*sinbeta*cos2gamma*V*Wc - 0.5*sinbeta*sinbeta*sin2gamma*V*Wd -cosbeta*(costheta*cospsi + sintheta*sinpsi*sqrt(RR))*We;
	omega = ga1/fa1;
	// if((float)ga1/fa1 > 1.05) omega =1.0;
	// if((float)ga1/fa1 < -1.05) omega =-1.0;

	// if(omega > 1 or omega < -1){
	//    std::cout<<" omega  "<< omega<<std::endl;
	//    std::cout<<" a1   "<< a.M()<<std::endl;

	//	std::cout<<"  in Omega "<< Wa*Wa << "  =  "<< Wc*Wc + Wd*Wd + We*We <<std::endl;
	//    std::cout<< " U  "<< U <<" V  "<<V <<" Wa "<<Wa <<" Wc  "<<Wc  <<" Wd  "<< Wd <<" We  "<< We  <<" f "<< fa1 <<" g "<< ga1 <<std::endl;
	//    std::cout<< " cospsi  "<< cospsi <<" costheta  "<<costheta <<" sin2psi "<<sin2psi <<" sintheta  "<<sintheta  <<" RR  "<< RR <<std::endl;
	// }

	//	PVC().Print();
	//	std::cout<<"---------- omega 1 "<< omega<<std::endl;
	if(omega > 1) omega == 1;  // in some events omega might exceed 1, this is when a1.M() - a1.M0( )~ Gamma_a1
	if(omega <-1) omega ==-1;
	if(omega > 0 or omega < 0) isValid_ = true;
	//	std::cout<<"--------   omega 2 "<< omega<<std::endl;
	return omega;
}

double
TauHelper::getA1omegaIntegratedOverTheta(){



        isValid_ = false;
        double omega(-999.);
	TLorentzVector pi1 = SSPion1ZFrame_;
	TLorentzVector pi2 = SSPion2ZFrame_;
	TLorentzVector pi3 = OSPionZFrame_;
	TLorentzVector a = A1ZFrame_;
	float mtau = 1.777;
	float cospsi  = CosPsi();
	float sinpsi  = sqrt(1 - cospsi*cospsi);
	float sin2psi = 2*sinpsi*cospsi;
	
	float sin2gamma = Sin2Cos2Gamma(pi1,pi2,pi3).at(0);
	float cos2gamma = Sin2Cos2Gamma(pi1,pi2,pi3).at(1);
	
	
	float cosbeta=CosBeta();
	float sinbeta = sqrt(1  - cosbeta*cosbeta);
	
	  
	float costheta=costheta2();
	float sintheta = sqrt(1 - costheta*costheta);

	double RR  = mtau*mtau/a.M()/a.M();
	float U = 0.5*(3*cospsi*cospsi - 1)*(1 - RR);
	float V = 0.5*(3*cospsi*cospsi - 1)*(1 + RR)*costheta + 0.5*3*sin2psi*sintheta*sqrt(RR);

	double U_integrated = (1-RR)*0.5*(3*Integral_squared_cospsi()   - 2 );
	double V_integrated = 3*0.5*(1+RR)*Integral_squared_cospsi_costheta() + 3*0.5*sqrt(RR)*Integral_sin2psi_sintheta();

	 // float Wa =WA();
	 // float Wc =WC();
	 // float Wd =WD();
	 // float We =WE();

	   float Wa =StructureFunctionKuhn().at(0);
	   float Wc =StructureFunctionKuhn().at(1);
	   float Wd =StructureFunctionKuhn().at(2);
	   float We =StructureFunctionKuhn().at(3);

	 // float Wa =StructureFunction().at(0);
	 // float Wc =StructureFunction().at(1);
	 // float Wd =StructureFunction().at(2);
	 // float We =StructureFunction().at(3);



	// float Wa =WA(pi1,pi2,pi3,a.M()*a.M());
	// float Wc =WC(pi1,pi2,pi3,a.M()*a.M());
	// float Wd =WD(pi1,pi2,pi3,a.M()*a.M());
	// float We =WE(pi1,pi2,pi3,a.M()*a.M());

	float fa1 = ( 2*(2  + RR) + 0.5*(3*cosbeta*cosbeta - 1)*U_integrated)*Wa/3 - 0.5*sinbeta*sinbeta*cos2gamma*U_integrated*Wc
	  + 0.5*sinbeta*sinbeta*sin2gamma*U_integrated*Wd + Integral_cospsi()*cosbeta*We;

	float ga1 = (0*costheta*(RR -2) - 0.5*(3*cosbeta*cosbeta - 1)*V_integrated)*Wa/3 +
	  0.5*sinbeta*sinbeta*cos2gamma*V_integrated*Wc - 
	  0.5*sinbeta*sinbeta*sin2gamma*V_integrated*Wd -
	  cosbeta*(Integral_cospsi_costheta() + Integral_sinpsi_sintheta()*sqrt(RR))*We;





//   double QQ=_Q*_Q;
//   double RR  = mtau*mtau/QQ;

//   double U_integrated = (1-RR)*0.5*(3*Integral_squared_cospsi()   - 2 );
//   double V_integrated = 3*0.5*(1+RR)*Integral_squared_cospsi_costheta() + 3*0.5*sqrt(RR)*Integral_sin2psi_sintheta();
//   double Beta = 0.5*(3*cosbetaLF()*cosbetaLF()- 1);
//   float U = 0.5*(3*cospsiLF()*cospsiLF() - 1)*(1 - RR);
//   float V = 0.5*(3*cospsiLF()*cospsiLF() - 1)*(1 + RR)*costhetaLF() + 0.5*3*2*cospsiLF()* sinpsiLF()*sinthetaLF()*sqrt(RR);
  
//   double f_beta_gamma = (  (2  + RR) + Beta*U_integrated)*WA()/3
//     - 0.5*sinbetaLF()*sinbetaLF()*cos2gammaLF()*U_integrated*WC() + 0.5*sinbetaLF()*sinbetaLF()*sin2gammaLF()*U_integrated*WD() + Integral_cospsi()*cosbetaLF()*WE();



//   double g_beta_gamma = ( (RR -2) - Beta*V_integrated)*WA()/3 + 0.5*sinbetaLF()*sinbetaLF()*cos2gammaLF()*V_integrated*WC() - 
//     0.5*sinbetaLF()*sinbetaLF()*sin2gammaLF()*V_integrated*WD() -cosbetaLF()*( Integral_cospsi_costheta() + Integral_sinpsi_sintheta()*sqrt(RR))*WE();


//   double omega = g_beta_gamma/f_beta_gamma;

//   //  if(omega > 1)
// {
//    // std::cout<<"*********************  "<<std::endl;
//    // std::cout<<"  f  "<< f_beta_gamma  <<"   g   "<<  g_beta_gamma  <<"   omega   "<< g_beta_gamma/f_beta_gamma <<std::endl;

//    // std::cout<< "INtegrals summaary "<< std::endl;
//    // std::cout<<"  RR   "<< RR << std::endl;
//    // std::cout<<" U_integrated  "<< U_integrated<< "   U   " <<  U <<std::endl;
//    // std::cout<<" V_integrated  "<<V_integrated << "  V   "<<  V <<std::endl;
//    // std::cout<<" Integral_squared_cospsi  "<< Integral_squared_cospsi()<< std::endl;
//    // std::cout<<" Integral_squared_cospsi_costheta  "<< Integral_squared_cospsi_costheta()<< std::endl;
//    // std::cout<<" Integral_sin2psi_sintheta  "<< Integral_sin2psi_sintheta()<< std::endl;
//    // std::cout<<" Integral_cospsi  "<< Integral_cospsi()<< std::endl;
//    // std::cout<<" Integral_sinpsi_sintheta  "<< Integral_sinpsi_sintheta()<< std::endl;
//    // std::cout<<" Integral_cospsi_costheta  "<< Integral_cospsi_costheta()<< std::endl;
//   }
	omega = ga1/fa1;



	// if(omega > 1 or omega < -1){
	//    std::cout<<" omega  "<< omega<<std::endl;
	//    std::cout<<" a1   "<< a.M()<<std::endl;

	//    std::cout<<"  "<< Wa*Wa << "  =  "<< Wc*Wc + Wd*Wd + We*We <<std::endl;
	//    std::cout<< " U  "<< U <<" V  "<<V <<" Wa "<<Wa <<" Wc  "<<Wc  <<" Wd  "<< Wd <<" We  "<< We  <<" f "<< fa1 <<" g "<< ga1 <<std::endl;
	//    std::cout<< " cospsi  "<< cospsi <<" costheta  "<<costheta <<" sin2psi "<<sin2psi <<" sintheta  "<<sintheta  <<" RR  "<< RR <<std::endl;
	// }

	if(omega > 0 or omega < 0) 	return omega;
	return -999;
}




double
TauHelper::I_squared_cospsi(double limit){

  double QQ=_Q*_Q;
  double beta_q = (mtau*mtau - QQ)/(mtau*mtau + QQ);
  return  (limit + 2 * (beta_q*beta_q -1)*log(limit) - pow(beta_q*beta_q-1,2)/limit  )/pow(beta_q,3);
}

double
TauHelper::Integral_squared_cospsi(){

  double QQ=_Q*_Q;
  double beta_q = (mtau*mtau - QQ)/(mtau*mtau + QQ);

  return I_squared_cospsi(1+beta_q) -I_squared_cospsi(1-beta_q);
}



double
TauHelper::I_cospsi(double limit){

  double QQ=_Q*_Q;
  double RR  = mtau*mtau/QQ;
  double beta_q = (mtau*mtau - QQ)/(mtau*mtau + QQ);

  return  ( limit + (beta_q*beta_q - 1)*log(limit) )/pow(beta_q,2);
}


double
TauHelper::Integral_cospsi(){

  double QQ=_Q*_Q;
  double beta_q = (mtau*mtau - QQ)/(mtau*mtau + QQ);

  return  I_cospsi(1+beta_q) - I_cospsi(1-beta_q);
}


double
TauHelper::I_squared_cospsi_costheta(double limit){

  double QQ=_Q*_Q;
  double RR  = mtau*mtau/QQ;
  double beta_q = (mtau*mtau - QQ)/(mtau*mtau + QQ);

  return  (limit*limit*0.5 + limit*(2*beta_q  -3) + ( beta_q*beta_q - 4 *beta_q +3)*log(limit) + pow(beta_q-1,2)/limit)/pow(beta_q,4);
}

double
TauHelper::Integral_squared_cospsi_costheta(){

  double QQ=_Q*_Q;
  double beta_q = (mtau*mtau - QQ)/(mtau*mtau + QQ);

  return  I_squared_cospsi_costheta(1+beta_q)  - I_squared_cospsi_costheta(1-beta_q);
}



double
TauHelper::I_cospsi_costheta(double limit){

  double QQ=_Q*_Q;
  double beta_q = (mtau*mtau - QQ)/(mtau*mtau + QQ);

  return  (limit*limit*0.5 + limit*(beta_q*beta_q  -2) + (1-  beta_q*beta_q )*log(limit))/pow(beta_q,3);

}


double
TauHelper::Integral_cospsi_costheta(){

  double QQ=_Q*_Q;
  double beta_q = (mtau*mtau - QQ)/(mtau*mtau + QQ);

  return  I_cospsi_costheta(1+beta_q) - I_cospsi_costheta(1-beta_q);

}


double
TauHelper::I_sinpsi_sintheta(double limit){

  double QQ=_Q*_Q;
  double beta_q = (mtau*mtau - QQ)/(mtau*mtau + QQ);

  return  sqrt(1-beta_q*beta_q)*(limit*2 - limit*limit*0.5 +  (  beta_q*beta_q -1  )*log(limit))/pow(beta_q,3);
}

double
TauHelper::Integral_sinpsi_sintheta(){

  double QQ=_Q*_Q;
  double beta_q = (mtau*mtau - QQ)/(mtau*mtau + QQ);

  return  I_sinpsi_sintheta(1+beta_q) - I_sinpsi_sintheta(1-beta_q);
}


double
TauHelper::I_sin2psi_sintheta(double limit){

  double QQ=_Q*_Q;
  double beta_q = (mtau*mtau - QQ)/(mtau*mtau + QQ);

  return  sqrt(1-beta_q*beta_q)*( -pow(limit,3)/3 + 3*0.5*limit*limit + limit*(beta_q*beta_q - 3) + (1-beta_q*beta_q)*log(limit)    )/pow(beta_q,4);
}

double
TauHelper::Integral_sin2psi_sintheta(){

  double QQ=_Q*_Q;
  double beta_q = (mtau*mtau - QQ)/(mtau*mtau + QQ);

  return  I_sin2psi_sintheta(1+beta_q)  -  I_sin2psi_sintheta(1-beta_q);
}



double
TauHelper::VV1(){ //  this is -V1^{2}
  double QQ = _Q*_Q;
  return  _s2 - 4*mpi*mpi + pow(_s3 - _s1,2)/4/QQ;
}



double
TauHelper::VV2(){ //  this is -V2^{2}
  double QQ = _Q*_Q;
  return  _s1 - 4*mpi*mpi + pow(_s3 - _s2,2)/4/QQ;
}

double
TauHelper::V1V2(){  // this is -V1V2
  double QQ = _Q*_Q;
  return  (QQ/2 - _s3 - 0.5*mpi*mpi) + (_s3 - _s1)*(_s3 - _s2)/4/QQ;
}


double
TauHelper::h0(){ // this is -3sqrt{h0}/2
  double QQ = _Q*_Q;
  return -4*mpi*mpi + pow(2*mpi*mpi - _s1 - _s2,2)/QQ;
}

double
TauHelper::h(){
  double QQ = _Q*_Q;
  return (_s1*_s2*_s3 - mpi*mpi*pow(QQ - mpi*mpi,2))/h0()/QQ;  // this is sqrt{h}
}

double 
TauHelper::WA(){
  //  std::cout<<" WA  F1  " << F1()<< std::endl;

  return  VV1()*F1().Rho2() + VV2()*F2().Rho2()  + 2*V1V2()*( F1()*Conjugate(F2()) ).Re();
 }

 double 
TauHelper::WC(){
   return  -(-VV1() + 2*h() )*F1().Rho2() - (-VV2() + 2*h())*F2().Rho2()   -   (-2*V1V2() - 4*h())*( F1()*Conjugate(F2()) ).Re();
 } 

 double
 TauHelper::WD(){
   double QQ = _Q*_Q;
   double undersqrt1 = VV1()  -h();
   double undersqrt2 = VV2()  -h();
   return  -sqrt(h()) * ( 2 * sqrt(undersqrt1) * F1().Rho2() - 2*sqrt(undersqrt2)*F2().Rho2()  
			  + (QQ - mpi*mpi + _s3)*(_s1 - _s2 )*( F1()*Conjugate(F2()) ).Re()/QQ/sqrt(h0() ) );


 }

 double
 TauHelper::WE(){
  return  3*sqrt(h()*h0())*( F1()*Conjugate(F2()) ).Im();
 }






TComplex 
TauHelper::F1(){
  TComplex scale(0, -2*sqrt(2)/3/fpi);
  TComplex res = scale*BreitWigner(_Q,"a1")*BRho(sqrt(_s2));
  return res;
}

TComplex 
TauHelper::F2(){
  TComplex scale(0, -2*sqrt(2)/3/fpi);
  TComplex res = scale*BreitWigner(_Q,"a1")*BRho(sqrt(_s1));
  return res;
}

TComplex 
TauHelper::F4(){
  TComplex scale(0, -gpiprimerhopi*grhopipi*fpiprime/2/pow(mrho,4)/pow(mpiprime,3));
  TComplex res = scale*BreitWigner(_Q,"piprime")*(_s1*(_s2-_s3)*BRho(sqrt(_s1)) + _s2*(_s1-_s3)*BRho(sqrt(_s2)));
  return res;
}

TComplex 
TauHelper::BRho(double Q){
   return (BreitWigner(Q) + beta*BreitWigner(Q,"rhoprime"))/(1+beta);
}

TComplex 
TauHelper::BreitWigner(double Q, string type){
  double QQ=Q*Q;
  double re,im;
  double m = Mass(type);
  double g  = Widths(Q,type);
  re = (m*m*(m*m - QQ))/(pow(m*m - QQ,2) + m*m*g*g);
  im = Q*g/(pow(m*m - QQ,2) + m*m*g*g);
 

  TComplex out(re,im);
  return out;
}

double
TauHelper::Widths(double Q, string type){
  double QQ = Q*Q;
  double Gamma;
  Gamma = Gamma0rho*mrho*pow( ppi(QQ)  / ppi(mrho*mrho), 3) /sqrt(QQ);
  if(type == "rhoprime"){
    Gamma=Gamma0rhoprime*QQ/mrhoprime/mrhoprime;
 }
  if(type == "a1"){
    Gamma=Gamma0a1*ga1(Q)/ga1(ma1);
    //    Gamma=Gamma0a1*ma1*ga1(Q)/ga1(ma1)/Q;
 }
  if(type == "piprime"){
    Gamma = Gamma0piprime*pow( sqrt(QQ)/mpiprime  ,5)*pow( (1-mrho*mrho/QQ)/(1-mrho*mrho/mpiprime/mpiprime) ,3);
  }
  return Gamma;
}
double TauHelper::ga1(double  Q){
  double QQ = Q*Q;
  return (QQ > pow(mrho + mpi,2)) ?  QQ*(1.623 + 10.38/QQ - 9.32/QQ/QQ   + 0.65/QQ/QQ/QQ)  : 4.1*pow(QQ - 9*mpi*mpi,3)*(  1 - 3.3*(QQ - 9*mpi*mpi)  + 5.8*pow(QQ - 9*mpi*mpi,2)  );
}
double
TauHelper::Mass(string type){
  double m = mrho;
  if(type == "rhoprime") return mrhoprime; 
  if(type == "a1") return ma1;
  if(type == "piprime") return mpiprime;
  return m;
}
double TauHelper::ppi(double QQ){  if(QQ < 4*mpi*mpi) std::cout<<"Warning! Can not compute ppi(Q); root square <0 ; return nan  "; return 0.5*sqrt(QQ - 4*mpi*mpi);}


TComplex 
TauHelper::Conjugate(TComplex a){
  return TComplex(a.Re(), -a.Im());
}






//------- L-wave BreightWigner for rho
TComplex 
TauHelper::BWIGML(double S, double M,  double G, double m1, double m2, int L){
  int IPOW;
  double MP = pow(m1+m2,2);
  double MM = pow(m1-m2,2);
  double MSQ = M*M;
  double W = sqrt(S);
  double WGS =0.0;
  double QS,QM;
  if(W > m1+m2){
    QS = sqrt(std::fabs( (S  - MP)*(S  - MM)))/W;
    QM = sqrt(std::fabs( (MSQ - MP)*(MSQ - MM)))/M;
    IPOW = 2*L +1;
    WGS=G*(MSQ/W)*pow(QS/QM, IPOW);
  }

 TComplex out;
 out = TComplex(MSQ,0)/TComplex(MSQ - S, -WGS) ;
 return out;
}

TComplex
TauHelper::FPIKM(double W, double XM1, double XM2){
  double ROM  = 0.773;
  double ROG  = 0.145;
  double ROM1 = 1.370;
  double ROG1 = 0.510;
  double BETA1=-0.145;
  
  double S=W*W;
  int L =1; // P-wave
  TComplex out = (BWIGML(S,ROM,ROG,XM1,XM2,L) + BETA1*BWIGML(S,ROM1,ROG1,XM1,XM2,L))/(1+BETA1);
  return out;
} 


TComplex
TauHelper::F3PI(double IFORM,double QQ,double SA,double SB){
  double MRO = 0.7743;
  double GRO = 0.1491;
  double MRP = 1.370 ;
  double GRP = 0.386 ;
  double MF2 = 1.275;
  double GF2 = 0.185;
  double MF0 = 1.186;
  double GF0 = 0.350;
  double MSG = 0.860;
  double GSG = 0.880;
  double MPIZ = mpi0;
  double MPIC = mpi;
  double M1;
  double M2;
  double M3;
  int IDK =1;  // --------- it is 3pi
  if(IDK ==1){
    M1=MPIZ;
    M2=MPIZ;
    M3=MPIC;
  }
  if(IDK==2){
    M1=MPIC;
    M2=MPIC;
    M3=MPIC;
  }

 
  double M1SQ = M1*M1;
  double M2SQ = M2*M2;
  double M3SQ = M3*M3;
  
  // parameter varioation for
  // systematics   from, https://arxiv.org/pdf/hep-ex/9902022.pdf

  double db2 = 0.094;   double dph2 = 0.253;
  double db3 = 0.094;   double dph3 = 0.104;
  double db4 = 0.296;   double dph4 = 0.170;
  double db5 = 0.167;   double dph5 = 0.104;
  double db6 = 0.284;   double dph6 = 0.036;
  double db7 = 0.148;   double dph7 = 0.063;

  double scale(0.);
  // if(doSystematic)
  //   {
  //     if(systType=="UP")
  // 	{
  // 	  scale =  1;
  // 	}
      
      
  //     if(systType=="DOWN")
  // 	{
  // 	  scale = -1;
  // 	}
  //   } 
      
  TComplex  BT1 = TComplex(1.,0.);
  TComplex  BT2 = TComplex(0.12  + scale*db2,0.)*TComplex(1, (0.99   +  scale*dph2)*TMath::Pi(), true);//  TComplex(1, 0.99*TMath::Pi(), true);   Real part must be equal to one, stupid polar implemenation in root
  TComplex  BT3 = TComplex(0.37  + scale*db3,0.)*TComplex(1, (-0.15  +  scale*dph3)*TMath::Pi(), true);
  TComplex  BT4 = TComplex(0.87  + scale*db4,0.)*TComplex(1, (0.53   +  scale*dph4)*TMath::Pi(), true);
  TComplex  BT5 = TComplex(0.71  + scale*db5,0.)*TComplex(1, (0.56   +  scale*dph5)*TMath::Pi(), true);
  TComplex  BT6 = TComplex(2.10  + scale*db6,0.)*TComplex(1, (0.23   +  scale*dph6)*TMath::Pi(), true);
  TComplex  BT7 = TComplex(0.77  + scale*db7,0.)*TComplex(1, (-0.54  +  scale*dph7)*TMath::Pi(), true);

  TComplex  F3PIFactor(0.,0.); // initialize to zero
  if(IDK == 2){
    if(IFORM == 1 || IFORM == 2 ){
      double S1 = SA;
      double S2 = SB;
      double S3 = QQ-SA-SB+M1SQ+M2SQ+M3SQ;
      //Lorentz invariants for all the contributions:
      double F134 = -(1./3.)*((S3-M3SQ)-(S1-M1SQ));
      double F15A = -(1./2.)*((S2-M2SQ)-(S3-M3SQ));
      double F15B = -(1./18.)*(QQ-M2SQ+S2)*(2.*M1SQ+2.*M3SQ-S2)/S2;
      double F167 = -(2./3.);
      
      // Breit Wigners for all the contributions:
      
      
      TComplex  FRO1 = BWIGML(S1,MRO,GRO,M2,M3,1);
      TComplex  FRP1 = BWIGML(S1,MRP,GRP,M2,M3,1);
      TComplex  FRO2 = BWIGML(S2,MRO,GRO,M3,M1,1);
      TComplex  FRP2 = BWIGML(S2,MRP,GRP,M3,M1,1);
      TComplex  FF21 = BWIGML(S1,MF2,GF2,M2,M3,2);
      TComplex  FF22 = BWIGML(S2,MF2,GF2,M3,M1,2);
      TComplex  FSG2 = BWIGML(S2,MSG,GSG,M3,M1,0);
      TComplex  FF02 = BWIGML(S2,MF0,GF0,M3,M1,0);
      
      
      F3PIFactor = BT1*FRO1+BT2*FRP1+
	BT3*TComplex(F134,0.)*FRO2+BT4*TComplex(F134,0.)*FRP2
	-BT5*TComplex(F15A,0.)*FF21-BT5*TComplex(F15B,0.)*FF22
	-BT6*TComplex(F167,0.)*FSG2-BT7*TComplex(F167,0.)*FF02;
      
    } else if (IFORM == 3 ){
      
      double S3 = SA;
      double S1 = SB;
      double S2 = QQ-SA-SB+M1SQ+M2SQ+M3SQ;
      
      double F34A = (1./3.)*((S2-M2SQ)-(S3-M3SQ));
      double F34B = (1./3.)*((S3-M3SQ)-(S1-M1SQ));
      double F35A = -(1./18.)*(QQ-M1SQ+S1)*(2.*M2SQ+2.*M3SQ-S1)/S1;
      double F35B =  (1./18.)*(QQ-M2SQ+S2)*(2.*M3SQ+2.*M1SQ-S2)/S2;
      double F36A = -(2./3.);
      double F36B =  (2./3.);
      
      //C Breit Wigners for all the contributions:
      TComplex  FRO1 = BWIGML(S1,MRO,GRO,M2,M3,1);
      TComplex  FRP1 = BWIGML(S1,MRP,GRP,M2,M3,1);
      TComplex  FRO2 = BWIGML(S2,MRO,GRO,M3,M1,1);
      TComplex  FRP2 = BWIGML(S2,MRP,GRP,M3,M1,1);
      TComplex  FF21 = BWIGML(S1,MF2,GF2,M2,M3,2);
      TComplex  FF22 = BWIGML(S2,MF2,GF2,M3,M1,2);
      TComplex  FSG1 = BWIGML(S1,MSG,GSG,M2,M3,0);
      TComplex  FSG2 = BWIGML(S2,MSG,GSG,M3,M1,0);
      TComplex  FF01 = BWIGML(S1,MF0,GF0,M2,M3,0);
      TComplex  FF02 = BWIGML(S2,MF0,GF0,M3,M1,0);
      
      F3PIFactor = 
	BT3*(TComplex(F34A,0.)*FRO1+TComplex(F34B,0.)*FRO2)+
	BT4*(TComplex(F34A,0.)*FRP1+TComplex(F34B,0.)*FRP2)
	-BT5*(TComplex(F35A,0.)*FF21+TComplex(F35B,0.)*FF22)
	-BT6*(TComplex(F36A,0.)*FSG1+TComplex(F36B,0.)*FSG2)
	-BT7*(TComplex(F36A,0.)*FF01+TComplex(F36B,0.)*FF02);
      
      // F3PIFactor = TComplex(0.,0.);
    }
  }
  
  if(IDK==1){
   if(IFORM == 1 || IFORM == 2 ){
     double  S1 = SA;
     double  S2 = SB;
     double  S3 = QQ-SA-SB+M1SQ+M2SQ+M3SQ;

// C it is 2pi0pi-
// C Lorentz invariants for all the contributions:
    double   F134 = -(1./3.)*((S3-M3SQ)-(S1-M1SQ));
    double   F150 =  (1./18.)*(QQ-M3SQ+S3)*(2.*M1SQ+2.*M2SQ-S3)/S3;
    double   F167 =  (2./3.);

    //C Breit Wigners for all the contributions:
    TComplex FRO1 = BWIGML(S1,MRO,GRO,M2,M3,1);
    TComplex FRP1 = BWIGML(S1,MRP,GRP,M2,M3,1);
    TComplex FRO2 = BWIGML(S2,MRO,GRO,M3,M1,1);
    TComplex FRP2 = BWIGML(S2,MRP,GRP,M3,M1,1);
    TComplex FF23 = BWIGML(S3,MF2,GF2,M1,M2,2);
    TComplex FSG3 = BWIGML(S3,MSG,GSG,M1,M2,0);
    TComplex FF03 = BWIGML(S3,MF0,GF0,M1,M2,0);

    F3PIFactor = BT1*FRO1+BT2*FRP1+
      BT3*TComplex(F134,0.)*FRO2+BT4*TComplex(F134,0.)*FRP2+
      BT5*TComplex(F150,0.)*FF23+
      BT6*TComplex(F167,0.)*FSG3+BT7*TComplex(F167,0.)*FF03;
   }
   else if (IFORM == 3 ){
     double   S3 = SA;
     double   S1 = SB;
     double   S2 = QQ-SA-SB+M1SQ+M2SQ+M3SQ;
      
     double F34A = (1./3.)*((S2-M2SQ)-(S3-M3SQ));
     double F34B = (1./3.)*((S3-M3SQ)-(S1-M1SQ));
     double F35  =-(1./2.)*((S1-M1SQ)-(S2-M2SQ));

     //C Breit Wigners for all the contributions:
     TComplex FRO1 = BWIGML(S1,MRO,GRO,M2,M3,1);
     TComplex FRP1 = BWIGML(S1,MRP,GRP,M2,M3,1);
     TComplex FRO2 = BWIGML(S2,MRO,GRO,M3,M1,1);
     TComplex FRP2 = BWIGML(S2,MRP,GRP,M3,M1,1);
     TComplex FF23 = BWIGML(S3,MF2,GF2,M1,M2,2);

     F3PIFactor = 
       BT3*(TComplex(F34A,0.)*FRO1+TComplex(F34B,0.)*FRO2)+
       BT4*(TComplex(F34A,0.)*FRP1+TComplex(F34B,0.)*FRP2)+
       BT5*TComplex(F35,0.)*FF23;
     
   }
  }
  TComplex FORMA1 = FA1A1P(QQ);
  return  F3PIFactor*FORMA1;
} 


TComplex
TauHelper::FA1A1P(double XMSQ){
  double  XM1 = 1.275000;
  double  XG1 =0.700 ; 
  double  XM2 = 1.461000 ;
  double  XG2 = 0.250; 
  TComplex BET = TComplex(0.00,0.);

  double GG1 = XM1*XG1/(1.3281*0.806);
  double GG2 = XM2*XG2/(1.3281*0.806);
  double XM1SQ = XM1*XM1;
  double XM2SQ = XM2*XM2;

  double GF = WGA1(XMSQ);
  double FG1 = GG1*GF;
  double FG2 = GG2*GF;
  TComplex F1 = TComplex(-XM1SQ,0.0)/TComplex(XMSQ-XM1SQ,FG1);
  TComplex F2 = TComplex(-XM2SQ,0.0)/TComplex(XMSQ-XM2SQ,FG2);
  TComplex FA1A1P = F1+BET*F2;

  return FA1A1P;
}


double 
TauHelper::WGA1(double QQ){
// C mass-dependent M*Gamma of a1 through its decays to 
// C.   [(rho-pi S-wave) + (rho-pi D-wave) + 
// C.    (f2 pi D-wave) + (f0pi S-wave)]
// C.  AND simple K*K S-wave
  double  MKST = 0.894;
  double  MK = 0.496;
  double  MK1SQ = (MKST+MK)*(MKST+MK);
  double  MK2SQ = (MKST-MK)*(MKST-MK);
  //C coupling constants squared:
  double   C3PI = 0.2384*0.2384;
  double   CKST = 4.7621*4.7621*C3PI;
// C Parameterization of numerical integral of total width of a1 to 3pi.
// C From M. Schmidtler, CBX-97-64-Update.
  double  S = QQ;
  double  WG3PIC = WGA1C(S);
  double  WG3PIN = WGA1N(S);

  //C Contribution to M*Gamma(m(3pi)^2) from S-wave K*K, if above threshold
  double  GKST = 0.0;
  if(S > MK1SQ) GKST = sqrt((S-MK1SQ)*(S-MK2SQ))/(2.*S);

  return C3PI*(WG3PIC+WG3PIN)+CKST*GKST;
}


double
TauHelper::WGA1C(double S){
  double STH,Q0,Q1,Q2,P0,P1,P2,P3,P4,G1_IM;
  Q0 =   5.80900; Q1 =  -3.00980; Q2 =   4.57920;
  P0 = -13.91400; P1 =  27.67900; P2 = -13.39300;
  P3 =   3.19240; P4 =  -0.10487; STH=0.1753;

  if(S < STH){
    G1_IM=0.0;
  }else if(S > STH && S < 0.823){
    G1_IM = Q0*   pow(S-STH,3)   *(1. + Q1*(S-STH) + Q2*pow(S-STH,2));
  }
  else{
    G1_IM = P0 + P1*S + P2*S*S+ P3*S*S*S + P4*S*S*S*S;
  }
  return G1_IM;
}

double
TauHelper::WGA1N(double S){
  double STH,Q0,Q1,Q2,P0,P1,P2,P3,P4,G1_IM;
  Q0 =   6.28450;Q1 =  -2.95950;Q2 =   4.33550;
  P0 = -15.41100;P1 =  32.08800;P2 = -17.66600;
  P3 =   4.93550;P4 =  -0.37498;STH   = 0.1676;
  if(S < STH){
    G1_IM = 0.0;
  }else if(S > STH && S < 0.823){
    G1_IM = Q0*pow(S-STH,3)*(1. + Q1*(S-STH) + Q2*pow(S-STH,2));
  }
  else{
    G1_IM = P0 + P1*S + P2*S*S+ P3*S*S*S + P4*S*S*S*S;
  }
  return G1_IM;
}



vector<double>
TauHelper::StructureFunction(){ 
  std::vector<double> out;
  TLorentzVector q1= SSPion1ZFrame_;
  TLorentzVector q2= SSPion2ZFrame_;
  TLorentzVector q3= OSPionZFrame_;
  TLorentzVector a1 = q1+q2+q3;

   TLorentzVector q1_a1rf= q1;//Boost(SSPion1ZFrame_,a1);
   TLorentzVector q2_a1rf= q2;//Boost(SSPion2ZFrame_,a1);
   TLorentzVector q3_a1rf= q3;//Boost(OSPionZFrame_,a1);
   TLorentzVector a1_a1rf= a1;//Boost(q1+q2+q3,a1);

  // TLorentzVector q1_a1rf= Boost(SSPion1ZFrame_,a1);
  // TLorentzVector q2_a1rf= Boost(SSPion2ZFrame_,a1);
  // TLorentzVector q3_a1rf= Boost(OSPionZFrame_,a1);
  // TLorentzVector a1_a1rf= Boost(q1+q2+q3,a1);


  //    TLorentzVector N = _nuLV;
  //    TLorentzVector P = _tauLV;
  double s1 = (q2+q3).M2();
  double s2 = (q1+q3).M2();
  double s3 = (q1+q2).M2();


  TLorentzVector vec1 = q2_a1rf - q3_a1rf -  a1_a1rf* (a1*(q2_a1rf-q3_a1rf)/a1_a1rf.M2());
  TLorentzVector vec2 = q3_a1rf - q1_a1rf -  a1_a1rf* (a1*(q3_a1rf-q1_a1rf)/a1_a1rf.M2());
  TLorentzVector vec3 = q1_a1rf - q2_a1rf -  a1_a1rf* (a1*(q1_a1rf-q2_a1rf)/a1_a1rf.M2());
   
  TComplex F1 = TComplex(COEF1)*F3PI(1,a1_a1rf.M2(),s1,s2);
  TComplex F2 = TComplex(COEF2)*F3PI(2,a1_a1rf.M2(),s2,s1);
  TComplex F3 = TComplex(COEF3)*F3PI(3,a1_a1rf.M2(),s3,s1);



  std::vector<TComplex> HADCUR;
  std::vector<TComplex> HADCURC;

  HADCUR.push_back(TComplex(vec1.E())*F1  + TComplex(vec2.E())*F2   +   TComplex(vec3.E())*F3 ); // energy component goes first
  HADCUR.push_back(TComplex(vec1.Px())*F1 + TComplex(vec2.Px())*F2  +   TComplex(vec3.Px())*F3 );
  HADCUR.push_back(TComplex(vec1.Py())*F1 + TComplex(vec2.Py())*F2  +   TComplex(vec3.Py())*F3 );
  HADCUR.push_back(TComplex(vec1.Pz())*F1 + TComplex(vec2.Pz())*F2  +   TComplex(vec3.Pz())*F3 );
  

  HADCURC.push_back(Conjugate(TComplex(vec1.E())*F1  + TComplex(vec2.E())*F2   +   TComplex(vec3.E())*F3) ); // energy component goes first
  HADCURC.push_back(Conjugate(TComplex(vec1.Px())*F1 + TComplex(vec2.Px())*F2  +   TComplex(vec3.Px())*F3 ));
  HADCURC.push_back(Conjugate(TComplex(vec1.Py())*F1 + TComplex(vec2.Py())*F2  +   TComplex(vec3.Py())*F3 ) );
  HADCURC.push_back(Conjugate(TComplex(vec1.Pz())*F1 + TComplex(vec2.Pz())*F2  +   TComplex(vec3.Pz())*F3) );


  TComplex Wa = HADCUR.at(1)*HADCURC.at(1) + HADCUR.at(2)*HADCURC.at(2);
  TComplex Wc = HADCUR.at(1)*HADCURC.at(1) - HADCUR.at(2)*HADCURC.at(2);
  TComplex Wd = HADCUR.at(1)*HADCURC.at(2) + HADCUR.at(2)*HADCURC.at(1);
  TComplex We = HADCUR.at(1)*HADCURC.at(2) - HADCUR.at(2)*HADCURC.at(1);
		  
  // std::cout<<" Wa  " << Wa.Re()   <<std::endl;
  // std::cout<<" Wc  " <<  Wc.Re()   <<std::endl;
  // std::cout<<" Wd  " <<  Wd.Re()   <<std::endl;
  // std::cout<<" We  " <<  We.Im()  <<std::endl;


  out.push_back(Wa.Re());
  out.push_back(Wc.Re());
  out.push_back(Wd.Re());  
  out.push_back(We.Im());


// //    HADCURC.push_back(Conjugate(TComplex(vec1.E())*F1  + TComplex(vec2.E())*F2   +   TComplex(vec3.E())*F3) ); // energy component goes first
// //    HADCURC.push_back(Conjugate(TComplex(vec1.Px())*F1 + TComplex(vec2.Px())*F2  +   TComplex(vec3.Px())*F3 ));
// //    HADCURC.push_back(Conjugate(TComplex(vec1.Py())*F1 + TComplex(vec2.Py())*F2  +   TComplex(vec3.Py())*F3 ) );
// //    HADCURC.push_back(Conjugate(TComplex(vec1.Pz())*F1 + TComplex(vec2.Pz())*F2  +   TComplex(vec3.Pz())*F3) );

// //    TLorentzVector CLV =  CLVEC(HADCUR,HADCURC,N );
// //    TLorentzVector CLA =  CLAXI(HADCUR,HADCURC,N );

// //    TComplex BWProd1 = f3(a1.M())*BreitWigner(sqrt(s2),"rho");
// //    TComplex BWProd2 = f3(a1.M())*BreitWigner(sqrt(s1),"rho");
 
// //    double omega = P*CLV - P*CLA;
   return out;
  }

vector<double>
TauHelper::StructureFunctionKuhn(){ 
  std::vector<double> out;
  TLorentzVector q1= SSPion1ZFrame_;
  TLorentzVector q2= SSPion2ZFrame_;
  TLorentzVector q3= OSPionZFrame_;
  TLorentzVector a1 = q1+q2+q3;

   TLorentzVector q1_a1rf= q1;//Boost(SSPion1ZFrame_,a1);
   TLorentzVector q2_a1rf= q2;//Boost(SSPion2ZFrame_,a1);
   TLorentzVector q3_a1rf= q3;//Boost(OSPionZFrame_,a1);
   TLorentzVector a1_a1rf= a1;//Boost(q1+q2+q3,a1);

  // TLorentzVector q1_a1rf= Boost(SSPion1ZFrame_,a1);
  // TLorentzVector q2_a1rf= Boost(SSPion2ZFrame_,a1);
  // TLorentzVector q3_a1rf= Boost(OSPionZFrame_,a1);
  // TLorentzVector a1_a1rf= Boost(q1+q2+q3,a1);


  //    TLorentzVector N = _nuLV;
  //    TLorentzVector P = _tauLV;
  double s1 = (q2+q3).M2();
  double s2 = (q1+q3).M2();
  double s3 = (q1+q2).M2();


  TLorentzVector vec1 = q1_a1rf - q3_a1rf -  a1_a1rf* (a1*(q1_a1rf-q3_a1rf)/a1_a1rf.M2());
  TLorentzVector vec2 = q2_a1rf - q3_a1rf -  a1_a1rf* (a1*(q2_a1rf-q3_a1rf)/a1_a1rf.M2());
  //  TLorentzVector vec3 = q1_a1rf - q2_a1rf -  a1_a1rf* (a1*(q1_a1rf-q2_a1rf)/a1_a1rf.M2());
   
  // TComplex F1 = TComplex(COEF1)*F3PI(1,a1_a1rf.M2(),s1,s2);
  // TComplex F2 = TComplex(COEF2)*F3PI(2,a1_a1rf.M2(),s2,s1);
  // TComplex F3 = TComplex(COEF3)*F3PI(3,a1_a1rf.M2(),s3,s1);



  std::vector<TComplex> HADCUR;
  std::vector<TComplex> HADCURC;

  HADCUR.push_back(TComplex(vec1.E()) *F1() + TComplex(vec2.E())*F2()  ); // energy component goes first
  HADCUR.push_back(TComplex(vec1.Px())*F1() + TComplex(vec2.Px())*F2()  );
  HADCUR.push_back(TComplex(vec1.Py())*F1() + TComplex(vec2.Py())*F2()  );
  HADCUR.push_back(TComplex(vec1.Pz())*F1() + TComplex(vec2.Pz())*F2()  );
  //  std::cout<<" Kuhn  F1  " << F1()<< std::endl;

  HADCURC.push_back(Conjugate(TComplex(vec1.E())*F1()  + TComplex(vec2.E())*F2()   ) ); // energy component goes first
  HADCURC.push_back(Conjugate(TComplex(vec1.Px())*F1() + TComplex(vec2.Px())*F2()  ) );
  HADCURC.push_back(Conjugate(TComplex(vec1.Py())*F1() + TComplex(vec2.Py())*F2()  ) );
  HADCURC.push_back(Conjugate(TComplex(vec1.Pz())*F1() + TComplex(vec2.Pz())*F2()  ) );


  TComplex Wa = HADCUR.at(1)*HADCURC.at(1) + HADCUR.at(2)*HADCURC.at(2);
  TComplex Wc = HADCUR.at(1)*HADCURC.at(1) - HADCUR.at(2)*HADCURC.at(2);
  TComplex Wd = HADCUR.at(1)*HADCURC.at(2) + HADCUR.at(2)*HADCURC.at(1);
  TComplex We = HADCUR.at(1)*HADCURC.at(2) - HADCUR.at(2)*HADCURC.at(1);
		  
  // std::cout<<" Wa  " << Wa.Re()   <<std::endl;
  // std::cout<<" Wc  " <<  Wc.Re()   <<std::endl;
  // std::cout<<" Wd  " <<  Wd.Re()   <<std::endl;
  // std::cout<<" We  " <<  We.Im()  <<std::endl;

  // _s12 = SSPion1ZFrame_  + SSPion2ZFrame_;
  // _s13 = SSPion1ZFrame_  + OSPionZFrame_;
  // _s23 = SSPion2ZFrame_  + OSPionZFrame_;
  // _s1  =  _s23.M2(); 
  // _s2  =  _s13.M2();
  // _s3  =  _s12.M2();


  out.push_back(Wa.Re());
  out.push_back(Wc.Re());
  out.push_back(Wd.Re());  
  out.push_back(We.Im());


// //    HADCURC.push_back(Conjugate(TComplex(vec1.E())*F1  + TComplex(vec2.E())*F2   +   TComplex(vec3.E())*F3) ); // energy component goes first
// //    HADCURC.push_back(Conjugate(TComplex(vec1.Px())*F1 + TComplex(vec2.Px())*F2  +   TComplex(vec3.Px())*F3 ));
// //    HADCURC.push_back(Conjugate(TComplex(vec1.Py())*F1 + TComplex(vec2.Py())*F2  +   TComplex(vec3.Py())*F3 ) );
// //    HADCURC.push_back(Conjugate(TComplex(vec1.Pz())*F1 + TComplex(vec2.Pz())*F2  +   TComplex(vec3.Pz())*F3) );

// //    TLorentzVector CLV =  CLVEC(HADCUR,HADCURC,N );
// //    TLorentzVector CLA =  CLAXI(HADCUR,HADCURC,N );

// //    TComplex BWProd1 = f3(a1.M())*BreitWigner(sqrt(s2),"rho");
// //    TComplex BWProd2 = f3(a1.M())*BreitWigner(sqrt(s1),"rho");
 
// //    double omega = P*CLV - P*CLA;
   return out;
  }



TVector3
TauHelper::Rotate(TVector3 LVec, TVector3 Rot){
  TVector3 vec = LVec;
  vec.RotateZ(0.5*TMath::Pi() - Rot.Phi());  
  vec.RotateX(Rot.Theta());
  return vec;
}
