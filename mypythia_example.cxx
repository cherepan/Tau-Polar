/**
 * Example of use of tauola C++ interfate. Pythia events are
 * generated with a stable tau. Taus are subseuently decay via
 * tauola, plots of polatization observables for tau-> mununu and tau-> pinu
 * are produced.
 *
 * @author Vladimir Cherepanov
 * @date  01 MAY 2017
 */ 
  
#include "Tauola/Log.h"
#include "Tauola/Tauola.h"
#include "Tauola/TauolaHepMCEvent.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "UserCodes/MultiplyNumbers.h"
#include "UserCodes/a1Helper.h"

//pythia header files
#ifdef PYTHIA8180_OR_LATER
#include "Pythia8/Pythia.h"
#include "Pythia8/Pythia8ToHepMC.h"
#else 
#include "Pythia.h"
#include "HepMCInterface.h"
#endif

//MC-TESTER header files
#include "Generate.h"
#include "HepMCEvent.H"
#include "Setup.H"
 
#include "tauola_print_parameters.h"
using namespace std;
using namespace Pythia8; 
using namespace Tauolapp;

int NumberOfEvents = 1000000; 
int EventsToCheck=10;

// elementary test of HepMC typically executed before
// detector simulation based on http://home.fnal.gov/~mrenna/HCPSS/HCPSShepmc.html
// similar test was performed in Fortran
// we perform it before and after Tauola (for the first several events)
void checkMomentumConservationInEvent(HepMC::GenEvent *evt)
{
 	cout<<"List of stable particles: "<<endl;

	double px=0.0,py=0.0,pz=0.0,e=0.0;

	int barcodetau1vertex(0),barcodetau2vertex(0);
	for ( HepMC::GenEvent::particle_const_iterator p = evt->particles_begin();  p != evt->particles_end(); ++p )
	  {
	    if( (*p)->status() == 1 )
	      {
		HepMC::FourVector m = (*p)->momentum();
		HepMC::GenVertex *productProductionvertex = (*p)->production_vertex();
		px+=m.px();
		py+=m.py();
		pz+=m.pz();
		e +=m.e();
	      }
	  }
	cout.precision(6);
	cout.setf(ios_base::floatfield);
}

void redMinus(TauolaParticle *minus)
{
  //   
  // this method can be used to redefine branching ratios in decay of tau-
  // either generally, or specific to  tau- with pointer *minus.
  //
  // Pointer *minus can be used to define properties of decays for taus
  // at specific point(s) in the event tree. Example: 
  // vector<TauolaParticle*> x=minus->getMothers();
  // and define special versions depending on x. 
  //
  // Any combination of methods
  // Tauola::setTauBr(int mode, double br);
  // Tauola::setTaukle(double bra1, double brk0, double brk0b,double brks);
  // can be called here 


  for(unsigned int dec=0; dec <23; dec++){
     double br =0.0; 
    if(dec==2) br=0.99;
     Tauola::setTauBr(dec, br);
   } 

     // Tauola::setTauBr(0, 0.0);     Tauola::setTauBr(1, 0.0);      Tauola::setTauBr(2, 1.0);      Tauola::setTauBr(3, 0.0);      Tauola::setTauBr(4, 0.0);
     // Tauola::setTauBr(5, 0.0);     Tauola::setTauBr(6, 0.0);      Tauola::setTauBr(7, 0.0);      Tauola::setTauBr(8, 0.0);      Tauola::setTauBr(9, 0.0);
     // Tauola::setTauBr(10, 0.0);   Tauola::setTauBr(11, 0.0);    Tauola::setTauBr(12, 0.0);    Tauola::setTauBr(13, 0.0);    Tauola::setTauBr(14, 0.0);
     // Tauola::setTauBr(15, 0.0);   Tauola::setTauBr(16, 0.0);    Tauola::setTauBr(17, 0.0);    Tauola::setTauBr(18, 0.0);    Tauola::setTauBr(19, 0.0);
     // Tauola::setTauBr(20, 0.0);   Tauola::setTauBr(21, 0.0);    Tauola::setTauBr(22, 0.0);


}

void redPlus(TauolaParticle *plus)
{
  //   
  // this method can be used to redefine branching ratios in decay of tau+
  // either generally, or specific to  tau+ with pointer *plus.
  //
  // Pointer *plus can be used to define properties of decays for tau
  // at specific point(s) in the event tree. Example: 
  // vector<TauolaParticle*> x=plus->getMothers();
  // and define special versions depending on x. 
  //
  // Any combination of methods
  // Tauola::setTauBr(int mode, double br);
  // Tauola::setTaukle(double bra1, double brk0, double brk0b,double brks);
  // can be called here 
  for(unsigned int dec=0; dec <23; dec++){
     double br =0.0;
     //    if(dec==3 || dec ==4 || dec ==5) br=0.33;
    if(dec ==5) br=0.99;
     Tauola::setTauBr(dec, br);
   }


  // Tauola::setTauBr(0, 0.0);     Tauola::setTauBr(1, 0.0);      Tauola::setTauBr(2, 0.0);      Tauola::setTauBr(3, 0.5);      Tauola::setTauBr(4, 0.5);
  // Tauola::setTauBr(5, 0.0);     Tauola::setTauBr(6, 0.0);      Tauola::setTauBr(7, 0.0);      Tauola::setTauBr(8, 0.0);      Tauola::setTauBr(9, 0.0);
  // Tauola::setTauBr(10, 0.0);   Tauola::setTauBr(11, 0.0);    Tauola::setTauBr(12, 0.0);    Tauola::setTauBr(13, 0.0);    Tauola::setTauBr(14, 0.0);
  // Tauola::setTauBr(15, 0.0);   Tauola::setTauBr(16, 0.0);    Tauola::setTauBr(17, 0.0);    Tauola::setTauBr(18, 0.0);    Tauola::setTauBr(19, 0.0);
  // Tauola::setTauBr(20, 0.0);   Tauola::setTauBr(21, 0.0);    Tauola::setTauBr(22, 0.0);


}
void SortPions(    std::vector<HepMC::GenParticle > pionsvec)
{

  int npim(0),npip(0);
  int    OSMCPionIndex(0);
  int    SSMCPion1Index(0);
  int    SSMCPion2Index(0);

    HepMC::GenParticle os;
    HepMC::GenParticle ss1;
    HepMC::GenParticle ss2;
    for(int i=0; i<pionsvec.size(); i++){
      if( pionsvec.at(i).pdg_id()== 211) npip++;
      if( pionsvec.at(i).pdg_id()==-211) npim++;
    }
    if(npip == 1 && npim==2){
      int nss=0;
      for(int i=0; i<pionsvec.size(); i++){
	if(pionsvec.at(i).pdg_id()== 211){
	  OSMCPionIndex=i;
	}
	if(pionsvec.at(i).pdg_id()==-211 && nss ==0){
	  nss++;
	  SSMCPion1Index=i;
	}
	if(pionsvec.at(i).pdg_id()==-211 &&  nss == 1){
	  SSMCPion2Index=i;
	}
      }
    }
    if( npip== 2 && npim==1){
      int nss=0;
      for(int i=0; i<pionsvec.size(); i++){
	if(pionsvec.at(i).pdg_id()== -211){
	  
	  OSMCPionIndex=i;
	}
	if(pionsvec.at(i).pdg_id()==211 && nss ==0){
	  nss++;
	  SSMCPion1Index=i;
	}
	if(pionsvec.at(i).pdg_id()==211 &&  nss == 1){
	  SSMCPion2Index=i;
	}
      }
    }
    os=pionsvec.at(OSMCPionIndex);
    ss1=pionsvec.at(SSMCPion1Index);
    ss2=pionsvec.at(SSMCPion2Index);
    
    
    // std::cout<<" os.pdgId()  "<<os.pdg_id() <<std::endl;
    // std::cout<<" ss.pdgId()  "<<ss1.pdg_id() <<std::endl;
    // std::cout<<" ss.pdgId()  "<<ss2.pdg_id() <<std::endl;
    
    pionsvec.clear();
    pionsvec.push_back(os);
    pionsvec.push_back(ss1);
    pionsvec.push_back(ss2);
}

int main(int argc,char **argv){

  Log::SummaryAtExit();

  // Initialization of pythia
  Pythia pythia;
  Event& event = pythia.event;
  TFile *file = new TFile("output.root","RECREATE");

  TH1F *pi_plus= new TH1F("pi_plus","#pi  plus",50,0,1);
  TH1F *pi_minus= new TH1F("pi_minus","#pi minus",50,0,1);

  TH1F *rho_plus= new TH1F("rho_plus","#rho  plus",50,-1,1);
  TH1F *rho_minus= new TH1F("rho_minus","#rho minus",50,-1,1);


  TH1F *piom_plus= new TH1F("piom_plus","#pi  plus",50,-1,1);
  TH1F *piom_minus= new TH1F("piom_minus","#pi minus",50,-1,1);


  TH1F *mu_plus= new TH1F("mu_plus","#mu plus",50,0,1);
  TH1F *mu_minus= new TH1F("mu_minus","#mu  minus",50,0,1);
  
  TH1F *om_plus= new TH1F("om_plus","#omega #pi #mu  plus",50,-1,1);
  TH1F *om_minus= new TH1F("om_minus","#omega  #pi #mu  minus",50,-1,1);
 
  TH1F *Omegapipi_plus= new TH1F("Omegapipi_plus","#omega #pi#pi plus",50,-1,1);
  TH1F *Omegapipi_minus= new TH1F("Omegapipi_minus","#omega #pi#pi  minus",50,-1,1);
  
  TH1F *Omegapirho_plus= new TH1F("Omegapirho_plus","#omega #pi#rho plus",50,-1,1);
  TH1F *Omegapirho_minus= new TH1F("Omegapirho_minus","#omega #pi#rho  minus",50,-1,1);
  

 
 TH1F *ommurho_plus= new TH1F("ommurho_plus","#omega #mu #rho plus",50,-1,1);
 TH1F *ommurho_minus= new TH1F("ommurho_minus","#omega  #mu #rho minus",50,-1,1);


 TH1F *omega_a1_plus= new TH1F("omega_a1_plus","#omega a1",50,-2,2);
 TH1F *omega_a1_minus= new TH1F("omega_a1_minus","#omega a1",50,-2,2);

 TH1F *omegabar_a1_plus= new TH1F("omegabar_a1_plus","#omega a1",50,-2,2);
 TH1F *omegabar_a1_minus= new TH1F("omegabar_a1_minus","#omega a1",50,-2,2);
 

 TH1F *TRFomegabar_a1_plus= new TH1F("TRFomegabar_a1_plus","#omega a1",50,-2,2);
 TH1F *TRFomegabar_a1_minus= new TH1F("TRFomegabar_a1_minus","#omega a1",50,-2,2);
 


 TH2F *cosbetacostheta_plus= new TH2F("cosbetacostheta_plus","cos#beta  cos#theta",50,-1,1,50,-1,1);
 TH2F *cosbetacostheta_minus= new TH2F("cosbetacostheta_minus","cos#beta  cos#theta",50,-1,1,50,-1,1);

 TH2F *TRFcosbetacostheta_plus= new TH2F("TRFcosbetacostheta_plus","cos#beta  cos#theta",50,-1,1,50,-1,1);
 TH2F *TRFcosbetacostheta_minus= new TH2F("TRFcosbetacostheta_minus","cos#beta  cos#theta",50,-1,1,50,-1,1);


  // Pythia8 HepMC interface depends on Pythia8 version
#ifdef PYTHIA8180_OR_LATER
  HepMC::Pythia8ToHepMC ToHepMC;
#else
  HepMC::I_Pythia8 ToHepMC;
#endif

  //pythia.readString("HadronLevel:all = off");
  pythia.readString("HadronLevel:Hadronize = off");
  pythia.readString("SpaceShower:QEDshowerByL = off");
  pythia.readString("SpaceShower:QEDshowerByQ = off");
  pythia.readString("PartonLevel:ISR = off");
  pythia.readString("PartonLevel:FSR = off");

  // Tauola is currently set to undecay taus. Otherwise, uncomment this line.
  // Uncommenting it will speed up the test significantly
  //  pythia.particleData.readString("15:mayDecay = off");

  pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
  pythia.readString("23:onMode = off"); 
  pythia.readString("23:onIfAny = 15");
  //  pythia.readString("23:onIfMatch = 15 -15");

  pythia.init( 11, -11, 92.);          //electron positron collisions

  // Set up Tauola

  // Set Tauola decay mode (if needed)
  //Tauola::setSameParticleDecayMode(2);     //19 and 22 contains K0 
    //   Tauola::setOppositeParticleDecayMode(3); // 20 contains eta

  // Set Higgs scalar-pseudoscalar mixing angle
  //  Tauola::setHiggsScalarPseudoscalarMixingAngle(0.7853);
  //  Tauola::setHiggsScalarPseudoscalarPDG(25);

  Tauola::initialize();

  tauola_print_parameters(); // Prints TAUOLA  parameters (residing inside its library): e.g. to test user interface

  // Our default units are GEV and MM, that will be outcome  units after TAUOLA
  // if HepMC unit variables  are correctly set. 
  // with the following coice you can fix the units for final outcome:
  //  Tauola::setUnits(Tauola::GEV,Tauola::MM); 
  //  Tauola::setUnits(Tauola::MEV,Tauola::CM); 

  // Other usefull settings:
  //  Tauola::setEtaK0sPi(0,0,0);  // switches to decay eta K0_S and pi0 1/0 on/off. 
  //  Tauola::setTauLifetime(0.0); //new tau lifetime in mm
    Tauola::spin_correlation.setAll(true);

  //  Log::LogDebug(true);

    Tauola::setRedefineTauMinus(redMinus);  // activates execution of routine redMinus in TAUOLA interface
    Tauola::setRedefineTauPlus(redPlus);    // activates execution of routine redPlus  in TAUOLA interface

  MC_Initialize();

  // Begin event loop. Generate event.
  for (int iEvent = 0; iEvent < NumberOfEvents; ++iEvent){

    if(iEvent%1000==0) Log::Info()<<"Event: "<<iEvent<<endl;
    if (!pythia.next()) continue;

    // Convert event record to HepMC
    HepMC::GenEvent * HepMCEvt = new HepMC::GenEvent();

    // Conversion needed if HepMC use different momentum units
    // than Pythia. However, requires HepMC 2.04 or higher.
    HepMCEvt->use_units(HepMC::Units::GEV,HepMC::Units::MM);

    ToHepMC.fill_next_event(event, HepMCEvt);

    if(iEvent<EventsToCheck)
    {
      cout<<"                                          "<<endl;
      cout<<"Momentum conservation chceck BEFORE/AFTER Tauola"<<endl;
      checkMomentumConservationInEvent(HepMCEvt);
    }

    // Run TAUOLA on the event
    TauolaHepMCEvent * t_event = new TauolaHepMCEvent(HepMCEvt);

    // Since we let Pythia decay taus, we have to undecay them first.
    t_event->undecayTaus();
    t_event->decayTaus();
    delete t_event; 

    if(iEvent<EventsToCheck)
    {
      checkMomentumConservationInEvent(HepMCEvt);
    }

    int JAK1(0); int SubJAK1(0);
    HepMC::GenParticle *FirstTau;
    std::vector<HepMC::GenParticle > FirstTauProducts;
    int JAK2(0); int SubJAK2(0);
    HepMC::GenParticle *SecondTau;
    std::vector<HepMC::GenParticle > SecondTauProducts;
    std::vector<HepMC::GenParticle > A1Pions;
    std::vector<HepMC::GenParticle > SortA1Pions;  //os, ss1, ss2
    for ( HepMC::GenEvent::particle_const_iterator p =HepMCEvt->particles_begin();  p != HepMCEvt->particles_end(); ++p ){  
      if((*p)->pdg_id()==15){
	FirstTau = *p;
	for ( HepMC::GenEvent::particle_const_iterator d =HepMCEvt->particles_begin();  d != HepMCEvt->particles_end(); ++d ){  
	  if((*d)->pdg_id()!=15){
	    if((*p)->end_vertex() == (*d)->production_vertex()){
	      FirstTauProducts.push_back(**d);
	      if(abs((*d)->pdg_id()) ==  12) {JAK1 =1; 
	      }else if(abs((*d)->pdg_id()) ==  14){ JAK1=2;
	      }else if(abs((*d)->pdg_id())==  211){ JAK1 = 3;}
	      if( abs((*d)->pdg_id())==213 ){
		JAK1 = 4;
		for ( HepMC::GenEvent::particle_const_iterator dd =HepMCEvt->particles_begin();  dd != HepMCEvt->particles_end(); ++dd ){  
		  if( abs((*dd)->pdg_id())!=213  ){
		    if((*d)->end_vertex() == (*dd)->production_vertex()){
		      FirstTauProducts.push_back(**dd);
		    }
		  }
		}
	      }
		if( abs((*d)->pdg_id())==20213 ){
		  JAK1 = 5; int npi(0);
		  for ( HepMC::GenEvent::particle_const_iterator dd =HepMCEvt->particles_begin();  dd != HepMCEvt->particles_end(); ++dd ){  
		    if( abs((*dd)->pdg_id())!=20213  ){
		      if((*d)->end_vertex() == (*dd)->production_vertex()){
		       
			FirstTauProducts.push_back(**dd);
			if(abs((*dd)->pdg_id())==  211)   {
			  A1Pions.push_back(**dd); npi++;
			}
		      }
		    }
		  }
		  if(npi==3) SubJAK1=51; else SubJAK1=52;
		}
	    }
	  } 
	}
      }
    
      
      if((*p)->pdg_id()==-15){
	SecondTau = *p;
	for ( HepMC::GenEvent::particle_const_iterator d =HepMCEvt->particles_begin();  d != HepMCEvt->particles_end(); ++d )
	  {  
	    if((*d)->pdg_id()!=-15){
	      if((*p)->end_vertex() == (*d)->production_vertex()){
		SecondTauProducts.push_back(**d);
		if(abs((*d)->pdg_id()) ==  12) {JAK2 =1; 
		}else if(abs((*d)->pdg_id()) ==  14){ JAK2=2;
		}else if(abs((*d)->pdg_id())==  211){ JAK2 = 3;}

		if( abs((*d)->pdg_id())==213 ){
		  JAK2 = 4;
		  for ( HepMC::GenEvent::particle_const_iterator dd =HepMCEvt->particles_begin();  dd != HepMCEvt->particles_end(); ++dd ){  
		    if( abs((*dd)->pdg_id())!=213  ){
		      if((*d)->end_vertex() == (*dd)->production_vertex()){
			SecondTauProducts.push_back(**dd);
		      }
		    }
		  }
		}
		if( abs((*d)->pdg_id())==20213 ){
		  JAK2 = 5; int npi(0);
		  for ( HepMC::GenEvent::particle_const_iterator dd =HepMCEvt->particles_begin();  dd != HepMCEvt->particles_end(); ++dd ){  
		    if( abs((*dd)->pdg_id())!=20213  ){
		      if((*d)->end_vertex() == (*dd)->production_vertex()){
			SecondTauProducts.push_back(**dd);
			if(abs((*dd)->pdg_id())==  211){
			  A1Pions.push_back(**dd); npi++;
			}
		      }
		    }
		  }
		  if(npi==3) SubJAK2=51; else SubJAK2=52;
		}
		
	      }
	    } 
	  }
      }
    }
  
    //    cout<<" JAK1 " << JAK1 << "  JAK2  "<< JAK2 <<endl;
 
    TLorentzVector tau1(0,0,0,0);
    TLorentzVector mu1(0,0,0,0);
    TLorentzVector numu1(0,0,0,0);
    TLorentzVector nutau1(0,0,0,0);
    TLorentzVector nutau2(0,0,0,0);
    TLorentzVector pi2(0,0,0,0);
    TLorentzVector pi1(0,0,0,0);
    TLorentzVector tau2(0,0,0,0);

    TLorentzVector rhopi2(0,0,0,0);
    TLorentzVector rhopi02(0,0,0,0);

    TLorentzVector a1ospi(0,0,0,0);
    TLorentzVector a1ss1pi(0,0,0,0);
    TLorentzVector a1ss2pi(0,0,0,0);
    TLorentzVector a1(0,0,0,0);

    tau1.SetPxPyPzE(FirstTau->momentum().px(), FirstTau->momentum().py(), FirstTau->momentum().pz(), FirstTau->momentum().e());
    tau2.SetPxPyPzE(SecondTau->momentum().px(), SecondTau->momentum().py(), SecondTau->momentum().pz(), SecondTau->momentum().e());
    
    for(std::vector<HepMC::GenParticle>::const_iterator a = FirstTauProducts.begin(); a!=FirstTauProducts.end(); ++a){
      if(JAK1==2){
	if(abs(a->pdg_id())==13)mu1.SetPxPyPzE(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e());
	if(abs(a->pdg_id())==14)numu1.SetPxPyPzE(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e());
	if(abs(a->pdg_id())==16)nutau1.SetPxPyPzE(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e());
      }

      if(JAK1==3){
	//	if(abs(a->pdg_id())==16)nutau1.SetPxPyPzE(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e());
	if(abs(a->pdg_id())==211)pi1.SetPxPyPzE(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e());
     }





    }

    for(std::vector<HepMC::GenParticle>::const_iterator a = SecondTauProducts.begin(); a!=SecondTauProducts.end(); ++a){
      if(JAK2==3){
	if(abs(a->pdg_id())==16)nutau2.SetPxPyPzE(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e());
	if(abs(a->pdg_id())==211)pi2.SetPxPyPzE(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e());
      }
    
      if(JAK2==4){
	if(abs(a->pdg_id())==111)rhopi02.SetPxPyPzE(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e());
	if(abs(a->pdg_id())==211)rhopi2.SetPxPyPzE(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e());
      }
      if(JAK2==5){
	if(abs(a->pdg_id())==20213)a1.SetPxPyPzE(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e());
      }
    }
   
    if(JAK2==5 && SubJAK2==51){
     SortPions(A1Pions);
     a1ospi.SetPxPyPzE(A1Pions.at(0).momentum().px(), A1Pions.at(0).momentum().py(), A1Pions.at(0).momentum().pz(), A1Pions.at(0).momentum().e());
     a1ss1pi.SetPxPyPzE(A1Pions.at(1).momentum().px(), A1Pions.at(1).momentum().py(), A1Pions.at(1).momentum().pz(), A1Pions.at(1).momentum().e());
     a1ss2pi.SetPxPyPzE(A1Pions.at(2).momentum().px(), A1Pions.at(2).momentum().py(), A1Pions.at(2).momentum().pz(), A1Pions.at(2).momentum().e());
     // std::cout<<"---------"<<std::endl;
     // tau2.Print();
     // a1.Print();
     // a1ospi.Print();
     // a1ss1pi.Print();
     // a1ss2pi.Print();
     vector<TLorentzVector> particles;
     particles.push_back(tau1);
     particles.push_back(a1ospi);
     particles.push_back(a1ss1pi);
     particles.push_back(a1ss2pi);

  // particles.at(0).Print();
  // particles.at(1).Print();
  // particles.at(2).Print();
  // particles.at(3).Print();

     // a1Helper Helper(particles, a1ospi+a1ss1pi+a1ss2pi);
     //     cout<<     Helper.WA()*Helper.WA()<< "  = "<<Helper.WC()*Helper.WC()+  Helper.WD()*Helper.WD()+ Helper.WE()*Helper.WE()<<endl;	 
     



    }



    // int barcodetau1vertex(0),barcodetau2vertex(0);
    // for ( HepMC::GenEvent::particle_const_iterator p =HepMCEvt->particles_begin();  p != HepMCEvt->particles_end(); ++p )
    //   {

    // 	if((*p)->pdg_id()==15){
    // 	  HepMC::GenVertex * tau1endvertex =(*p)->end_vertex(); 
    // 	  barcodetau1vertex = tau1endvertex->barcode();  
    // 	  tau1.SetPxPyPzE((*p)->momentum().px(), (*p)->momentum().py(), (*p)->momentum().pz(), (*p)->momentum().e());
    // 	}
    // 	if((*p)->pdg_id()==-15){
    // 	  HepMC::GenVertex * tau2endvertex =(*p)->end_vertex(); 
    // 	  tau2.SetPxPyPzE((*p)->momentum().px(), (*p)->momentum().py(), (*p)->momentum().pz(), (*p)->momentum().e());
    // 	  barcodetau2vertex = tau2endvertex->barcode();  
    // 	}
    // 	if( (*p)->status() == 1 )
    // 	  {
    // 	    HepMC::FourVector m = (*p)->momentum();
    // 	    HepMC::GenVertex *productProductionvertex = (*p)->production_vertex();
    // 	    if(productProductionvertex->barcode()==barcodetau1vertex ){
    // 	      if(abs((*p)->pdg_id())==13)mu1.SetPxPyPzE((*p)->momentum().px(), (*p)->momentum().py(), (*p)->momentum().pz(), (*p)->momentum().e());
    // 	      if(abs((*p)->pdg_id())==14)numu1.SetPxPyPzE((*p)->momentum().px(), (*p)->momentum().py(), (*p)->momentum().pz(), (*p)->momentum().e());
    // 	      if(abs((*p)->pdg_id())==16)nutau1.SetPxPyPzE((*p)->momentum().px(), (*p)->momentum().py(), (*p)->momentum().pz(), (*p)->momentum().e());
    // 	    }
    // 	    if(productProductionvertex->barcode()==barcodetau2vertex ){
    // 	      if(abs((*p)->pdg_id())==16)nutau2.SetPxPyPzE((*p)->momentum().px(), (*p)->momentum().py(), (*p)->momentum().pz(), (*p)->momentum().e());
    // 	      if(abs((*p)->pdg_id())==211)pi2.SetPxPyPzE((*p)->momentum().px(), (*p)->momentum().py(), (*p)->momentum().pz(), (*p)->momentum().e());
    // 	    }
    // 	  }
    //   }
    Log::Debug(5)<<"helicites =  "<<Tauola::getHelPlus()<<" "<<Tauola::getHelMinus()
                 <<" electroweak wt= "<<Tauola::getEWwt()<<endl;

    bool HelPlus=false;
    bool HelMinus=false;
    if(Tauola::getHelPlus() ==1 )HelPlus=true;;
    if(Tauola::getHelPlus() ==-1)HelMinus=true;;

    double x1p(0),x2p(0),x3p(0);
    double x1m(0),x2m(0),x3m(0);
    double omp(0),omm(0),OmegaPMuPi(0), OmegaMMuPi(0),OmegaPMuRho(0), OmegaMMuRho(0);
    double ompipi_plus(0), ompipi_minus(0);
    double ompirho_plus(0), ompirho_minus(0);
    double x1pim(0), x1pip(0);
    //   if(Tauola::getHelPlus() ==1){
    if(HelPlus){
      if(tau2.E()!=0 && pi2.E()!=0) {x1p=2*pi2.E()/tau2.E() - 1; pi_plus->Fill(pi2.E()/tau2.E()); piom_plus->Fill(x1p);}
      if(tau1.E()!=0 && mu1.E()!=0){x2p =2*mu1.E()/tau1.E()-1;  mu_plus->Fill(mu1.E()/tau1.E());}

      if(tau1.E()!=0 && pi1.E()!=0){x1pip =2*pi1.E()/tau1.E()-1;}
    

      if(rhopi2.E()!=0 && rhopi02.E()!=0){x3p = (rhopi2.E() -   rhopi02.E())/(rhopi2.E() +   rhopi02.E());  rho_plus->Fill(x3p); }

      if(x1p!=0 && x2p!=0){OmegaPMuPi = (x1p+x2p)/(1+x1p*x2p);om_plus->Fill(OmegaPMuPi);}
      if(x2p!=0 && x3p!=0){OmegaPMuRho = (x2p+x3p)/(1+x2p*x3p);ommurho_plus->Fill(OmegaPMuRho);}
     
      if(x1p!=0 && x1pip!=0){ompipi_plus = (x1p+x1pip)/(1+x1p*x1pip);Omegapipi_plus->Fill(ompipi_plus);}
      if(x1pip!=0 && x3p!=0){ompirho_plus = (x1pip+x3p)/(1+x1pip*x3p);Omegapirho_plus->Fill(ompirho_plus);}

      if(JAK2==5 && SubJAK2==51){
	vector<TLorentzVector> particles;
	particles.push_back(tau1);
	particles.push_back(a1ospi);
	particles.push_back(a1ss1pi);
	particles.push_back(a1ss2pi);
	a1Helper Helper(particles, a1ospi+a1ss1pi+a1ss2pi);
 	if(Helper.getA1omega()> 0 or Helper.getA1omega() < 0)	omega_a1_plus->Fill(Helper.getA1omega());


	if(Helper.getg()!=0)omegabar_a1_plus->Fill(Helper.vgetA1omega("bar"));
	if(Helper.TRF_vgetA1omega()>0 or Helper.TRF_vgetA1omega()<=0)TRFomegabar_a1_plus->Fill(Helper.TRF_vgetA1omega());
	// std::cout<<" Plus:   scalar:     "<< Helper.vgetfscalar() << std::endl;
	// std::cout<<"  Mminus:   scalar  g:     "<< Helper.vgetgscalar() << std::endl;
	 // std::cout<<"Plus  :   scalar  bar  f:     "<< Helper.vgetfscalar("bar"); 	std::cout<<" Plus :   scalar    f:     "<< Helper.vgetfscalar() << std::endl;
	 // std::cout<<"Plus  :   scalar  bar g:     "<< Helper.vgetgscalar("bar");	 std::cout<<"Plus  :   scalar   g:     "<< Helper.vgetgscalar()<<std::endl;


	//	std::cout<<"HelPlus:   costhetaLF()   "<< Helper.costhetaLF() <<"   get omega   "<< Helper.getA1omega()<<  "  omega bar   " << Helper.getA1omegaBar()<<std::endl;
	//	std::cout<<"HelPlus:   costhetaLF()   "<< Helper.costhetaLF() <<"   get omega   "<< Helper.getA1omega()<<  "  omega v   " << Helper.vgetA1omega(1)<<std::endl;
	 std::cout<<"Plus:   --------------   " << "    costheta    " << Helper.costhetaLF()  <<"   get omega   "<< Helper.getA1omega()<<  "  omega bar   " << Helper.vgetA1omega("bar")<<  "    omega bar scalar    " << Helper.vgetA1omegascalar("bar")<<"   TRF  Omega   " <<  Helper.TRF_vgetA1omega()<<std::endl;
	// std::cout<<Helper.WA()*Helper.WA() << "   =====    "<<   Helper.WC()*Helper.WC() + Helper.WD()*Helper.WD() + Helper.WE()*Helper.WE() <<std::endl;

	 if(Helper.TRF_vgetg() > Helper.TRF_vgetf() ){
	   std::cout<<"   g > f !!!!!!!!!!!";  std::cout<<Helper.WA()*Helper.WA() << "      equaolyt    "<<   Helper.WC()*Helper.WC() + Helper.WD()*Helper.WD() + Helper.WE()*Helper.WE() <<std::endl;

	   std::cout<<" TRF_vgetg()  "<<  Helper.TRF_vgetg() << " TRF_vgetf()   " << Helper.TRF_vgetf() <<std::endl; 
	   std::cout<<"  vgetg (bar) "<< Helper.vgetg("bar")  << "  vgetf(bar)   " << Helper.vgetg("bar")<< std::endl; 
	   std::cout<<"  vgetg () "<< Helper.vgetg()  << "  vgetf()   " << Helper.vgetg()<< std::endl; 

	 }
	 if(Helper.WA()*Helper.WA()  <  Helper.WC()*Helper.WC() + Helper.WD()*Helper.WD() + Helper.WE()*Helper.WE() ){

	   std::cout<<"  WA   "<<Helper.WA()<< std::endl;
	   std::cout<<"  WC   "<<Helper.WC()<< std::endl;
          std::cout<<"  WD   "<<Helper.WD()<< std::endl;
	   std::cout<<"  WE   "<<Helper.WE()<< std::endl;

 
	 }
       cosbetacostheta_plus->Fill(Helper.cosbeta(),Helper.costhetaLF());
       TRFcosbetacostheta_plus->Fill(Helper.TRF_cosbeta(),Helper.costhetaLF());
      }
    }

    //   if(Tauola::getHelPlus() ==-1 ){
    if(HelMinus ){
      if(tau2.E()!=0 && pi2.E()!=0) { x1m = 2*pi2.E()/tau2.E()-1; pi_minus->Fill(pi2.E()/tau2.E()); piom_minus->Fill(x1m);}
      if(tau1.E()!=0 && mu1.E()!=0) {x2m = 2*mu1.E()/tau1.E()-1; mu_minus->Fill(mu1.E()/tau1.E());}

      if(tau1.E()!=0 && pi1.E()!=0){x1pim =2*pi1.E()/tau1.E()-1;}
      if(rhopi2.E()!=0 && rhopi02.E()!=0){x3m = (rhopi2.E() -   rhopi02.E())/(rhopi2.E() +   rhopi02.E());  rho_minus->Fill(x3m); }

      if(x1m!=0 && x2m!=0){OmegaMMuPi = (x1m+x2m)/(1+x1m*x2m); om_minus->Fill(OmegaMMuPi);}
      if(x2m!=0 && x3m!=0){OmegaMMuRho = (x2m+x3m)/(1+x2m*x3m); ommurho_minus->Fill(OmegaMMuRho);}


      if(x1m!=0 && x1pim!=0){ompipi_minus = (x1m+x1pim)/(1+x1m*x1pim);Omegapipi_minus->Fill(ompipi_minus);}

      if(x1pim!=0 && x3m!=0){ompirho_minus = (x1pim+x3m)/(1+x1pim*x3m);Omegapirho_minus->Fill(ompirho_minus);}
      if(JAK2==5 && SubJAK2==51){
	vector<TLorentzVector> particles;
	particles.push_back(tau1);
	particles.push_back(a1ospi);
	particles.push_back(a1ss1pi);
	particles.push_back(a1ss2pi);
	a1Helper Helper(particles, a1ospi+a1ss1pi+a1ss2pi);

 	if(Helper.getA1omega()> 0 or Helper.getA1omega() < 0)	omega_a1_minus->Fill(Helper.getA1omega());
	if(Helper.TRF_vgetA1omega()>0 or Helper.TRF_vgetA1omega()<=0)TRFomegabar_a1_minus->Fill(Helper.TRF_vgetA1omega());

	if(Helper.getg()!=0)omegabar_a1_minus->Fill(Helper.vgetA1omega("bar"));
	 // std::cout<<"  Mminus:   scalar  bar  f:     "<< Helper.vgetfscalar("bar"); 	std::cout<<"  Mminus:   scalar    f:     "<< Helper.vgetfscalar() << std::endl;
	 // std::cout<<"  Mminus:   scalar  bar g:     "<< Helper.vgetgscalar("bar");	 std::cout<<"  Mminus:   scalar   g:     "<< Helper.vgetgscalar()<<std::endl;
	std::cout<<"Minus:   --------------   " << "    costheta    " << Helper.costhetaLF()  <<"   get omega   "<< Helper.getA1omega()<<  "  omega bar   " << Helper.vgetA1omega("bar")<<  "    omega bar scalar    " << Helper.vgetA1omegascalar("bar")<<"   TRF  Omega   " <<  Helper.TRF_vgetA1omega()<<std::endl;
	 if(Helper.TRF_vgetg() > Helper.TRF_vgetf() ){
	   std::cout<<"   g > f !!!!!!!!!!!";  std::cout<<Helper.WA()*Helper.WA() << "      equaolyt    "<<   Helper.WC()*Helper.WC() + Helper.WD()*Helper.WD() + Helper.WE()*Helper.WE() <<std::endl;

	   //	   Helper.debug();
	   // std::cout<<" TRF_vgetg()  "<<  Helper.TRF_vgetg() << " TRF_vgetf()   " << Helper.TRF_vgetf() <<std::endl; 
	   // std::cout<<"  vgetg (bar) "<< Helper.vgetg("bar")  << "  vgetf(bar)   " << Helper.vgetg("bar")<< std::endl; 
	   // std::cout<<"  vgetg () "<< Helper.vgetg()  << "  vgetf()   " << Helper.vgetg()<< std::endl; 

	 }

	//	std::cout<<Helper.WA()*Helper.WA() << "   =====    "<<   Helper.WC()*Helper.WC() + Helper.WD()*Helper.WD() + Helper.WE()*Helper.WE() <<std::endl;
	//	std::cout<<"HelMinus:   costhetaLF()    "<< Helper.costhetaLF() <<"   get omega   "<< Helper.getA1omega()<<  "  omega v   " << Helper.vgetA1omega(1)<<std::endl;
	cosbetacostheta_minus->Fill(Helper.cosbeta(),Helper.costhetaLF());
	TRFcosbetacostheta_minus->Fill(Helper.TRF_cosbeta(),Helper.costhetaLF());
      }
 //  double TRF_vgetf(TString type="");
 // double TRF_vgetg(TString type="");
 // double TRF_vgetA1omega(TString type="");
 // double TRF_cosbeta();      double  TRF_cosalpha();   double  TRF_cosgamma();  
 // double TRF_sinbeta();        double TRF_sinalpha();    double  TRF_singamma();  

    
    }
    
  
    
  
    // Run MC-TESTER on the event
    HepMCEvent temp_event(*HepMCEvt,false);
    MC_Analyze(&temp_event);

    // Print some events at the end of the run
    if(iEvent>=NumberOfEvents-5){  
      // pythia.event.list();
      HepMCEvt->print();
    }

    // Clean up HepMC event
    delete HepMCEvt;  
  }
 int a, b;
 a = 5;//atoi(argv[1]);
 b =4;// atoi(argv[2]);
 MultiplyNumbers ab;
  ab.setA(a);
  ab.setB(b);
  printf(" test  %d + %d = %d\n", ab.getA(), ab.getB(), ab.getProduct());
 
  pi_plus->Write();
  pi_minus->Write();

  Omegapipi_plus->Write();
  Omegapipi_minus->Write(); 

  Omegapirho_plus->Write();
  Omegapirho_minus->Write();



  rho_plus->Write();
  rho_minus->Write();
  ommurho_minus->Write();
  ommurho_plus->Write();
  mu_plus->Write();
  mu_minus->Write();
  om_plus->Write();
  om_minus->Write();
  piom_plus->Write();
  piom_minus->Write();
  omega_a1_minus->Write();
  omega_a1_plus->Write();
  omegabar_a1_minus->Write();
omegabar_a1_plus->Write();
cosbetacostheta_minus->Write();
cosbetacostheta_plus->Write();
  file->Write();
  file->Close();
  pythia.statistics();
  MC_Finalize();

  // This is an access to old FORTRAN info on generated tau sample. 
  // That is why it refers to old version number (eg. 2.7) for TAUOLA.
  //Tauola::summary();
}


