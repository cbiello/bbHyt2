//Fastjet CCP-Fortran interface for IFN algorithm
//It contains the routines for the antikt flavour-jet clustering algorithm
//of Ref.[2306.07314]. We provide a POWHEG interface by including the old
//routine for the standard flavour-blind algorithm.
//
//Author: Christian Biello
//Date: 19 Sep 24

#include <iostream>
#include <iomanip>
#include <cstdlib>

#include "/u/asankar/Pkgs/fastjet-install/include/fastjet/ClusterSequence.hh"
#include "/u/asankar/Pkgs/fastjet-install/include/fastjet/SISConePlugin.hh"
#include "/u/asankar/POWHEG-BOX-RES/MiNNLOPS_res/bbH4FS_local/IFNPlugin/IFNPlugin.hh"

namespace fj = fastjet;
namespace fjc = fastjet::contrib;
using namespace std;


extern "C" {   
// f77 interface to the pp generalised-kt (sequential recombination)
// algorithms, as defined in arXiv.org:0802.1189, which includes
// kt, Cambridge/Aachen and anti-kt as special cases.
void fastjetppgenkt_(const double * p, const int & npart,                   
                     const double & R, const double & palg, const double & ptmin,
                     double * f77jets, int & njets, int * f77jetvec) {

    // transfer p[4*ipart+0..3] -> input_particles[i]
    vector<fj::PseudoJet> input_particles;   
    for (int i=0; i<npart; i++) {
      valarray<double> mom(4); // mom[0..3]
      for (int j=0;j<=3; j++) {
         mom[j] = *(p++);
      }
      fj::PseudoJet psjet(mom);
      input_particles.push_back(psjet);    
      // label input_particles entries
      input_particles[i].set_user_index(i+1);
    }
    
    // prepare jet def and run fastjet
    fj::JetDefinition jet_def;
    if (palg == 1.0) {
      jet_def = fj::JetDefinition(fj::kt_algorithm, R);
    }  else if (palg == 0.0) {
      jet_def = fj::JetDefinition(fj::cambridge_algorithm, R);
    }  else if (palg == -1.0) {
      jet_def = fj::JetDefinition(fj::antikt_algorithm, R);
    } else {
      jet_def = fj::JetDefinition(fj::genkt_algorithm, R, palg);
    }

    
    // perform clustering
    fj::ClusterSequence cs(input_particles, jet_def);
    // extract jets (pt-ordered)
    vector<fj::PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(ptmin));
    njets = jets.size();

    // find particles inside i-th jet
    vector<fj::PseudoJet> *constit;
    constit=new vector<fj::PseudoJet>[njets];
    for(int i=0; i<njets; i++) {
      constit[i] = cs.constituents(jets[i]); 
      //cout<<"jet "<<i<<endl;
      //cout<<"mult "<<constit[i].size()<<endl;
      for(int j=0; j<constit[i].size(); j++) {
	*(f77jetvec + constit[i][j].user_index()-1) = i+1;
      }
    }

    // transfer jets -> f77jets[4*ijet+0..3]
    for (int i=0; i<njets; i++) {
      for (int j=0;j<=3; j++) {
        *f77jets = jets[i][j];
        f77jets++;
      } 
    }

    // clean up
    delete [] constit;
    
   }

//IFN routine
//You need to pass the flavour of the particles via * flavp,
//and you will obtain the flavour of the jets via * f77isabjet:
//isabjet=0 is a light-jet
//isabjet=1 is a b-jet
//according to the modulo2 definition.
  
  void fastjetifnclustering_(const double * p, const int * flavp, const int & npart,                   
                     const double & R, const double & ptmin,
			     double * f77jets, int * f77isabjet, int & njets, int * f77jetvec){

    // transfer p[4*ipart+0..3] -> input_particles[i]
    vector<fj::PseudoJet> input_particles;   
    for (int i=0; i<npart; i++) {
      valarray<double> mom(4); // mom[0..3]
      int pdgflav;
      for (int j=0;j<=3; j++) {
         mom[j] = *(p++);
      }
      pdgflav = *(flavp++);
      fj::PseudoJet psjet(mom);
      psjet.set_user_info(new fjc::FlavHistory(pdgflav));
      input_particles.push_back(psjet);
      // label input_particles entries
      input_particles[i].set_user_index(i+1);
    }
    
    //prepare a basic jet definition
    fj::JetDefinition base_jet_def(fj::antikt_algorithm, R);
    // enable it to track flavours (default is net flavour)
    fjc::FlavRecombiner flav_recombiner;
    base_jet_def.set_recombiner(&flav_recombiner);

    // And then we set up the IFNPlugin that builds on the base_jet_def
    // The main free parameter, alpha, in the uij distance, 
    //   uij = max(pt_i, pt_j)^alpha min(pt_i, pt_j)^(2-alpha) Omega_ij
    // See eq.(2) of Ref.[2306.07314]
    double alpha = 2.0;
    
    // The flavour summation scheme; should be one of 
    //   - FlavRecombiner::net
    //   - FlavRecombiner::modulo_2
    fjc::FlavRecombiner::FlavSummation flav_summation = fjc::FlavRecombiner::modulo_2;

    //define the jet clustering algo by using the IFNPlugin
    fj::JetDefinition::Plugin * pluginifn = new fjc::IFNPlugin( base_jet_def, alpha );
    fj::JetDefinition jet_def_ifn(pluginifn);
    
    vector<fj::PseudoJet> jets = jet_def_ifn(input_particles);
    vector<fj::PseudoJet> basejets = base_jet_def(input_particles);

    njets = jets.size();

    //create the flav vector
    std::vector<int> is_a_bjet(njets);
    
    // ----------------------------------------------------
    // some prints for understanding the flow and
    // a routine for saving modulo2 flav
    for (unsigned int ijet = 0; ijet < jets.size() && ijet < njets; ijet++) {

      // first print out the original anti-kt jets and then IFN jets
      //anti-kt jet
      //      const auto & base_jet  = basejets [ijet];
      //      cout << "base jet  " << ijet << ": ";
      //      cout << "pt=" << base_jet.pt() << " rap=" << base_jet.rap() << " phi=" << base_jet.phi();
      //      cout << ", flav = " << fjc::FlavHistory::current_flavour_of(base_jet).description() << endl;
      
      const auto & IFN_jet  = jets [ijet];
      //            cout << "IFN jet  " << ijet << ": ";
      //            cout << "pt=" << IFN_jet.pt() << " rap=" << IFN_jet.rap() << " phi=" << IFN_jet.phi();
      //            cout << ", flav = " << fjc::FlavHistory::current_flavour_of(IFN_jet).description() << endl;

      auto flavmodulo2=fjc::FlavHistory::current_flavour_of(IFN_jet);
      flavmodulo2.reset_all_but_flav(5);
      flavmodulo2.apply_modulo_2();
      //            cout << ", modulo2 flav = " << flavmodulo2.description() << endl;
      
      if (flavmodulo2.is_multiflavoured()){
	cout << "CB: error in IFN modulo2 should not have a multiflavoured jet" << endl;
	exit(1);
      } else if (flavmodulo2.is_flavourless()) {
	is_a_bjet[ijet]=0;
      } else {
	is_a_bjet[ijet]=1;
      }

      //      cout << "is_a_bjet= " << is_a_bjet[ijet];
      //      cout << endl;
      
      }

    // transfer jets -> f77jets[4*ijet+0..3]
    for (int i=0; i<njets; i++) {
      for (int j=0;j<=3; j++) {
        *f77jets = jets[i][j];
        f77jets++;
      }
    }

    // transfer is_a_bjet -> f77isabjet[ijet]
    for (int i=0; i<njets; i++) {
      *f77isabjet = is_a_bjet[i];
      f77isabjet++;
    }
    
    }

}
