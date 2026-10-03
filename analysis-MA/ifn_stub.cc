// IFN clustering wrapper.
//
// The real thing needs the IFNPlugin (fjcontrib, Caola-Grabarczyk-Hutt-
// Salam-Scyboz-Thaler).  Point IFNPLUGIN at it in the Makefile and this
// file forwards to it; with IFNPLUGIN empty the -ifn histograms simply
// stay empty and a warning is printed once.
#include <iostream>
#ifdef HAVE_IFN
#include "IFNPlugin.hh"
#include "fastjet/ClusterSequence.hh"
namespace fj = fastjet;
namespace fjc = fastjet::contrib;
#endif

extern "C" {
void fastjetifn_(const double * p, const int * flavp, const int & npart,
                 const double & R, const double & ptmin,
                 double * f77jets, int * f77isabjet, int & njets,
                 int * f77jetvec) {
#ifndef HAVE_IFN
  static bool warned = false;
  if (!warned) {
    std::cout << " ***********************************************\n"
              << "  analysis-MA: built without the IFN plugin.\n"
              << "  The  -ifn  histograms will stay empty.\n"
              << "  Set IFNPLUGIN in the Makefile to enable them.\n"
              << " ***********************************************"
              << std::endl;
    warned = true;
  }
  njets = 0;
  for (int i = 0; i < npart; i++) f77jetvec[i] = 0;
  (void)p; (void)flavp; (void)R; (void)ptmin; (void)f77jets; (void)f77isabjet;
#else
  std::vector<fj::PseudoJet> input;
  for (int i = 0; i < npart; i++) {
    fj::PseudoJet pj(p[4*i], p[4*i+1], p[4*i+2], p[4*i+3]);
    pj.set_user_info(new fjc::FlavHistory(flavp[i]));
    input.push_back(pj);
  }
  double alpha = 2.0, omega = 3.0 - alpha;
  fj::JetDefinition base(fj::antikt_algorithm, R);
  fjc::FlavRecombiner flav_recombiner;
  base.set_recombiner(&flav_recombiner);
  fjc::IFNPlugin * plugin =
      new fjc::IFNPlugin(base, alpha, omega, fjc::FlavRecombiner::net);
  fj::JetDefinition jd(plugin);
  fj::ClusterSequence cs(input, jd);
  std::vector<fj::PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(ptmin));
  njets = jets.size();
  for (int i = 0; i < njets; i++) {
    f77jets[4*i+0] = jets[i].px(); f77jets[4*i+1] = jets[i].py();
    f77jets[4*i+2] = jets[i].pz(); f77jets[4*i+3] = jets[i].E();
    auto fl = fjc::FlavHistory::current_flavour_of(jets[i]);
    fl.reset_all_but_flav(5);
    fl.apply_modulo_2();
    f77isabjet[i] = fl.is_flavourless() ? 0 : 1;
  }
  std::vector<int> cl = cs.particle_jet_indices(jets);
  for (int i = 0; i < npart; i++) f77jetvec[i] = cl[i] + 1;
  delete plugin;
#endif
}
}
