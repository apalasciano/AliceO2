// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file taskBPlus.cxx
/// \brief B+ analysis task
///
/// \author Gian Michele Innocenti <gian.michele.innocenti@cern.ch>, CERN
/// \author Antonio Palasciano <antonio.palasciano@cern.ch>,
/// \author Deepa Thomas <deepa.thomas@cern.ch>, UT Austin

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "AnalysisDataModel/HFSecondaryVertex.h"
#include "AnalysisCore/HFSelectorCuts.h"
#include "AnalysisDataModel/HFCandidateSelectionTables.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::aod::hf_cand_prong2;
using namespace o2::aod::hf_cand_bplus;
using namespace o2::analysis::hf_cuts_bplus_tod0pi;
using namespace o2::framework::expressions;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  ConfigParamSpec optionDoMC{"doMC", VariantType::Bool, false, {"Fill MC histograms."}};
  workflowOptions.push_back(optionDoMC);
}

#include "Framework/runDataProcessing.h"

/// BPlus analysis task
struct TaskBPlus {
  HistogramRegistry registry{
    "registry",
    {{"hPtProng0", "2-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
     {"hPtProng1", "2-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
     {"hPtCand", "2-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}}}};

  // Configurable<int> selectionFlagBPlus{"selectionFlagBPlus", 1, "Selection Flag for BPlus"};
  Configurable<int> d_selectionFlagBPlus{"d_selectionFlagBPlus", 1, "Selection Flag for B+"};
  Configurable<double> cutYCandMax{"cutYCandMax", -1., "max. cand. rapidity"};
  Configurable<std::vector<double>> bins{"pTBins", std::vector<double>{hf_cuts_bplus_tod0pi::pTBins_v}, "pT bin limits"};

  void init(o2::framework::InitContext&)
  {
    registry.add("hMass", "2-prong candidates;inv. mass (D0  #pi+ #pi-) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{500, 0., 10.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hdeclength", "2-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{200, 0., 0.4}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hdeclengthxy", "2-prong candidates;decay length xy (cm);entries", {HistType::kTH2F, {{200, 0., 0.4}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0", "2-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -0.05, 0.05}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1", "2-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -0.05, 0.05}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPA", "2-prong candidates;cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEta", "2-prong candidates;candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpParErr", "2-prong candidates;impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLenErr", "2-prong candidates;decay length error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLenXYErr", "2-prong candidates;decay length xy error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
  }

  Filter filterSelectCandidates = (aod::hf_selcandidate_bplus::isSelBPlusToD0Pi >= d_selectionFlagBPlus);

  void process(soa::Filtered<soa::Join<aod::HfCandBPlus, aod::HFSelBPlusToD0PiCandidate>> const& candidates)
  {
    for (auto& candidate : candidates) {
      if (!(candidate.hfflag() & 1 << BPlusToD0Pi)) {
        continue;
      }
      if (cutYCandMax >= 0. && std::abs(YBplus(candidate)) > cutYCandMax) {
        continue;
      }

      registry.fill(HIST("hMass"), InvMassBplus(candidate), candidate.pt());
      registry.fill(HIST("hPtCand"), candidate.pt());
      registry.fill(HIST("hPtProng0"), candidate.ptProng0());
      registry.fill(HIST("hPtProng1"), candidate.ptProng1());
      registry.fill(HIST("hdeclength"), candidate.decayLength(), candidate.pt());
      registry.fill(HIST("hdeclengthxy"), candidate.decayLengthXY(), candidate.pt());
      registry.fill(HIST("hd0Prong0"), candidate.impactParameter0(), candidate.pt());
      registry.fill(HIST("hd0Prong1"), candidate.impactParameter1(), candidate.pt());
      registry.fill(HIST("hCPA"), candidate.cpa(), candidate.pt());
      registry.fill(HIST("hEta"), candidate.eta(), candidate.pt());
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter0(), candidate.pt());
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter1(), candidate.pt());
      registry.fill(HIST("hDecLenErr"), candidate.errorDecayLength(), candidate.pt());
      registry.fill(HIST("hDecLenXYErr"), candidate.errorDecayLengthXY(), candidate.pt());
    } // candidate loop
  }   // process
};    // struct

struct TaskBPlusMC {
  HistogramRegistry registry{
    "registry",
    {{"hPtRecSig", "2-prong candidates (rec. matched);#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
     {"hPtRecBg", "2-prong candidates (rec. unmatched);#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
     {"hPtGen", "2-prong candidates (gen. matched);#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
     {"hPtGenSig", "2-prong candidates (rec. matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
     {"hMassRecSigD0", "2-prong candidates (rec. matched);inv. mass (D0 #pi+) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 1.5, 2.}}}},
     {"hMassRecBgD0", "2-prong candidates (rec. matched);inv. mass (D0 #pi+) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 1.5, 2.}}}}}};


  Configurable<int> d_selectionFlagBPlus{"d_selectionFlagBPlus", 1, "Selection Flag for B+"};
  Configurable<double> cutYCandMax{"cutYCandMax", -1., "max. cand. rapidity"};
  Configurable<std::vector<double>> bins{"pTBins", std::vector<double>{hf_cuts_bplus_tod0pi::pTBins_v}, "pT bin limits"};

    void init(o2::framework::InitContext&)
  {
    registry.add("hCPARecSig", "2-prong candidates (rec. matched);cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPARecBg", "2-prong candidates (rec. unmatched);cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEtaRecSig", "2-prong candidates (rec. matched);#it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEtaRecBg", "2-prong candidates (rec. unmatched);#it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEtaGen", "2-prong candidates (gen. matched);#it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});

    registry.add("hPtProng0RecSig", "2-prong candidates (rec. matched);prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng1RecSig", "2-prong candidates (rec. matched);prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng0RecBg", "2-prong candidates (rec. unmatched);prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng1RecBg", "2-prong candidates (rec. unmatched);prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtGenProng0", "2-prong candidates (gen. matched);prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtGenProng1", "2-prong candidates (gen. matched);prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});

    registry.add("hMassRecSig", "2-prong candidates (rec. matched);inv. mass (D0 #pi+) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{500, 5.0, 5.40}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassRecBg", "2-prong candidates (rec. unmatched);inv. mass (D0  #pi+) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{500, 0., 10.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0RecSig", "2-prong candidates (rec. matched);prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.02, 0.02}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1RecSig", "2-prong candidates (rec. matched);prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.02, 0.02}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0RecBg", "2-prong candidates (rec. unmatched);prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.02, 0.02}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1RecBg", "2-prong candidates (rec. unmatched);prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.02, 0.02}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDeclengthRecSig", "2-prong candidates (rec. matched);decay length (cm);entries", {HistType::kTH2F, {{200, 0., 0.04}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDeclengthRecBg", "2-prong candidates (rec. unmatched);decay length (cm);entries", {HistType::kTH2F, {{200, 0., 0.04}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
  }

  Filter filterSelectCandidates = (aod::hf_selcandidate_bplus::isSelBPlusToD0Pi >= d_selectionFlagBPlus);

  void process(soa::Filtered<soa::Join<aod::HfCandBPlus, aod::HFSelBPlusToD0PiCandidate, aod::HfCandBPMCRec>> const& candidates, 
               soa::Join<aod::McParticles, aod::HfCandBPMCGen> const& particlesMC, aod::BigTracksMC const& tracks)
  {
    // MC rec.
    for (auto& candidate : candidates) {
      if (!(candidate.hfflag() & 1 << BPlusToD0Pi)) {
        continue;
      }
      if (cutYCandMax >= 0. && std::abs(YBplus(candidate)) > cutYCandMax) {
        continue;
      }
      if (std::abs(candidate.flagMCMatchRec()) == 1 << BPlusToD0Pi) {
        auto indexMother = RecoDecay::getMother(particlesMC, candidate.index1_as<aod::BigTracksMC>().mcParticle_as<soa::Join<aod::McParticles, aod::HfCandBPMCGen>>(), 521, true);
        auto particleMother = particlesMC.iteratorAt(indexMother);
        registry.fill(HIST("hPtGenSig"), particleMother.pt());
        registry.fill(HIST("hPtRecSig"), candidate.pt());
        registry.fill(HIST("hCPARecSig"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("hEtaRecSig"), candidate.eta(), candidate.pt());
        //registry.fill(HIST("hMassRecSigD0"), InvMassD0(prong0D0));
        registry.fill(HIST("hDeclengthRecSig"), candidate.decayLength(), candidate.pt());
        registry.fill(HIST("hMassRecSig"), InvMassBplus(candidate), candidate.pt());
        registry.fill(HIST("hd0Prong0RecSig"), candidate.impactParameter0(), candidate.pt());
        registry.fill(HIST("hd0Prong1RecSig"), candidate.impactParameter1(), candidate.pt());
        registry.fill(HIST("hPtProng0RecSig"), candidate.ptProng0(), candidate.pt());
        registry.fill(HIST("hPtProng1RecSig"), candidate.ptProng1(), candidate.pt());
      } else {
        registry.fill(HIST("hPtRecBg"), candidate.pt());
        registry.fill(HIST("hCPARecBg"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("hEtaRecBg"), candidate.eta(), candidate.pt());
        registry.fill(HIST("hDeclengthRecBg"), candidate.decayLength(), candidate.pt());
        registry.fill(HIST("hMassRecBg"), InvMassBplus(candidate), candidate.pt());
        registry.fill(HIST("hd0Prong0RecBg"), candidate.impactParameter0(), candidate.pt());
        registry.fill(HIST("hd0Prong1RecBg"), candidate.impactParameter1(), candidate.pt());
        registry.fill(HIST("hPtProng0RecBg"), candidate.ptProng0(), candidate.pt());
        registry.fill(HIST("hPtProng1RecBg"), candidate.ptProng1(), candidate.pt());
      }
    } // rec
    // MC gen. level
    //Printf("MC Particles: %d", particlesMC.size());
    for (auto& particle : particlesMC) {
      if (std::abs(particle.flagMCMatchGen()) == 1  << BPlusToD0Pi) {
        if (cutYCandMax >= 0. && std::abs(RecoDecay::Y(array{particle.px(), particle.py(), particle.pz()}, RecoDecay::getMassPDG(521))) > cutYCandMax) {
        continue;
        }
        registry.fill(HIST("hPtGen"), particle.pt());
        registry.fill(HIST("hEtaGen"), particle.eta(), particle.pt());

        float ptProngs[2];
        int counter = 0;
        for (int iD = particle.daughter0(); iD <= particle.daughter1(); ++iD) {
          ptProngs[counter] = particlesMC.iteratorAt(iD).pt();
          counter++;
        }
        registry.fill(HIST("hPtGenProng0"), ptProngs[0], particle.pt());
        registry.fill(HIST("hPtGenProng1"), ptProngs[1], particle.pt());
      }
    } //gen
  }   // process
};    // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<TaskBPlus>(cfgc, TaskName{"hf-task-bplus"})};
  const bool doMC = cfgc.options().get<bool>("doMC");
  if (doMC) {
    workflow.push_back(adaptAnalysisTask<TaskBPlusMC>(cfgc, TaskName{"hf-task-bplus-mc"}));
  }
  return workflow;
}