// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file HFXToJpsiPiPiCandidateSelector.cxx
/// \brief X(3872) selection task.
/// \note Adapted from HFJpsiToEECandidateSelector.cxx
/// \author Rik Spijkers <r.spijkers@students.uu.nl>, Utrecht University

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "AnalysisCore/HFSelectorCuts.h"
#include "AnalysisDataModel/HFSecondaryVertex.h"
#include "AnalysisDataModel/HFCandidateSelectionTables.h"
using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::aod::hf_cand_bplus;
using namespace o2::analysis;
using namespace o2::analysis::hf_cuts_bplus_tod0pi;

/// Struct for applying B+ selection cuts
struct HFBPlusToD0PiCandidateSelector {

  Produces<aod::HFSelBPlusToD0PiCandidate> hfSelBPlusToD0PiCandidate;

  Configurable<double> d_pTCandMin{"d_pTCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> d_pTCandMax{"d_pTCandMax", 50., "Upper bound of candidate pT"};

  Configurable<double> d_pidTPCMinpT{"d_pidTPCMinpT", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<double> d_pidTPCMaxpT{"d_pidTPCMaxpT", 10., "Upper bound of track pT for TPC PID"};
  Configurable<double> d_pidTOFMinpT{"d_pidTOFMinpT", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<double> d_pidTOFMaxpT{"d_pidTOFMaxpT", 10., "Upper bound of track pT for TOF PID"};

  Configurable<double> d_TPCNClsFindablePIDCut{"d_TPCNClsFindablePIDCut", 70., "Lower bound of TPC findable clusters for good PID"};
  Configurable<double> d_nSigmaTPC{"d_nSigmaTPC", 3., "Nsigma cut on TPC only"};
  Configurable<double> d_nSigmaTOF{"d_nSigmaTOF", 3., "Nsigma cut on TOF only"};

  Configurable<std::vector<double>> pTBins{"pTBins", std::vector<double>{hf_cuts_bplus_tod0pi::pTBins_v}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"BPlus_to_d0pi_cuts", {hf_cuts_bplus_tod0pi::cuts[0], npTBins, nCutVars, pTBinLabels, cutVarLabels}, "B+ candidate selection per pT bin"};

  /// Selection on goodness of daughter tracks
  /// \note should be applied at candidate selection
  /// \param track is daughter track
  /// \return true if track is good
  template <typename T>
  bool daughterSelection(const T& track)
  {
    return true;
  }

  /// Conjugate independent toplogical cuts
  /// \param hfCandBplus is candidate
  /// \param trackNeutral is the track compatible with the D0(bar) hypothesis
  /// \param trackPos is the track with the pi+ hypothesis
  /// \return true if candidate passes all cuts
  template <typename T1, typename T2, typename T3>
  bool selectionTopol(const T1& hfCandBPlus, const T2& hfCandD0, const T3& trackPos)
  {
    auto candpT = hfCandBPlus.pt();
    int pTBin = findBin(pTBins, candpT);
    if (pTBin == -1) {
      // Printf("B+ topol selection failed at getpTBin");
      return false;
    }

    if (candpT < d_pTCandMin || candpT >= d_pTCandMax) {
      // Printf("B+ topol selection failed at cand pT check");
      return false; //check that the candidate pT is within the analysis range
    }

    if (TMath::Abs(InvMassBplus(hfCandBPlus) - RecoDecay::getMassPDG(521)) > cuts->get(pTBin, "m")) {
      // Printf("B+ topol selection failed at mass diff check");
      return false; //check that mass difference is within bounds
    }

    if ((hfCandD0.pt() < cuts->get(pTBin, "pT D0")) || (trackPos.pt() < cuts->get(pTBin, "pT Pi")) ) {
      // Printf("B+ topol selection failed at daughter pT check");
      return false; //cut on daughter pT
    }

    if (TMath::Abs(hfCandBPlus.cpa()) < cuts->get(pTBin, "CPA")) {
      return false; // CPA check
    }

    if ((TMath::Abs(hfCandBPlus.impactParameter0()) > cuts->get(pTBin, "d0 D0")) ||
        (TMath::Abs(hfCandBPlus.impactParameter1()) > cuts->get(pTBin, "d0 Pi"))) { 
      return false; // DCA check on daughters
    } 

    if (TMath::Abs(hfCandBPlus.decayLength()) < cuts->get(pTBin, "B decLen")) {
      return false; 
    }//"B decLen"

    if (TMath::Abs(InvMassD0(hfCandD0) - RecoDecay::getMassPDG(421)) > cuts->get(pTBin, "DeltaMD0")) {
      // Printf("D0 topol selection failed at mass diff check");
      return false; //check that mass difference is within bounds
    }

    if (TMath::Abs(hfCandD0.cpa()) < cuts->get(pTBin, "CPA D0")) {
      return false; // CPA check
    }


    // More cuts have to be addee here: do of both the decay products, some information aout the D0( pointing angle for sure, mass, theta*, DecayLength?)

    return true;
  }

  /// Check if track is ok for TPC PID
  /// \param track is the track
  /// \note function to be expanded
  /// \return true if track is ok for TPC PID
  template <typename T>
  bool validTPCPID(const T& track)
  {
    if (TMath::Abs(track.pt()) < d_pidTPCMinpT || TMath::Abs(track.pt()) >= d_pidTPCMaxpT) {
      return false;
    }
    //if (track.TPCNClsFindable() < d_TPCNClsFindablePIDCut) return false;
    return true;
  }

  /// Check if track is ok for TOF PID
  /// \param track is the track
  /// \note function to be expanded
  /// \return true if track is ok for TOF PID
  template <typename T>
  bool validTOFPID(const T& track)
  {
    if (TMath::Abs(track.pt()) < d_pidTOFMinpT || TMath::Abs(track.pt()) >= d_pidTOFMaxpT) {
      return false;
    }
    return true;
  }

  /// Check if track is compatible with given TPC Nsigma cut for the pion hypothesis
  /// \param track is the track
  /// \param nSigmaCut is the nsigma threshold to test against
  /// \return true if track satisfies TPC pion hypothesis for given Nsigma cut
  template <typename T>
  bool selectionPIDTPC(const T& track, int nSigmaCut)
  {
    if (nSigmaCut > 999.) {
      return true;
    }
    return track.tpcNSigmaPi() < nSigmaCut;
  }

  /// Check if track is compatible with given TOF NSigma cut for the pion hypothesis
  /// \param track is the track
  /// \param nSigmaCut is the nSigma threshold to test against
  /// \return true if track satisfies TOF pion hypothesis for given NSigma cut
  template <typename T>
  bool selectionPIDTOF(const T& track, double nSigmaCut)
  {
    if (nSigmaCut > 999.) {
      return true;
    }
    return track.tofNSigmaPi() < nSigmaCut;
  }

  /// PID selection on daughter track
  /// \param track is the daughter track
  /// \return 1 if successful PID match, 0 if successful PID rejection, -1 if no PID info
  template <typename T>
  int selectionPID(const T& track)
  { // use both TPC and TOF here; in run5 only TOF makes sense. add some flag for run3/run5 data later?
    if (validTOFPID(track)) {
      if (!selectionPIDTOF(track, d_nSigmaTOF)) {
        return 0; //rejected by PID
      } else {
        return 1; //positive PID
      }
    } else {
      return -1; //no PID info
    }

    // no tpc in run5, so for now comment it out
    // if (validTPCPID(track)) {
    //   if (!selectionPIDTPC(track, d_nSigmaTPC)) {
    //     return 0; //rejected by PID
    //   } else {
    //     return 1; //positive PID
    //   }
    // } else {
    //   return -1; //no PID info
    // }
  }
  void process(aod::HfCandBPlus const& hfCandBs, aod::HfCandProng2, aod::BigTracksPID const& tracks)
  {
    for (auto& hfCandB : hfCandBs) { //looping over Bplus candidates
      // D0 is always index0 by default
      auto candD0 = hfCandB.index0(); // not a track, is by default a D0 particle
      auto trackPos = hfCandB.index1_as<aod::BigTracksPID>(); //positive daughter

      int statusBplus= 0;

      // check if flagged as B+ --> D0bar Pi
      if (!(hfCandB.hfflag() & 1 << BPlusToD0Pi)) {
        hfSelBPlusToD0PiCandidate(statusBplus);
        // Printf("B+ candidate selection failed at hfflag check");
        continue;
      }

      // daughter track validity selection
      if (!daughterSelection(trackPos)) {
        hfSelBPlusToD0PiCandidate(statusBplus);
        // Printf("B+ candidate selection failed at daughter selection");
        continue;
      }

      if (!selectionTopol(hfCandB, candD0, trackPos)) {
        hfSelBPlusToD0PiCandidate(statusBplus);
        // Printf("B+ candidate selection failed at selection topology");
        continue;
      }

      if (selectionPID(trackPos) == 0) {
        hfSelBPlusToD0PiCandidate(statusBplus);
        // Printf("B+ candidate selection failed at selection PID");
        continue;
      }

      hfSelBPlusToD0PiCandidate(1);
      // Printf("B+ candidate selection successful, candidate should be selected");
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HFBPlusToD0PiCandidateSelector>(cfgc, TaskName{"hf-bplus-tod0pi-candidate-selector"})};
}