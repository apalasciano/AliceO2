// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file preclusters-to-clusters-original-workflow.cxx
/// \brief This is an executable that runs the original MLEM cluster finder via DPL.
///
/// \author Philippe Pillot, Subatech

#include "Framework/runDataProcessing.h"

#include "MCHWorkflow/ClusterFinderOriginalSpec.h"

using namespace o2::framework;

WorkflowSpec defineDataProcessing(const ConfigContext&)
{
  return WorkflowSpec{o2::mch::getClusterFinderOriginalSpec()};
}
