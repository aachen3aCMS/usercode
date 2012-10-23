#!/bin/bash
######################################################################
# skim.sh - Skim many datasets at once
# (C) 2011 - 2012 Martin Weber
######################################################################

# single required parameter: the dataset, and any optional parameters
submit_mc()
{
    # Remove starting slash and ending, replace slashes by pluses
    TAG=`echo $1 | sed -e s:/AOD::g -e s:^/::g -e s:/:_:g `
    DS=$1
    shift
    ./prepare.sh -f -m -s $DS -t ${TAG} -g "${GLOBALTAG}::All" -c SkimmingCuts.cfg ${*}
}

# two obligatory parameters: the dataset and the json file, and any optional parameters
submit_data()
{
    # Remove starting slash and ending, replace slashes by pluses
    TAG=`echo $1 | sed -e s:/AOD::g -e s:^/::g -e s:/:+:g `
    DS=$1
    JSON=$2
    shift
    shift
    ./prepare.sh -f -d -s $DS -t ${TAG} -g "${GLOBALTAG}::All" -j $JSON -c SkimmingCuts.cfg $*
}

######################################################################
# Data 2012

# json files can be found in /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV
# just take DCSonly JSON file and filter events later in analysis (RunLumiRanges.cpp)
JSON="20121032_json_DCSONLY.txt"

# 2012 A+B ReReco, CMSSW_5_3_2_patch4
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions#2012_A_and_B_data_re_reco_with_C
GLOBALTAG="FT_53_V6_AN2"
submit_data /DoubleMu/Run2012A-13Jul2012-v1/AOD $JSON
#submit_data /DoubleMu/Run2012B-13Jul2012-v4/AOD $JSON

# 2012 special ECAL recovery (A+B), CMSSW_5_3_3_patch1
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions#2012_A_and_B_data_re_reco_with_C
GLOBALTAG="FT_53_V6C_AN1"
#submit_data /DoubleMu/Run2012A-recover-06Aug2012-v1/AOD $JSON

# 2012 special PFPhoton-PFEelectron bug fix rereco (ReReco24Aug)
GLOBALTAG="FT_53_V10_AN1"

# 2012 C PromptReco analysis, CMSSW_5_3_?????
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions#Global_tags_for_prompt_reco_data
GLOBALTAG="GR_P_V40_AN1"
#submit_data /DoubleMu/Run2012C-PromptReco-v1/AOD $JSON
#submit_data /DoubleMu/Run2012C-PromptReco-v2/AOD $JSON

######################################################################
# MC Summer 12 - July ReReco

# CMSSW_5_3_2_patch4
GLOBALTAG="START53_V7F"

#submit_mc /DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM

#submit_mc /soft332_m0eq1000_m12eq500_tanbeq20_sgnMueq1_A0eq0_lp211eq0.01_hw4000_CMSSW523p3_noPU_AOD_20.08.2012/sonnen-soft332_m0eq1000_m12eq500_tanbeq20_sgnMueq1_A0eq0_lp211eq0.01_hw4000_CMSSW523p3_noPU_AOD_20.08.2012-1c2ba989717b25223297e99c10c41dbd/USER
