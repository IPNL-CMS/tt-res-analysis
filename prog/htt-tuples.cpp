/**
 * This program produces ROOT trees with input variables for the H->tt analysis. Systematic
 * variations are supported.
 */

#include <BasicObservables.hpp>
#include <DumpWeights.hpp>
#include <LOSystWeights.hpp>
#include <TopPtWeight.hpp>
#include <TTObservables.hpp>
#include <TTSemilepRecoRochester.hpp>

#include <mensura/core/BTagWPService.hpp>
#include <mensura/core/Dataset.hpp>
#include <mensura/core/FileInPath.hpp>
#include <mensura/core/RunManager.hpp>
#include <mensura/core/SystService.hpp>

#include <mensura/extensions/BTagEffService.hpp>
#include <mensura/extensions/BTagSFService.hpp>
#include <mensura/extensions/BTagWeight.hpp>
#include <mensura/extensions/DatasetBuilder.hpp>
#include <mensura/extensions/GenWeightSyst.hpp>
#include <mensura/extensions/JetFilter.hpp>
#include <mensura/extensions/LeptonFilter.hpp>
#include <mensura/extensions/LeptonSFWeight.hpp>
#include <mensura/extensions/MetFilter.hpp>
#include <mensura/extensions/PileUpWeight.hpp>
#include <mensura/extensions/TFileService.hpp>
#include <mensura/extensions/WeightCollector.hpp>

#include <mensura/PECReader/PECGeneratorReader.hpp>
 #include <mensura/PECReader/PECGenParticleReader.hpp>
#include <mensura/PECReader/PECInputData.hpp>
#include <mensura/PECReader/PECJetMETReader.hpp>
#include <mensura/PECReader/PECLeptonReader.hpp>
#include <mensura/PECReader/PECPileUpReader.hpp>
#include <mensura/PECReader/PECTriggerFilter.hpp>

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

#include <iostream>
#include <list>
#include <regex>
#include <sstream>


using namespace std;
namespace po = boost::program_options;


enum class Channel
{
    Electron,
    Muon
};

enum class SampleGroup
{
    Data,
    TT,
    OtherMC
};


int main(int argc, char **argv)
{
    // Parse arguments
    po::options_description options("Allowed options");
    options.add_options()
      ("help,h", "Prints help message")
      ("channel", po::value<string>(), "Lepton channel (required argument)")
      ("samples", po::value<string>(), "Group of input samples (required argument)")
      ("syst,s", po::value<string>(), "Systematic shift");
    
    po::positional_options_description positionalOptions;
    positionalOptions.add("channel", 1);
    positionalOptions.add("samples", 1);
    
    po::variables_map optionsMap;
    po::store(
      po::command_line_parser(argc, argv).options(options).positional(positionalOptions).run(),
      optionsMap);
    po::notify(optionsMap);
    
    if (optionsMap.count("help"))
    {
        cerr << "Produces ROOT trees with input variables for the H->tt analysis.\n";
        cerr << "Usage: htt-tuples channel samples [options]\n";
        cerr << options << endl;
        return EXIT_FAILURE;
    }
    
    
    if (not optionsMap.count("channel"))
    {
        cerr << "Required argument \e[1mchannel\e[0m is missing.\n";
        return EXIT_FAILURE;
    }
    
    string channelText(optionsMap["channel"].as<string>());
    Channel channel;
    
    if (channelText == "mu")
        channel = Channel::Muon;
    else if (channelText == "e")
        channel = Channel::Electron;
    else
    {
        cerr << "Cannot recognize channel \"" << channelText << "\".\n";
        return EXIT_FAILURE;
    }
    
    
    if (not optionsMap.count("samples"))
    {
        cerr << "Required argument \e[1msamples\e[0m is missing.\n";
        return EXIT_FAILURE;
    }
    
    string sampleGroupText(optionsMap["samples"].as<string>());
    SampleGroup sampleGroup;
    
    if (sampleGroupText == "data")
        sampleGroup = SampleGroup::Data;
    else if (sampleGroupText == "tt")
        sampleGroup = SampleGroup::TT;
    else if (sampleGroupText == "other")
        sampleGroup = SampleGroup::OtherMC;
    else
    {
        cerr << "Cannot recognize sample group \"" << sampleGroupText << "\". Supported groups "
          "are \"data\", \"tt\", and \"other\" (for other simulation)\n";
        return EXIT_FAILURE;
    }
    
    
    string systType("None");
    SystService::VarDirection systDirection = SystService::VarDirection::Undefined;
    
    if (optionsMap.count("syst"))
    {
        if (sampleGroup == SampleGroup::Data)
        {
            cerr << "Cannot perform systematic variations in collision data.\n";
            return EXIT_FAILURE;
        }
        
        
        string systArg(optionsMap["syst"].as<string>());
        boost::to_lower(systArg);
        
        std::regex systRegex("(jec|jer|metuncl)[-_]?(up|down)", std::regex::extended);
        std::smatch matchResult;
        
        if (not std::regex_match(systArg, matchResult, systRegex))
        {
            cerr << "Cannot recognize systematic variation \"" << systArg << "\".\n";
            return EXIT_FAILURE;
        }
        else
        {
            if (matchResult[1] == "jec")
                systType = "JEC";
            else if (matchResult[1] == "jer")
                systType = "JER";
            else if (matchResult[1] == "metuncl")
                systType = "METUncl";
            
            if (matchResult[2] == "up")
                systDirection = SystService::VarDirection::Up;
            else if (matchResult[2] == "down")
                systDirection = SystService::VarDirection::Down;
        }
    }
    
    
    // Add a new search path
    string const installPath(getenv("TTRES_ANALYSIS_INSTALL"));
    FileInPath::AddLocation(installPath + "/data/");
    FileInPath::AddLocation(installPath + "/config/");
    
    
    // Input datasets
    list<Dataset> datasets;
    DatasetBuilder datasetBuilder("/gridgroup/cms/popov/PECData/2016Charlie/samples_v3.json");
    
    if (sampleGroup == SampleGroup::Data)
    {
        if (channel == Channel::Muon)
            datasets = datasetBuilder({"SingleMuon-Run2016_330_all"});
        else
            datasets = datasetBuilder({"SingleElectron-Run2016_330_all"});
    }
    else if (sampleGroup == SampleGroup::TT)
    {
        datasets = datasetBuilder({"ttbar-pw_330_all"});
        
        if (systType == "None")
            datasets.splice(datasets.end(), datasetBuilder({
              "ttbar-pw-isrup_330_ojF", "ttbar-pw-isrdown_330_all",
              "ttbar-pw-fsrup_330_all", "ttbar-pw-fsrdown_330_all",
              "ttbar-pw-hdampup_330_rAN", "ttbar-pw-hdampdown_330_HTS",
              "ttbar-pw-m1755_330_all", "ttbar-pw-m1695_330_all",
              "ttbar-pw-ueup_330_tmX", "ttbar-pw-uedown_330_ctt"
            }));
    }
    else
    {
        datasets = datasetBuilder({
          "t-tchan-pw_330rc1_OBU", "tbar-tchan-pw_330rc1_LKv", "t-tWchan-pw_330rc1_ZYV",
          "tbar-tWchan-pw_330rc1_Gll", "t-schan-amcnlo_330_MNE",
          "Wjets-1j-mg_330rc1_feB", "Wjets-2j-mg_330rc1_aKV", "Wjets-3j-mg_330rc1_all",
          "Wjets-4j-mg_330rc1_all", "DY-mg_330rc1_all",
          "WW_330_gaS", "WZ_330_jyL", "ZZ_330_EqH",
          "ttW-lep_330_oSc", "ttW-had_330_Abh", "ttZ-lep_330_WGU", "ttZ-had_330_hfg"});
        
        datasets.splice(datasets.end(), datasetBuilder({
          "A-res-semilep-m400-relW2p5_hSj", "A-res-dilep-m400-relW2p5_iWO",
          "A-res-semilep-m400-relW5_Tdt", "A-res-dilep-m400-relW5_SZC",
          "A-res-semilep-m400-relW10_HvC", "A-res-dilep-m400-relW10_dAW",
          "A-res-semilep-m400-relW25_zMT", "A-res-dilep-m400-relW25_whQ",
          "A-res-semilep-m400-relW50_FFv", "A-res-dilep-m400-relW50_OpI",
          "A-res-semilep-m500-relW2p5_hXS", "A-res-dilep-m500-relW2p5_wyD",
          "A-res-semilep-m500-relW5_UUA", "A-res-dilep-m500-relW5_AwY",
          "A-res-semilep-m500-relW10_HEU", "A-res-dilep-m500-relW10_Yip",
          "A-res-semilep-m500-relW25_Rwg", "A-res-dilep-m500-relW25_TaZ",
          "A-res-semilep-m500-relW50_mIo", "A-res-dilep-m500-relW50_jyT",
          "A-res-semilep-m600-relW2p5_Xbg", "A-res-dilep-m600-relW2p5_NAO",
          "A-res-semilep-m600-relW5_xAs", "A-res-dilep-m600-relW5_AuN",
          "A-res-semilep-m600-relW10_unA", "A-res-dilep-m600-relW10_jfm",
          "A-res-semilep-m600-relW25_ZIp", "A-res-dilep-m600-relW25_Kps",
          "A-res-semilep-m600-relW50_sJU", "A-res-dilep-m600-relW50_oWP",
          "A-res-semilep-m750-relW2p5_dSV", "A-res-dilep-m750-relW2p5_uth",
          "A-res-semilep-m750-relW5_jAk", "A-res-dilep-m750-relW5_qth",
          "A-res-semilep-m750-relW10_gWH", "A-res-dilep-m750-relW10_hOd",
          "A-res-semilep-m750-relW25_ZaF", "A-res-dilep-m750-relW25_iYk",
          "A-res-semilep-m750-relW50_lUC", "A-res-dilep-m750-relW50_Vlv",
          "A-int-semilep-m400-relW2p5_ibc", "A-int-dilep-m400-relW2p5_Joq",
          "A-int-semilep-m400-relW5_NPp", "A-int-dilep-m400-relW5_YLE",
          "A-int-semilep-m400-relW10_jgF", "A-int-dilep-m400-relW10_TPI",
          "A-int-semilep-m400-relW25_Xmd", "A-int-dilep-m400-relW25_JKV",
          "A-int-semilep-m400-relW50_kVF", "A-int-dilep-m400-relW50_BdM",
          "A-int-semilep-m500-relW2p5_YDF", "A-int-dilep-m500-relW2p5_JUJ",
          "A-int-semilep-m500-relW5_uyS", "A-int-dilep-m500-relW5_OyZ",
          "A-int-semilep-m500-relW10_Ijy", "A-int-dilep-m500-relW10_fRi",
          "A-int-semilep-m500-relW25_ZKN", "A-int-dilep-m500-relW25_umJ",
          "A-int-semilep-m500-relW50_mbt", "A-int-dilep-m500-relW50_Hgv",
          "A-int-semilep-m600-relW2p5_ZDm", "A-int-dilep-m600-relW2p5_FBJ",
          "A-int-semilep-m600-relW5_vQy", "A-int-dilep-m600-relW5_MWQ",
          "A-int-semilep-m600-relW10_ugm", "A-int-dilep-m600-relW10_Ghx",
          "A-int-semilep-m600-relW25_HYu", "A-int-dilep-m600-relW25_ySC",
          "A-int-semilep-m600-relW50_lpK", "A-int-dilep-m600-relW50_ZXh",
          "A-int-semilep-m750-relW2p5_Oqn", "A-int-dilep-m750-relW2p5_wNo",
          "A-int-semilep-m750-relW5_IgW", "A-int-dilep-m750-relW5_DxN",
          "A-int-semilep-m750-relW10_Hve", "A-int-dilep-m750-relW10_Aiv",
          "A-int-semilep-m750-relW25_EOV", "A-int-dilep-m750-relW25_DhG",
          "A-int-semilep-m750-relW50_xRp", "A-int-dilep-m750-relW50_MaU",
          "H-res-semilep-m400-relW2p5_332_Msw", "H-res-dilep-m400-relW2p5_332_gEi",
          "H-res-semilep-m400-relW5_332_CTk", "H-res-dilep-m400-relW5_332_pOM",
          "H-res-semilep-m400-relW10_332_jTl", "H-res-dilep-m400-relW10_332_pOY",
          "H-res-semilep-m400-relW25_332_CQc", "H-res-dilep-m400-relW25_332_eTZ",
          "H-res-semilep-m400-relW50_332_nwU", "H-res-dilep-m400-relW50_332_XhS",
          "H-res-semilep-m500-relW2p5_332_clu", "H-res-dilep-m500-relW2p5_332_KXG",
          "H-res-semilep-m500-relW5_332_aQQ", "H-res-dilep-m500-relW5_332_Gkt",
          "H-res-semilep-m500-relW10_332_fFV", "H-res-dilep-m500-relW10_332_fGi",
          "H-res-semilep-m500-relW25_332_HEf", "H-res-dilep-m500-relW25_332_cFM",
          "H-res-semilep-m500-relW50_332_npg", "H-res-dilep-m500-relW50_332_tNA",
          "H-res-semilep-m600-relW2p5_332_Nco", "H-res-dilep-m600-relW2p5_332_DYi",
          "H-res-semilep-m600-relW5_332_BYp", "H-res-dilep-m600-relW5_332_wsF",
          "H-res-semilep-m600-relW10_332_Fqm", "H-res-dilep-m600-relW10_332_FTU",
          "H-res-semilep-m600-relW25_332_RbK", "H-res-dilep-m600-relW25_332_lZC",
          "H-res-semilep-m600-relW50_332_mBi", "H-res-dilep-m600-relW50_332_ljX",
          "H-res-semilep-m750-relW2p5_332_SIz", "H-res-dilep-m750-relW2p5_332_KrU",
          "H-res-semilep-m750-relW5_332_RKo", "H-res-dilep-m750-relW5_332_QMO",
          "H-res-semilep-m750-relW10_332_iMW", "H-res-dilep-m750-relW10_332_oEz",
          "H-res-semilep-m750-relW25_332_liv", "H-res-dilep-m750-relW25_332_agK",
          "H-res-semilep-m750-relW50_332_JeQ", "H-res-dilep-m750-relW50_332_Xws",
          "H-int-semilep-m400-relW2p5_332_Uyv", "H-int-dilep-m400-relW2p5_332_NyN",
          "H-int-semilep-m400-relW5_332_TSO", "H-int-dilep-m400-relW5_332_umr",
          "H-int-semilep-m400-relW10_332_AsU", "H-int-dilep-m400-relW10_332_Rtd",
          "H-int-semilep-m400-relW25_332_MpC", "H-int-dilep-m400-relW25_332_kqh",
          "H-int-semilep-m400-relW50_332_MSk", "H-int-dilep-m400-relW50_332_eux",
          "H-int-semilep-m500-relW2p5_332_pvi", "H-int-dilep-m500-relW2p5_332_Epv",
          "H-int-semilep-m500-relW5_332_zsm", "H-int-dilep-m500-relW5_332_NYZ",
          "H-int-semilep-m500-relW10_332_xBb", "H-int-dilep-m500-relW10_332_Ksa",
          "H-int-semilep-m500-relW25_332_lpS", "H-int-dilep-m500-relW25_332_Mcg",
          "H-int-semilep-m500-relW50_332_DAa", "H-int-dilep-m500-relW50_332_lNS",
          "H-int-semilep-m600-relW2p5_332_rXQ", "H-int-dilep-m600-relW2p5_332_kCu",
          "H-int-semilep-m600-relW5_332_GHN", "H-int-dilep-m600-relW5_332_Kfj",
          "H-int-semilep-m600-relW10_332_nvS", "H-int-dilep-m600-relW10_332_DZL",
          "H-int-semilep-m600-relW25_332_zGL", "H-int-dilep-m600-relW25_332_YIG",
          "H-int-semilep-m600-relW50_332_BJB", "H-int-dilep-m600-relW50_332_kXR",
          "H-int-semilep-m750-relW2p5_332_PUI", "H-int-dilep-m750-relW2p5_332_nWR",
          "H-int-semilep-m750-relW5_332_APZ", "H-int-dilep-m750-relW5_332_CfS",
          "H-int-semilep-m750-relW10_332_uzR", "H-int-dilep-m750-relW10_332_LFo",
          "H-int-semilep-m750-relW25_332_BAb", "H-int-dilep-m750-relW25_332_xrC",
          "H-int-semilep-m750-relW50_332_jUO", "H-int-dilep-m750-relW50_332_wbr"
        }));
        
        if (channel == Channel::Muon)
            datasets.splice(datasets.end(), datasetBuilder({"QCD-mu-15-20_330_JhB",
              "QCD-mu-20-30_330_gEW", "QCD-mu-30-50_330_fUg", "QCD-mu-50-80_330_xQZ",
              "QCD-mu-80-120_330_all", "QCD-mu-120-170_330_all", "QCD-mu-170-300_330_all",
              "QCD-mu-300-470_330_all", "QCD-mu-470-600_330_all", "QCD-mu-600-800_330_all",
              "QCD-mu-800-1000_330_all", "QCD-mu-1000-inf_330_all"}));
        else
            datasets.splice(datasets.end(), datasetBuilder({"QCD-em-20-30_330_Ctz",
              "QCD-em-30-50_330_all", "QCD-em-50-80_330_all", "QCD-em-80-120_330_all",
              "QCD-em-120-170_330_all", "QCD-em-170-300_330_PDg", "QCD-em-300-inf_330_vQG",
              "QCD-bce-20-30_330_YWN", "QCD-bce-30-80_330_uuz", "QCD-bce-80-170-bkp_330_jpV",
              "QCD-bce-170-250_330_MfK", "QCD-bce-250-inf_330_MWY"}));
    }
    
    
    // Triggers
    list<TriggerRange> triggerRanges;
    
    if (channel == Channel::Muon)
        triggerRanges.emplace_back(TriggerRange(0, -1, {"IsoMu24", "IsoTkMu24"}, 35861.523,
          {"IsoMu24", "IsoTkMu24"}));
    else
        triggerRanges.emplace_back(0, -1, "Ele27_WPTight_Gsf", 35861.523,
          "Ele27_WPTight_Gsf");

    
    // Common definition of b-tagging that will be used everywhere
    BTagger const bTagger(BTagger::Algorithm::CMVA, BTagger::WorkingPoint::Medium);
    
    
    // Construct the run manager
    RunManager manager(datasets.begin(), datasets.end());
    
    
    // Register services
    if (sampleGroup != SampleGroup::Data)
        manager.RegisterService(new SystService(systType, systDirection));
    
    BTagWPService *bTagWPService = new BTagWPService("BTagWP_80Xv2.json");
    manager.RegisterService(bTagWPService);
    
    BTagEffService *bTagEffService = new BTagEffService("BTagEff_80Xv3.root");
    bTagEffService->SetDefaultEffLabel("ttbar");
    manager.RegisterService(bTagEffService);
    
    BTagSFService *bTagSFService = new BTagSFService(bTagger, "BTagSF_cMVAv2_80Xv3.csv");
    bTagSFService->SetMeasurement(BTagSFService::Flavour::Bottom, "ttbar");
    bTagSFService->SetMeasurement(BTagSFService::Flavour::Charm, "ttbar");
    bTagSFService->SetMeasurement(BTagSFService::Flavour::Light, "incl");
    manager.RegisterService(bTagSFService);
    
    ostringstream outputNameStream;
    outputNameStream << "output/" << channelText;
    
    if (systType != "None")
        outputNameStream << "_" << systType << "_" <<
          ((systDirection == SystService::VarDirection::Up) ? "up" : "down");
    
    outputNameStream << "/%";
    
    manager.RegisterService(new TFileService(outputNameStream.str()));
    
    
    // Register plugins
    manager.RegisterPlugin(new PECInputData);
    manager.RegisterPlugin(
      BuildPECTriggerFilter((sampleGroup == SampleGroup::Data), triggerRanges));
    
    manager.RegisterPlugin(new PECLeptonReader);
    
    if (channel == Channel::Muon)
        manager.RegisterPlugin(new LeptonFilter("LeptonFilter", Lepton::Flavour::Muon, 26., 2.4));
    else
        manager.RegisterPlugin(new LeptonFilter("LeptonFilter", Lepton::Flavour::Electron,
          30., 2.5));
    
    PECJetMETReader *jetReader = new PECJetMETReader;
    jetReader->SetSelection(20., 2.4);
    manager.RegisterPlugin(jetReader);
    
    JetFilter *jetFilter = new JetFilter(20., bTagger);
    jetFilter->AddSelectionBin(4, -1, 2, -1);
    manager.RegisterPlugin(jetFilter);
    
    manager.RegisterPlugin(new MetFilter(MetFilter::Mode::MtW, 50.));
    manager.RegisterPlugin(new PECPileUpReader);
    
    if (sampleGroup != SampleGroup::Data)
    {
        manager.RegisterPlugin(new PileUpWeight((channel == Channel::Muon) ?
          "Run2016_SingleMuon_v1_finebin.root" : "Run2016_SingleElectron_v1_finebin.root",
          "simPUProfiles_80Xv2.root", 0.05));
        
        if (channel == Channel::Muon)
            manager.RegisterPlugin(new LeptonSFWeight(Lepton::Flavour::Muon,
              "MuonSF_2016_80Xv2.root",
              {"Track", "ID_Tight", "Iso_Tight", "IsoMu24_OR_IsoTkMu24"}));
        else
            manager.RegisterPlugin(new LeptonSFWeight(Lepton::Flavour::Electron,
              "ElectronSF_2016_80Xv2.root",
              {"Track", "CutBasedID_Tight", "Ele27_WPTight_Gsf"}));
        
        BTagWeight *bTagReweighter = new BTagWeight(bTagger);
        bTagReweighter->RequestSystematics();
        manager.RegisterPlugin(bTagReweighter);
        
        
        PECGeneratorReader *generatorReader = new PECGeneratorReader;
        
        if (sampleGroup == SampleGroup::TT)
            generatorReader->RequestAltWeights();
        
        manager.RegisterPlugin(generatorReader);
        
        
        // Dedicated reweighting for signal
        LOSystWeights *scaleWeights = new LOSystWeights(2, "NNPDF30_lo_as_0130");
        scaleWeights->SelectDatasets({"A-.+", "H-.+"});
        manager.RegisterPlugin(scaleWeights);
        
        
        // For SM tt use additional weights
        if (sampleGroup == SampleGroup::TT)
        {
            manager.RegisterPlugin(new PECGenParticleReader());
            
            GenWeightSyst *genWeightSyst = new GenWeightSyst("genWeightVars.json");
            genWeightSyst->NormalizeByMeanWeights(
              "/gridgroup/cms/popov/PECData/2016Charlie/lheWeights_v1.json");
            manager.RegisterPlugin(genWeightSyst);
            
            TopPtWeight *topPtWeights = new TopPtWeight();
            topPtWeights->SelectDatasets({"ttbar-pw_330_all"});
            manager.RegisterPlugin(topPtWeights);
            
            
            manager.RegisterPlugin(new WeightCollector({"LeptonSFWeight", "PileUpWeight",
              "BTagWeight", "GenWeightSyst", "TopPtWeight"}));
        }
        else
            manager.RegisterPlugin(new WeightCollector({"LeptonSFWeight", "PileUpWeight",
              "BTagWeight", "LOSystWeights"}));
    }
    
    
    // Plugin to calculate observables
    manager.RegisterPlugin(new BasicObservables(bTagger));
    
    
    // High-level reconstruction
    TTSemilepRecoRochester *ttRecoPlugin = new TTSemilepRecoRochester;
    ttRecoPlugin->SetLikelihood("TTRecoLikelihood_2016-pt20-v3.root");
    ttRecoPlugin->SetBTagSelection(BTagger::Algorithm::CMVA, bTagWPService->GetThreshold(bTagger),
      false /* both b-quark jets must be tagged */);
    manager.RegisterPlugin(ttRecoPlugin);
    
    
    // Observables exploiting reconstructed top quarks
    manager.RegisterPlugin(new TTObservables);
    
    
    // Event weights
    if (sampleGroup != SampleGroup::Data)
        manager.RegisterPlugin(new DumpWeights("EventWeights"));
    
    
    // Process the datasets
    manager.Process(1);
    
    
    return EXIT_SUCCESS;
}
