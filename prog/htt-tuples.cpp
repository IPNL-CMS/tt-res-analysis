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
    DatasetBuilder datasetBuilder("/gridgroup/cms/popov/PECData/2016Delta/samples_v1.json");
    
    if (sampleGroup == SampleGroup::Data)
    {
        if (channel == Channel::Muon)
            datasets = datasetBuilder({"SingleMuon-Run2016_333_all"});
        else
            datasets = datasetBuilder({"SingleElectron-Run2016_333_all"});
    }
    else if (sampleGroup == SampleGroup::TT)
    {
        datasets = datasetBuilder({"ttbar-pw_333_all"});
        
        if (systType == "None")
            datasets.splice(datasets.end(), datasetBuilder({
              "ttbar-pw-isrup_333_Jic", "ttbar-pw-isrdown_333_all",
              "ttbar-pw-fsrup_333_all", "ttbar-pw-fsrdown_333_all",
              "ttbar-pw-hdampup_333_all", "ttbar-pw-hdampdown_333_all",
              "ttbar-pw-m1755_333_all", "ttbar-pw-m1695_333_all",
              "ttbar-pw-ueup_333_all", "ttbar-pw-uedown_333_all"
            }));
    }
    else
    {
        datasets = datasetBuilder({
          "t-tchan-pw_333_ecs", "tbar-tchan-pw_333_MWZ", "t-tWchan-pw_333_WoS",
          "tbar-tWchan-pw_333_eGC", "t-schan-amcnlo_333_ErJ",
          "Wjets-1j-mg_333_JKN", "Wjets-2j-mg_333_QrW", "Wjets-3j-mg_333_all",
          "Wjets-4j-mg_333_all", "DY-mg_333_all",
          "WW_333_qpN", "WZ_333_qsl", "ZZ_333_ydJ",
          "ttW-lep_333_QUK", "ttW-had_333_mya", "ttZ-lep_333_QWi", "ttZ-had_333_Vpe"});
        
        datasets.splice(datasets.end(), datasetBuilder({
          "A-res-semilep-m400-relW2p5_333_ckX", "A-res-dilep-m400-relW2p5_333_AJR",
          "A-res-semilep-m400-relW5_333_Gcq", "A-res-dilep-m400-relW5_333_dEC",
          "A-res-semilep-m400-relW10_333_qpo", "A-res-dilep-m400-relW10_333_rqN",
          "A-res-semilep-m400-relW25_333_Jxx", "A-res-dilep-m400-relW25_333_UJg",
          "A-res-semilep-m400-relW50_333_vTE", "A-res-dilep-m400-relW50_333_hzH",
          "A-res-semilep-m500-relW2p5_333_Ory", "A-res-dilep-m500-relW2p5_333_kxt",
          "A-res-semilep-m500-relW5_333_LYW", "A-res-dilep-m500-relW5_333_LFD",
          "A-res-semilep-m500-relW10_333_FXK", "A-res-dilep-m500-relW10_333_wph",
          "A-res-semilep-m500-relW25_333_aZp", "A-res-dilep-m500-relW25_333_sKi",
          "A-res-semilep-m500-relW50_333_gEA", "A-res-dilep-m500-relW50_333_gmN",
          "A-res-semilep-m600-relW2p5_333_ZBC", "A-res-dilep-m600-relW2p5_333_ApB",
          "A-res-semilep-m600-relW5_333_hOz", "A-res-dilep-m600-relW5_333_qln",
          "A-res-semilep-m600-relW10_333_Ibj", "A-res-dilep-m600-relW10_333_alD",
          "A-res-semilep-m600-relW25_333_mlb", "A-res-dilep-m600-relW25_333_xbk",
          "A-res-semilep-m600-relW50_333_Wto", "A-res-dilep-m600-relW50_333_Qpk",
          "A-res-semilep-m750-relW2p5_333_IWi", "A-res-dilep-m750-relW2p5_333_TGM",
          "A-res-semilep-m750-relW5_333_DKc", "A-res-dilep-m750-relW5_333_huq",
          "A-res-semilep-m750-relW10_333_iPi", "A-res-dilep-m750-relW10_333_GIC",
          "A-res-semilep-m750-relW25_333_ijc", "A-res-dilep-m750-relW25_333_AsE",
          "A-res-semilep-m750-relW50_333_OEj", "A-res-dilep-m750-relW50_333_UxA",
          "A-int-semilep-m400-relW2p5_333_LEx", "A-int-dilep-m400-relW2p5_333_ryF",
          "A-int-semilep-m400-relW5_333_JyF", "A-int-dilep-m400-relW5_333_XQz",
          "A-int-semilep-m400-relW10_333_WZQ", "A-int-dilep-m400-relW10_333_FLD",
          "A-int-semilep-m400-relW25_333_Sqb", "A-int-dilep-m400-relW25_333_SRA",
          "A-int-semilep-m400-relW50_333_qPD", "A-int-dilep-m400-relW50_333_NzY",
          "A-int-semilep-m500-relW2p5_333_HRX", "A-int-dilep-m500-relW2p5_333_xIs",
          "A-int-semilep-m500-relW5_333_oix", "A-int-dilep-m500-relW5_333_itg",
          "A-int-semilep-m500-relW10_333_EMU", "A-int-dilep-m500-relW10_333_tcn",
          "A-int-semilep-m500-relW25_333_DzU", "A-int-dilep-m500-relW25_333_pQh",
          "A-int-semilep-m500-relW50_333_gvO", "A-int-dilep-m500-relW50_333_VzJ",
          "A-int-semilep-m600-relW2p5_333_dYG", "A-int-dilep-m600-relW2p5_333_Ffd",
          "A-int-semilep-m600-relW5_333_AHd", "A-int-dilep-m600-relW5_333_QBA",
          "A-int-semilep-m600-relW10_333_FXb", "A-int-dilep-m600-relW10_333_Yow",
          "A-int-semilep-m600-relW25_333_kjf", "A-int-dilep-m600-relW25_333_AgQ",
          "A-int-semilep-m600-relW50_333_vdk", "A-int-dilep-m600-relW50_333_dSy",
          "A-int-semilep-m750-relW2p5_333_lCs", "A-int-dilep-m750-relW2p5_333_ADQ",
          "A-int-semilep-m750-relW5_333_jJe", "A-int-dilep-m750-relW5_333_LvY",
          "A-int-semilep-m750-relW10_333_yPO", "A-int-dilep-m750-relW10_333_Trw",
          "A-int-semilep-m750-relW25_333_Qdt", "A-int-dilep-m750-relW25_333_lzj",
          "A-int-semilep-m750-relW50_333_tST", "A-int-dilep-m750-relW50_333_ZRR",
          "H-res-semilep-m400-relW2p5_333_ORF", "H-res-dilep-m400-relW2p5_333_xvF",
          "H-res-semilep-m400-relW5_333_yri", "H-res-dilep-m400-relW5_333_zzO",
          "H-res-semilep-m400-relW10_333_rwm", "H-res-dilep-m400-relW10_333_hDn",
          "H-res-semilep-m400-relW25_333_Jkz", "H-res-dilep-m400-relW25_333_HkF",
          "H-res-semilep-m400-relW50_333_VLG", "H-res-dilep-m400-relW50_333_oGj",
          "H-res-semilep-m500-relW2p5_333_klL", "H-res-dilep-m500-relW2p5_333_oFU",
          "H-res-semilep-m500-relW5_333_vui", "H-res-dilep-m500-relW5_333_Ahe",
          "H-res-semilep-m500-relW10_333_yLr", "H-res-dilep-m500-relW10_333_dzy",
          "H-res-semilep-m500-relW25_333_hby", "H-res-dilep-m500-relW25_333_jVD",
          "H-res-semilep-m500-relW50_333_rpY", "H-res-dilep-m500-relW50_333_jPK",
          "H-res-semilep-m600-relW2p5_333_Ubo", "H-res-dilep-m600-relW2p5_333_Mrl",
          "H-res-semilep-m600-relW5_333_eFQ", "H-res-dilep-m600-relW5_333_PLt",
          "H-res-semilep-m600-relW10_333_anO", "H-res-dilep-m600-relW10_333_GDC",
          "H-res-semilep-m600-relW25_333_FYm", "H-res-dilep-m600-relW25_333_Lgk",
          "H-res-semilep-m600-relW50_333_tVh", "H-res-dilep-m600-relW50_333_yjt",
          "H-res-semilep-m750-relW2p5_333_hCO", "H-res-dilep-m750-relW2p5_333_brb",
          "H-res-semilep-m750-relW5_333_zmf", "H-res-dilep-m750-relW5_333_RgV",
          "H-res-semilep-m750-relW10_333_xoP", "H-res-dilep-m750-relW10_333_ybn",
          "H-res-semilep-m750-relW25_333_eWp", "H-res-dilep-m750-relW25_333_VQw",
          "H-res-semilep-m750-relW50_333_uwo", "H-res-dilep-m750-relW50_333_AzS",
          "H-int-semilep-m400-relW2p5_333_DaG", "H-int-dilep-m400-relW2p5_333_PyC",
          "H-int-semilep-m400-relW5_333_ebU", "H-int-dilep-m400-relW5_333_oCG",
          "H-int-semilep-m400-relW10_333_NOl", "H-int-dilep-m400-relW10_333_tlu",
          "H-int-semilep-m400-relW25_333_JvM", "H-int-dilep-m400-relW25_333_RBI",
          "H-int-semilep-m400-relW50_333_OcP", "H-int-dilep-m400-relW50_333_Caw",
          "H-int-semilep-m500-relW2p5_333_QZi", "H-int-dilep-m500-relW2p5_333_MLm",
          "H-int-semilep-m500-relW5_333_gSL", "H-int-dilep-m500-relW5_333_roN",
          "H-int-semilep-m500-relW10_333_bOG", "H-int-dilep-m500-relW10_333_OhM",
          "H-int-semilep-m500-relW25_333_tEG", "H-int-dilep-m500-relW25_333_FbC",
          "H-int-semilep-m500-relW50_333_CsU", "H-int-dilep-m500-relW50_333_YPp",
          "H-int-semilep-m600-relW2p5_333_LOk", "H-int-dilep-m600-relW2p5_333_CVz",
          "H-int-semilep-m600-relW5_333_cae", "H-int-dilep-m600-relW5_333_bUg",
          "H-int-semilep-m600-relW10_333_jZu", "H-int-dilep-m600-relW10_333_cOO",
          "H-int-semilep-m600-relW25_333_Wyh", "H-int-dilep-m600-relW25_333_XEx",
          "H-int-semilep-m600-relW50_333_Urx", "H-int-dilep-m600-relW50_333_Qbn",
          "H-int-semilep-m750-relW2p5_333_sYp", "H-int-dilep-m750-relW2p5_333_yRA",
          "H-int-semilep-m750-relW5_333_iIp", "H-int-dilep-m750-relW5_333_yUA",
          "H-int-semilep-m750-relW10_333_sIm", "H-int-dilep-m750-relW10_333_bwh",
          "H-int-semilep-m750-relW25_333_yUk", "H-int-dilep-m750-relW25_333_DPk",
          "H-int-semilep-m750-relW50_333_BxK", "H-int-dilep-m750-relW50_333_tEs"
        }));
        
        // if (channel == Channel::Muon)
        //     datasets.splice(datasets.end(), datasetBuilder({"QCD-mu-15-20_330_JhB",
        //       "QCD-mu-20-30_330_gEW", "QCD-mu-30-50_330_fUg", "QCD-mu-50-80_330_xQZ",
        //       "QCD-mu-80-120_330_all", "QCD-mu-120-170_330_all", "QCD-mu-170-300_330_all",
        //       "QCD-mu-300-470_330_all", "QCD-mu-470-600_330_all", "QCD-mu-600-800_330_all",
        //       "QCD-mu-800-1000_330_all", "QCD-mu-1000-inf_330_all"}));
        // else
        //     datasets.splice(datasets.end(), datasetBuilder({"QCD-em-20-30_330_Ctz",
        //       "QCD-em-30-50_330_all", "QCD-em-50-80_330_all", "QCD-em-80-120_330_all",
        //       "QCD-em-120-170_330_all", "QCD-em-170-300_330_PDg", "QCD-em-300-inf_330_vQG",
        //       "QCD-bce-20-30_330_YWN", "QCD-bce-30-80_330_uuz", "QCD-bce-80-170-bkp_330_jpV",
        //       "QCD-bce-170-250_330_MfK", "QCD-bce-250-inf_330_MWY"}));
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
        {
            manager.RegisterPlugin(new LeptonSFWeight("TriggerSFWeight", Lepton::Flavour::Muon,
              "MuonSF_2016_80Xv2.root", {"IsoMu24_OR_IsoTkMu24"}));
            manager.RegisterPlugin(new LeptonSFWeight("LeptonSFWeight", Lepton::Flavour::Muon,
              "MuonSF_2016_80Xv2.root", {"Track", "ID_Tight", "Iso_Tight"}));
        }
        else
        {
            manager.RegisterPlugin(new LeptonSFWeight("TriggerSFWeight", Lepton::Flavour::Electron,
              "ElectronSF_2016_80Xv2.root", {"Ele27_WPTight_Gsf"}));
            manager.RegisterPlugin(new LeptonSFWeight("LeptonSFWeight", Lepton::Flavour::Electron,
              "ElectronSF_2016_80Xv2.root", {"Track", "CutBasedID_Tight"}));
        }
        
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
              "/gridgroup/cms/popov/PECData/2016Delta/lheWeights_v1.json");
            manager.RegisterPlugin(genWeightSyst);
            
            TopPtWeight *topPtWeights = new TopPtWeight();
            topPtWeights->SelectDatasets({"ttbar-pw[-_].*"});
            manager.RegisterPlugin(topPtWeights);
            
            
            manager.RegisterPlugin(new WeightCollector({"TriggerSFWeight", "LeptonSFWeight",
              "PileUpWeight", "BTagWeight", "GenWeightSyst", "TopPtWeight"}));
        }
        else
            manager.RegisterPlugin(new WeightCollector({"TriggerSFWeight", "LeptonSFWeight",
              "PileUpWeight", "BTagWeight", "LOSystWeights"}));
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
    manager.Process(16);
    
    
    return EXIT_SUCCESS;
}
