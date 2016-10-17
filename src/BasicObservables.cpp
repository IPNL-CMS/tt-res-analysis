#include <BasicObservables.hpp>

#include <mensura/core/BTagWPService.hpp>
#include <mensura/core/LeptonReader.hpp>
#include <mensura/core/JetMETReader.hpp>
#include <mensura/core/PileUpReader.hpp>
#include <mensura/core/Processor.hpp>
#include <mensura/core/PhysicsObjects.hpp>
#include <mensura/core/ROOTLock.hpp>

#include <mensura/extensions/TFileService.hpp>
#include <mensura/extensions/WeightCollector.hpp>

#include <mensura/PECReader/PECGeneratorReader.hpp>
#include <mensura/PECReader/PECTriggerFilter.hpp>


BasicObservables::BasicObservables(BTagger const &bTagger_):
    AnalysisPlugin("BasicObservables"),
    bTagger(bTagger_),
    fileServiceName("TFileService"), fileService(nullptr),
    bTagWPServiceName("BTagWP"), bTagWPService(nullptr),
    leptonPluginName("Leptons"), leptonPlugin(nullptr),
    jetmetPluginName("JetMET"), jetmetPlugin(nullptr),
    puPluginName("PileUp"), puPlugin(nullptr),
    triggerFilterName("TriggerFilter"), triggerFilter(nullptr),
    generatorPluginName("Generator"), generatorPlugin(nullptr),
    weightCollectorName("EventWeights"), weightCollector(nullptr)
{}


BasicObservables::BasicObservables(BasicObservables const &src):
    AnalysisPlugin(src),
    bTagger(src.bTagger),
    fileServiceName(src.fileServiceName), fileService(nullptr),
    bTagWPServiceName(src.bTagWPServiceName), bTagWPService(nullptr),
    leptonPluginName(src.leptonPluginName), leptonPlugin(nullptr),
    jetmetPluginName(src.jetmetPluginName), jetmetPlugin(nullptr),
    puPluginName(src.puPluginName), puPlugin(nullptr),
    triggerFilterName(src.triggerFilterName), triggerFilter(nullptr),
    generatorPluginName(src.generatorPluginName), generatorPlugin(nullptr),
    weightCollectorName(src.weightCollectorName), weightCollector(nullptr)
{}


void BasicObservables::BeginRun(Dataset const &dataset)
{
    isMC = dataset.IsMC();
    
    
    // Save pointers to services and readers
    fileService = dynamic_cast<TFileService const *>(GetMaster().GetService(fileServiceName));
    bTagWPService = dynamic_cast<BTagWPService const *>(GetMaster().GetService(bTagWPServiceName));
    
    leptonPlugin = dynamic_cast<LeptonReader const *>(GetDependencyPlugin(leptonPluginName));
    jetmetPlugin = dynamic_cast<JetMETReader const *>(GetDependencyPlugin(jetmetPluginName));
    puPlugin = dynamic_cast<PileUpReader const *>(GetDependencyPlugin(puPluginName));
    
    
    // Save pointer to reweighting plugins
    if (isMC)
    {
        triggerFilter =
          dynamic_cast<PECTriggerFilter const *>(GetDependencyPlugin(triggerFilterName));
        generatorPlugin =
          dynamic_cast<PECGeneratorReader const *>(GetDependencyPlugin(generatorPluginName));
        weightCollector =
          dynamic_cast<WeightCollector const *>(GetDependencyPlugin(weightCollectorName));
    }
    
    
    // Create output tree
    tree = fileService->Create<TTree>("", "BasicVars", "Basic observables");
    
    
    // Assign branch addresses
    ROOTLock::Lock();
    
    tree->Branch("nJet30", &nJet30);
    tree->Branch("nBJet30", &nBJet30);
    
    tree->Branch("Pt_Lep", &Pt_Lep);
    tree->Branch("Eta_Lep", &Eta_Lep);
    
    tree->Branch("Pt_J1", &Pt_J1);
    tree->Branch("Eta_J1", &Eta_J1);
    tree->Branch("Pt_J2", &Pt_J2);
    tree->Branch("Eta_J2", &Eta_J2);
    tree->Branch("Pt_J3", &Pt_J3);
    tree->Branch("Pt_J4", &Pt_J4);
    tree->Branch("Pt_BJ1", &Pt_BJ1);
    
    tree->Branch("CSV_J1", &CSV_J1);
    tree->Branch("CSV_J2", &CSV_J2);
    
    tree->Branch("M_J1J2", &M_J1J2);
    tree->Branch("DR_J1J2", &DR_J1J2);
    tree->Branch("Ht", &Ht);
    tree->Branch("St", &St);
    
    tree->Branch("MET", &MET);
    tree->Branch("Phi_MET", &Phi_MET);
    tree->Branch("MtW", &MtW);
    tree->Branch("nPV", &nPV);
    
    if (isMC)
        tree->Branch("weights", weights, "weights[7]/F");
    
    ROOTLock::Unlock();
    
    
    // Common event weight in this dataset
    auto const &firstFile = dataset.GetFiles().front();
    weightDataset = firstFile.GetWeight();
}


Plugin *BasicObservables::Clone() const
{
    return new BasicObservables(*this);
}


bool BasicObservables::ProcessEvent()
{
    auto const &leptons = leptonPlugin->GetLeptons();
    auto const &met = jetmetPlugin->GetMET();
    
    if (leptons.size() > 0)
    {
        auto const &lep = leptons.front();
        
        Pt_Lep = lep.Pt();
        Eta_Lep = lep.Eta();
        MtW = sqrt(pow(lep.Pt() + met.Pt(), 2) - pow(lep.P4().Px() + met.P4().Px(), 2) -
         pow(lep.P4().Py() + met.P4().Py(), 2));
    }
    else
        Pt_Lep = Eta_Lep = MtW = 0.;
    
    
    auto const &jets = jetmetPlugin->GetJets();
    Pt_J1 = Eta_J1 = Pt_J2 = Eta_J2 = Pt_J3 = Pt_J4 = M_J1J2 = DR_J1J2 = 0.;
    
    if (jets.size() > 0)
    {
        Pt_J1 = jets[0].Pt();
        Eta_J1 = jets[0].Eta();
        CSV_J1 = jets[0].BTag(BTagger::Algorithm::CMVA);
    }
    
    if (jets.size() > 1)
    {
        Pt_J2 = jets[1].Pt();
        Eta_J2 = jets[1].Eta();
        CSV_J2 = jets[1].BTag(BTagger::Algorithm::CMVA);
        
        M_J1J2 = (jets[0].P4() + jets[1].P4()).M();
        DR_J1J2 = jets[0].P4().DeltaR(jets[1].P4());
    }
    
    if (jets.size() > 2)
        Pt_J3 = jets[2].Pt();
    
    if (jets.size() > 3)
        Pt_J4 = jets[3].Pt();
    
    
    nJet30 = nBJet30 = 0;
    Ht = 0.;
    
    for (auto const &j: jets)
    {
        Ht += j.Pt();
        
        if (j.Pt() < 30.)
            continue;
        
        ++nJet30;
        
        if (bTagWPService->IsTagged(bTagger, j))
            ++nBJet30;
    }
    
    Pt_BJ1 = 0.;
    
    for (auto const &j: jets)
        if (bTagWPService->IsTagged(bTagger, j))
        {
            Pt_BJ1 = j.Pt();
            break;
        }
    
    
    MET = met.Pt();
    Phi_MET = met.Phi();
    nPV = puPlugin->GetNumVertices();
    
    St = Ht + Pt_Lep + MET;
    
    
    if (isMC)
    {
        double const w = weightDataset * triggerFilter->GetWeight() *
          generatorPlugin->GetNominalWeight();
        
        weights[0] = w * weightCollector->GetWeight();
        weights[1] = w * weightCollector->GetWeightUp("PileUpWeight", 0);
        weights[2] = w * weightCollector->GetWeightDown("PileUpWeight", 0);
        weights[3] = w * weightCollector->GetWeightUp("BTagWeight", 0);
        weights[4] = w * weightCollector->GetWeightDown("BTagWeight", 0);
        weights[5] = w * weightCollector->GetWeightUp("BTagWeight", 1);
        weights[6] = w * weightCollector->GetWeightDown("BTagWeight", 1);
    }
    
    
    tree->Fill();
    return true;
}
