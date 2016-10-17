#include <TTObservables.hpp>

#include <mensura/core/Processor.hpp>
#include <mensura/core/ROOTLock.hpp>


TTObservables::TTObservables(std::string const name /*= "TTVars"*/):
    AnalysisPlugin(name),
    fileServiceName("TFileService"), fileService(nullptr),
    ttRecoPluginName("TTReco"), ttRecoPlugin(nullptr)
{}


void TTObservables::BeginRun(Dataset const &)
{
    // Save pointers to services and plugins
    fileService = dynamic_cast<TFileService const *>(GetMaster().GetService(fileServiceName));
    ttRecoPlugin = dynamic_cast<TTSemilepRecoBase const *>(GetDependencyPlugin(ttRecoPluginName));
    
    
    // Set up the output tree
    tree = fileService->Create<TTree>("", GetName().c_str(),
      "Observables relying on tt reconstruction");
    
    ROOTLock::Lock();
    
    tree->Branch("BestRank", &bfBestRank);
    tree->Branch("RecoStatus", &bfRecoStatus);
    
    tree->Branch("MassTopLep", &bfMassTopLep);
    tree->Branch("MassTopHad", &bfMassTopHad);
    tree->Branch("MassWHad", &bfMassWHad);
    tree->Branch("PtTT", &bfPtTT);
    
    tree->Branch("PtTopLep", &bfPtTopLep);
    tree->Branch("PtTopHad", &bfPtTopHad);
    
    tree->Branch("MassTT", &bfMassTT);
    tree->Branch("EtaTT", &bfEtaTT);
    tree->Branch("DRTT", &bfDRTT);
    
    ROOTLock::Unlock();
}


Plugin *TTObservables::Clone() const
{
    return new TTObservables(*this);
}


void TTObservables::SetRecoPluginName(std::string const &pluginName)
{
    ttRecoPluginName = pluginName;
}


bool TTObservables::ProcessEvent()
{
    bfRecoStatus = ttRecoPlugin->GetRecoStatus();
    
    if (bfRecoStatus == 0)
    {
        bfBestRank = ttRecoPlugin->GetRank();
        
        TLorentzVector const p4TopLep = ttRecoPlugin->GetTopLepP4();
        TLorentzVector const p4TopHad = ttRecoPlugin->GetTopHadP4();
        TLorentzVector const p4TT(p4TopLep + p4TopHad);
        
        bfMassTopLep = p4TopLep.M();
        bfMassTopHad = p4TopHad.M();
        bfMassWHad = (ttRecoPlugin->GetJet(TTSemilepRecoBase::DecayJet::q1TopHad).P4() +
          ttRecoPlugin->GetJet(TTSemilepRecoBase::DecayJet::q2TopHad).P4()).M();
        bfPtTT = p4TT.Pt();
        
        bfPtTopLep = p4TopLep.Pt();
        bfPtTopHad = p4TopHad.Pt();
        
        bfMassTT = p4TT.M();
        bfEtaTT = p4TT.Eta();
        bfDRTT = p4TopLep.DeltaR(p4TopHad);
    }
    else
    {
        // Reconstruction has been aborted. Fill variables with some dummy values
        bfBestRank = 0.;
        bfMassTopLep = bfMassTopHad = bfMassWHad = bfPtTT = 0.;
        bfPtTopLep = bfPtTopHad = 0.;
        bfMassTT = bfEtaTT = bfDRTT = 0.;
    }
    
    
    tree->Fill();
    return true;
}
