#include <TTSemilepRecoBase.hpp>

#include <mensura/core/JetMETReader.hpp>
#include <mensura/core/Processor.hpp>

#include <cmath>
#include <stdexcept>


TTSemilepRecoBase::TTSemilepRecoBase(std::string name /*= "TTReco"*/):
    AnalysisPlugin(name),
    jetmetPluginName("JetMET"), jetmetPlugin(nullptr),
    minPt(0.), maxAbsEta(std::numeric_limits<double>::infinity())
{}


TTSemilepRecoBase::TTSemilepRecoBase(TTSemilepRecoBase const &src) noexcept:
    AnalysisPlugin(src),
    jetmetPluginName(src.jetmetPluginName), jetmetPlugin(nullptr),
    minPt(src.minPt), maxAbsEta(src.maxAbsEta)
{}


void TTSemilepRecoBase::BeginRun(Dataset const &)
{
    // Save pointer to jet reader
    jetmetPlugin = dynamic_cast<JetMETReader const *>(GetDependencyPlugin(jetmetPluginName));
}


Jet const &TTSemilepRecoBase::GetJet(DecayJet type) const
{
    Jet const *jet = nullptr;
    
    switch (type)
    {
        case DecayJet::bTopLep:
            jet = bTopLep;
            break;
        
        case DecayJet::bTopHad:
            jet = bTopHad;
            break;
        
        case DecayJet::q1TopHad:
            jet = q1TopHad;
            break;
        
        case DecayJet::q2TopHad:
            jet = q2TopHad;
            break;
        
        default:
            // This must never happen
            throw std::runtime_error("TTSemilepRecoBase::GetJet: Unhandled jet type.");
    }
    
    
    if (not jet)
        throw std::runtime_error("TTSemilepRecoBase::GetJet: Requested jet is not available. "
          "This probably means that reconstruction for the current event has been aborted.");
    
    return *jet;
}


double TTSemilepRecoBase::GetRank() const
{
    return highestRank;
}


unsigned TTSemilepRecoBase::GetRecoStatus() const
{
    return recoStatus;
}


TLorentzVector TTSemilepRecoBase::GetTopLepP4() const
{
    return GetLepton().P4() + GetNeutrino().P4() + GetJet(DecayJet::bTopLep).P4();
}


TLorentzVector TTSemilepRecoBase::GetTopHadP4() const
{
    return GetJet(DecayJet::bTopHad).P4() + GetJet(DecayJet::q1TopHad).P4() +
      GetJet(DecayJet::q2TopHad).P4();
}


void TTSemilepRecoBase::SetJetSelection(double minPt_,
  double maxAbsEta_ /*= std::numeric_limits<double>::infinity()*/)
{
    minPt = minPt_;
    maxAbsEta = maxAbsEta_;
}


void TTSemilepRecoBase::PerformJetAssignment(std::vector<Jet> const &jets)
{
    // Reset data describing the current-best interpretation
    highestRank = -std::numeric_limits<double>::infinity();
    bTopLep = bTopHad = q1TopHad = q2TopHad = nullptr;
    
    
    // Save a pointer to the collection of jets and apply the selection to it
    selectedJetIndices.clear();
    
    for (unsigned i = 0; i < jets.size(); ++i)
    {
        if (std::abs(jets.at(i).Eta()) > maxAbsEta)
            continue;
        
        if (jets.at(i).Pt() < minPt)
            break;  // The jet collection is ordered in pt
        
        selectedJetIndices.push_back(i);
    }
    
    unsigned const nSelectedJets = selectedJetIndices.size();
    
    
    // Do not attempt reconstruction if there is not enough jets
    if (nSelectedJets < 4)
    {
        recoStatus = 1;
        return;
    }
    
    
    // Loop over all possible ways of jet assignment to find the best one
    for (unsigned iiBTopLepCand = 0; iiBTopLepCand < nSelectedJets; ++iiBTopLepCand)
        for (unsigned iiBTopHadCand = 0; iiBTopHadCand < nSelectedJets; ++iiBTopHadCand)
        {
            if (iiBTopLepCand == iiBTopHadCand)
                continue;
            
            for (unsigned iiQ1TopHadCand = 0; iiQ1TopHadCand < nSelectedJets; ++iiQ1TopHadCand)
            {
                if (iiQ1TopHadCand == iiBTopLepCand or iiQ1TopHadCand == iiBTopHadCand)
                    continue;
                
                // When looping for the subleading light-flavour jet, take into account that the
                //collection is still ordered in jet pt
                for (unsigned iiQ2TopHadCand = iiQ1TopHadCand + 1; iiQ2TopHadCand < nSelectedJets;
                  ++iiQ2TopHadCand)
                {
                    // An interpretation has been constructed. Evaluate it
                    double const rank = ComputeRank(jets.at(selectedJetIndices.at(iiBTopLepCand)),
                      jets.at(selectedJetIndices.at(iiBTopHadCand)),
                      jets.at(selectedJetIndices.at(iiQ1TopHadCand)),
                      jets.at(selectedJetIndices.at(iiQ2TopHadCand)));
                    
                    if (rank > highestRank)
                    {
                        highestRank = rank;
                        
                        bTopLep = &jets.at(selectedJetIndices.at(iiBTopLepCand));
                        bTopHad = &jets.at(selectedJetIndices.at(iiBTopHadCand));
                        q1TopHad = &jets.at(selectedJetIndices.at(iiQ1TopHadCand));
                        q2TopHad = &jets.at(selectedJetIndices.at(iiQ2TopHadCand));
                    }
                }
            }
        }
    
    
    recoStatus = 0;
}


void TTSemilepRecoBase::SetRecoFailure(unsigned code)
{
    recoStatus = code;
}


bool TTSemilepRecoBase::ProcessEvent()
{
    PerformJetAssignment(jetmetPlugin->GetJets());
    
    // Always return true since this plugin does not perform event filtering
    return true;
}
