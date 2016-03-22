#include <TTSemilepRecoBase.hpp>

#include <PECFwk/core/JetMETReader.hpp>
#include <PECFwk/core/Processor.hpp>

#include <cmath>
#include <stdexcept>


TTSemilepRecoBase::TTSemilepRecoBase(std::string name /*= "TTReco"*/):
    AnalysisPlugin(name),
    minPt(0.), maxAbsEta(std::numeric_limits<double>::infinity())
{}


TTSemilepRecoBase::~TTSemilepRecoBase() noexcept
{}


TTSemilepRecoBase::TTSemilepRecoBase(TTSemilepRecoBase const &src) noexcept:
    AnalysisPlugin(src),
    jetmetPluginName(src.jetmetPluginName), jetmetPlugin(nullptr),
    minPt(src.minPt), maxAbsEta(src.maxAbsEta)
{}


void TTSemilepRecoBase::BeginRun(Dataset const &)
{
    // Save pointer to jet reader
    jetmetPlugin = dynamic_cast<JetMETReader const *>(
      GetMaster().GetPluginBefore(jetmetPluginName, GetName()));
}


Jet const &TTSemilepRecoBase::GetJet(DecayJet type) const
{
    switch (type)
    {
        case DecayJet::bTopLep:
            return jets->at(iBTopLep);
        
        case DecayJet::bTopHad:
            return jets->at(iBTopHad);
        
        case DecayJet::q1TopHad:
            return jets->at(iQ1TopHad);
        
        case DecayJet::q2TopHad:
            return jets->at(iQ2TopHad);
    }
    
    // This must never happen
    throw std::runtime_error("TTSemilepRecoBase::GetJet: Unhandled jet type.");
}


double TTSemilepRecoBase::GetRank() const
{
    return highestRank;
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


bool TTSemilepRecoBase::ProcessEvent()
{
    // Reset data describing the current-best interpretation
    highestRank = -std::numeric_limits<double>::infinity();
    iBTopLep = iBTopHad = iQ1TopHad = iQ2TopHad = -1;
    
    
    // Obtain collection of jets and apply the selection to it
    jets = &jetmetPlugin->GetJets();
    selectedJetIndices.clear();
    
    for (unsigned i = 0; i < jets->size(); ++i)
    {
        if (std::abs(jets->at(i).Eta()) > maxAbsEta)
            continue;
        
        if (jets->at(i).Pt() < minPt)
            break;  // The jet collection is ordered in pt
        
        selectedJetIndices.push_back(i);
    }
    
    unsigned const nSelectedJets = selectedJetIndices.size();
    
    
    // Do not attempt reconstruction if there is not enough jets. However, do not reject the event
    if (selectedJetIndices.size() < 4)
        return true;
    
    
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
                    double const rank = ComputeRank(jets->at(selectedJetIndices.at(iiBTopLepCand)),
                      jets->at(selectedJetIndices.at(iiBTopHadCand)),
                      jets->at(selectedJetIndices.at(iiQ1TopHadCand)),
                      jets->at(selectedJetIndices.at(iiQ2TopHadCand)));
                    
                    if (rank > highestRank)
                    {
                        highestRank = rank;
                        
                        iBTopLep = selectedJetIndices.at(iiBTopLepCand);
                        iBTopHad = selectedJetIndices.at(iiBTopHadCand);
                        iQ1TopHad = selectedJetIndices.at(iiQ1TopHadCand);
                        iQ2TopHad = selectedJetIndices.at(iiQ2TopHadCand);
                    }
                }
            }
        }
    
    
    // Always return true since this plugin does not perform event filtering
    return true;
}
