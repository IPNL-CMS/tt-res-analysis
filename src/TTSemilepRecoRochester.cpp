#include <TTSemilepRecoRochester.hpp>

#include <NuRecoRochester.hpp>

#include <mensura/core/FileInPath.hpp>
#include <mensura/core/JetMETReader.hpp>
#include <mensura/core/LeptonReader.hpp>
#include <mensura/core/Processor.hpp>
#include <mensura/core/ROOTLock.hpp>

#include <TFile.h>

#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>


TTSemilepRecoRochester::TTSemilepRecoRochester(std::string name /*= "TTReco"*/):
    TTSemilepRecoBase(name),
    leptonPluginName("Leptons"), leptonPlugin(nullptr)
{}


TTSemilepRecoRochester::TTSemilepRecoRochester(TTSemilepRecoRochester const &src):
    TTSemilepRecoBase(src),
    leptonPluginName(src.leptonPluginName), leptonPlugin(nullptr),
    likelihoodNeutrino(src.likelihoodNeutrino), likelihoodMass(src.likelihoodMass)
{}


void TTSemilepRecoRochester::BeginRun(Dataset const &dataset)
{
    TTSemilepRecoBase::BeginRun(dataset);
    
    // Save pointers to additional readers
    leptonPlugin = dynamic_cast<LeptonReader const *>(GetDependencyPlugin(leptonPluginName));
    
    
    // Make sure histograms with likelihoods have been provided
    if (not likelihoodNeutrino or not likelihoodMass)
    {
        std::ostringstream message;
        message << "TTSemilepRecoRochester[\"" << GetName() <<
          "\"]::BeginRun: No histograms with likelihoods have been provided.";
        throw std::runtime_error(message.str());
    }
}


Plugin *TTSemilepRecoRochester::Clone() const
{
    return new TTSemilepRecoRochester(*this);
}


Lepton const &TTSemilepRecoRochester::GetLepton() const
{
    if (not lepton)
    {
        std::ostringstream message;
        message << "TTSemilepRecoRochester[\"" << GetName() <<
          "\"]::GetLepton: Current event contains no leptons.";
        throw std::runtime_error(message.str());
    }
    
    return *lepton;
}


Candidate const &TTSemilepRecoRochester::GetNeutrino() const
{
    return neutrino;
}


void TTSemilepRecoRochester::SetLikelihood(std::string const &path,
  std::string const histNeutrinoName /*= "nusolver_chi2_right"*/,
  std::string const histMassName /*= "mWhad_vs_mtophad_right"*/)
{
    // Open file defining likelihoods
    std::string const resolvedPath(FileInPath::Resolve(path));
    
    ROOTLock::Lock();
    TFile inputFile(resolvedPath.c_str());
    ROOTLock::Unlock();
    
    if (inputFile.IsZombie())
    {
        std::ostringstream message;
        message << "TTSemilepRecoRochester[\"" << GetName() << "\"]::SetLikelihood: File \"" <<
          resolvedPath << "\" is not a valid ROOT file.";
        throw std::runtime_error(message.str());
    }
    
    
    // Read the histograms
    ROOTLock::Lock();
    
    likelihoodNeutrino.reset(dynamic_cast<TH1 *>(inputFile.Get(histNeutrinoName.c_str())));
    
    if (not likelihoodNeutrino)
    {
        ROOTLock::Unlock();
        
        std::ostringstream message;
        message << "TTSemilepRecoRochester[\"" << GetName() << "\"]::SetLikelihood: File \"" <<
          resolvedPath << "\" does not contain histogram \"" << histNeutrinoName << "\".";
        throw std::runtime_error(message.str());
    }
    
    likelihoodMass.reset(dynamic_cast<TH2 *>(inputFile.Get(histMassName.c_str())));
    
    if (not likelihoodMass)
    {
        ROOTLock::Unlock();
        
        std::ostringstream message;
        message << "TTSemilepRecoRochester[\"" << GetName() << "\"]::SetLikelihood: File \"" <<
          resolvedPath << "\" does not contain histogram \"" << histMassName << "\".";
        throw std::runtime_error(message.str());
    }
    
    likelihoodNeutrino->SetDirectory(nullptr);
    likelihoodMass->SetDirectory(nullptr);
    
    ROOTLock::Unlock();
    
    
    // Make sure the histograms are normalized to describe probability density
    likelihoodNeutrino->Scale(1. / likelihoodNeutrino->Integral("width"));
    likelihoodMass->Scale(1. / likelihoodMass->Integral("width"));
}


double TTSemilepRecoRochester::ComputeRank(Jet const &bTopLep, Jet const &bTopHad,
  Jet const &q1TopHad, Jet const &q2TopHad)
{
    double logLikelihood = 0.;
    TLorentzVector p4Nu;
    
    
    // Check if the b-quark jet from t -> blv has changed since previous interpretation
    if (&bTopLep == cachedBTopLep)
    {
        // The jet is the same. No need to reconstruct the neutrino again
        logLikelihood += cachedLogLikelihoodNu;
        p4Nu = cachedP4Nu;
    }
    else
    {
        // Reconstruct the neutrino. Skip the current interpretation if it cannot be reconstructed
        NuRecoRochester nuBuilder(&lepton->P4(), &bTopLep.P4());
        
        if (not nuBuilder.IsReconstructable())
            return -std::numeric_limits<double>::infinity();
        
        double nuDistance;
        p4Nu = nuBuilder.GetBest(met->P4().Px(), met->P4().Py(), 1., 1., 0., nuDistance);
        nuDistance = std::sqrt(nuDistance);  // It is actually set to the squared value
        neutrinoReconstructed = true;
        
        
        // Compute (logarithm of) the likelihood for the neutrino distance. If the distance falls
        //into the overflow bin of the histogram, reject the current interpretation
        int bin = likelihoodNeutrino->FindFixBin(nuDistance);
        
        if (likelihoodNeutrino->IsBinOverflow(bin))
            return -std::numeric_limits<double>::infinity();
        
        neutrinoLikelihoodInRange = true;
        
        
        // Update the cached values
        cachedP4Nu = p4Nu;
        cachedLogLikelihoodNu = std::log(likelihoodNeutrino->GetBinContent(bin));
        
        
        // Include neutrino likelihood to the full one
        logLikelihood += cachedLogLikelihoodNu;
    }
    
    
    // Compute masses of hadronically decaying top quark and W boson
    TLorentzVector const p4W(q1TopHad.P4() + q2TopHad.P4());
    double const mW = p4W.M();
    double const mTop = (p4W + bTopHad.P4()).M();
    
    
    // Add (logarithm of) the likelihood for the masses. Reject the current interpretation is at
    //least one of the masses falls into overflow
    int bin = likelihoodMass->FindFixBin(mW, mTop);
    
    if (likelihoodMass->IsBinOverflow(bin))
        return -std::numeric_limits<double>::infinity();
    
    massLikelihoodInRange = true;
    logLikelihood += std::log(likelihoodMass->GetBinContent(bin));
    
    
    // Update best neutrino candidate if needed. At this point GetRank() returns log-likelihood of
    //the best interpretation found so far
    if (logLikelihood > GetRank())
        neutrino.SetP4(p4Nu);
    
    
    return logLikelihood;
}


bool TTSemilepRecoRochester::ProcessEvent()
{
    // Do not attempt reconstruction if the current event contains no leptons
    if (leptonPlugin->GetLeptons().size() == 0)
    {
        lepton = nullptr;
        SetRecoFailure(1);
        return true;
    }
    else
        lepton = &leptonPlugin->GetLeptons().front();
    
    
    // Per-event initialization
    met = &jetmetPlugin->GetMET();
    neutrino.SetPxPyPzE(0., 0., 0., 0.);
    
    neutrinoReconstructed = neutrinoLikelihoodInRange = massLikelihoodInRange = false;
    
    
    // Clear the cache
    cachedBTopLep = nullptr;
    
    
    // Perform jet assigment calling dedicated method from the base class
    PerformJetAssignment(jetmetPlugin->GetJets());
    
    
    // Declare failure of the reconstruction if the best rank is (-inf). This could have happend
    //only if all event interpretations have been rejected
    if (GetRank() == -std::numeric_limits<double>::infinity())
    {
        if (not neutrinoReconstructed)
            SetRecoFailure(2);
        else if (not neutrinoLikelihoodInRange)
            SetRecoFailure(3);
        else if (not massLikelihoodInRange)
            SetRecoFailure(4);
        else
            SetRecoFailure(5);
    }
    
    
    // Always return true since this plugin does not perform event filtering
    return true;
}
