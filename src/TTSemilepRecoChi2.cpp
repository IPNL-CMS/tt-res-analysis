#include <TTSemilepRecoChi2.hpp>

#include <NuRecoBase.hpp>

#include <mensura/core/LeptonReader.hpp>
#include <mensura/core/JetMETReader.hpp>
#include <mensura/core/Processor.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>



// Functions defining various types of chi2 terms
double exprMassTopLep(Lepton const &l, Candidate const &nu, Jet const &bTopLep,  Jet const &,
  Jet const &, Jet const &)
{
    return (l.P4() + nu.P4() + bTopLep.P4()).M();
}


double exprMassTopHad(Lepton const &, Candidate const &, Jet const &, Jet const &bTopHad,
  Jet const &q1TopHad, Jet const &q2TopHad)
{
    return (bTopHad.P4() + q1TopHad.P4() + q2TopHad.P4()).M();
}


double exprMassWHad(Lepton const &, Candidate const &, Jet const &, Jet const &,
  Jet const &q1TopHad, Jet const &q2TopHad)
{
    return (q1TopHad.P4() + q2TopHad.P4()).M();
}


double exprPtTT(Lepton const &l, Candidate const &nu, Jet const &bTopLep, Jet const &bTopHad,
  Jet const &q1TopHad, Jet const &q2TopHad)
{
    return (l.P4() + nu.P4() + bTopLep.P4() + bTopHad.P4() + q1TopHad.P4() + q2TopHad.P4()).Pt();
}



double TTSemilepRecoChi2::Chi2Term::Eval(Lepton const &l, Candidate const &nu, Jet const &bTopLep,
  Jet const &bTopHad, Jet const &q1TopHad, Jet const &q2TopHad) const
{
    double const x = expression(l, nu, bTopLep, bTopHad, q1TopHad, q2TopHad);
    return std::pow((x - mean) / variance, 2);
}



TTSemilepRecoChi2::TTSemilepRecoChi2(std::string name /*= "TTReco"*/):
    TTSemilepRecoBase(name),
    leptonPluginName("Leptons"), leptonPlugin(nullptr),
    nuRecoPluginName("NuReco"), nuRecoPlugin(nullptr)
{}


TTSemilepRecoChi2::~TTSemilepRecoChi2() noexcept
{}


TTSemilepRecoChi2::TTSemilepRecoChi2(TTSemilepRecoChi2 const &src) noexcept:
    TTSemilepRecoBase(src),
    leptonPluginName(src.leptonPluginName), leptonPlugin(nullptr),
    nuRecoPluginName(src.nuRecoPluginName), nuRecoPlugin(nullptr),
    chi2Terms(src.chi2Terms)
{}


void TTSemilepRecoChi2::AddChi2Term(Expression expression, double mean, double variance)
{
    switch (expression)
    {
        case Expression::MassTopLep:
            chi2Terms.emplace_back(Chi2Term{exprMassTopLep, mean, variance});
            break;
        
        case Expression::MassTopHad:
            chi2Terms.emplace_back(Chi2Term{exprMassTopHad, mean, variance});
            break;
        
        case Expression::MassWHad:
            chi2Terms.emplace_back(Chi2Term{exprMassWHad, mean, variance});
            break;
        
        case Expression::PtTT:
            chi2Terms.emplace_back(Chi2Term{exprPtTT, mean, variance});
            break;
        
        default:
            throw std::runtime_error("TTSemilepRecoChi2::AddChi2Term: Unhandled expression type.");
    }
}


void TTSemilepRecoChi2::BeginRun(Dataset const &)
{
    // Save pointers to readers and plugin for neutrino reconstruction
    leptonPlugin = dynamic_cast<LeptonReader const *>(
      GetMaster().GetPluginBefore(leptonPluginName, GetName()));
    jetmetPlugin = dynamic_cast<JetMETReader const *>(
      GetMaster().GetPluginBefore(jetmetPluginName, GetName()));
    nuRecoPlugin = dynamic_cast<NuRecoBase const *>(
      GetMaster().GetPluginBefore(nuRecoPluginName, GetName()));
}


Plugin *TTSemilepRecoChi2::Clone() const
{
    return new TTSemilepRecoChi2(*this);
}


Lepton const &TTSemilepRecoChi2::GetLepton() const
{
    auto const &leptons = leptonPlugin->GetLeptons();
    
    if (leptons.size() == 0)
        throw std::runtime_error("TTSemilepRecoChi2::GetLepton: Current event contains no "
          "leptons.");
    
    return leptons.front();
}


Candidate const &TTSemilepRecoChi2::GetNeutrino() const
{
    if (not bestNu)
        throw std::runtime_error("TTSemilepRecoChi2::GetNeutrino: No neutrino has been "
          "reconstructed in the current event.");
    
    return *bestNu;
}


double TTSemilepRecoChi2::ComputeRank(Jet const &bTopLep, Jet const &bTopHad, Jet const &q1TopHad,
  Jet const &q2TopHad)
{
    // There might be several solutions for neutrino. Loop over all of them and find the minimal
    //chi^2 for the probed jet assignment
    double minChi2CurInterp = std::numeric_limits<double>::infinity();
    
    for (auto const &nu: nuRecoPlugin->GetNeutrinos())
    {
        // Calculate the chi^2 for the current neutrino solution. Only the first lepton is used
        double chi2 = 0.;
        
        for (auto const &term: chi2Terms)
            chi2 += term.Eval(leptonPlugin->GetLeptons().front(), nu, bTopLep, bTopHad, q1TopHad,
              q2TopHad);
        
        
        // Update the minimal chi^2 in the current interpretation
        if (chi2 < minChi2CurInterp)
            minChi2CurInterp = chi2;
        
        
        // Update the best neutrino solution if needed
        if (chi2 < minChi2)
        {
            minChi2 = chi2;
            bestNu = &nu;
        }
    }
    
    
    // Most likely interpretation should have the largest rank, thus use negative chi^2 as the
    //rank
    return -minChi2CurInterp;
}


bool TTSemilepRecoChi2::ProcessEvent()
{
    // Reset data describing best solution for neutrino
    minChi2 = std::numeric_limits<double>::infinity();
    bestNu = nullptr;
    
    
    // Do not attempt reconstruction if the current event contains no leptons or no reconstructed
    // neutrinos. However, do not reject the event
    if (leptonPlugin->GetLeptons().size() == 0 or nuRecoPlugin->GetNeutrinos().size() == 0)
    {
        SetRecoFailure(1);
        return true;
    }
    
    
    // Perform jet assigment calling dedicated method from the base class
    PerformJetAssignment(jetmetPlugin->GetJets());
    
    
    // Always return true since this plugin does not filter events
    return true;
}
