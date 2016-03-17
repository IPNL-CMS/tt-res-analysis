#include <NuRecoRunI.hpp>

#include <PECFwk/core/LeptonReader.hpp>
#include <PECFwk/core/JetMETReader.hpp>

#include <cassert>
#include <cmath>


NuRecoRunI::NuRecoRunI(std::string name /*= "NuReco"*/):
    NuRecoBase(name)
{}


NuRecoRunI::~NuRecoRunI() noexcept
{}


Plugin *NuRecoRunI::Clone() const
{
    return new NuRecoRunI(*this);
}


bool NuRecoRunI::ProcessEvent()
{
    // Read leptons and MET. Only the leading tight lepton will be used to reconstruct neutrino.
    auto const &leptons = leptonPlugin->GetLeptons();
    
    // Cannot perform a reconstruction when there are no leptons in the event
    if (leptons.size() == 0)
        return true;
    
    auto const &leptonP4 = leptons.front().P4();
    auto const &metP4 = jetmetPlugin->GetMET().P4();
    
    
    // Reset the collection of neutrinos from the previous event
    neutrinos.clear();
    
    
    
    // Reconstruct neutrino. Code is copied from this method [1], with non-essential modifications
    //[1] https://github.com/IPNL-CMS/MttExtractorAnalysis/blob/a198a88bbaccd26c79fef9095ea558416eb2f9e9/plugins/SortingAlgorithm.cc#L9
    TLorentzVector nuP4(metP4);
    
    // Standard quadratic equation for neutrino pz:
    //  a * pz(nu)^2 + b * pz(nu) + c = 0,
    // which originates from the mass constraint m(l+nu) = m(W)
    const double m_w = 80.419;
    double const lambda = (m_w * m_w - leptonP4.M() * leptonP4.M() +
      2 * (nuP4.Px() * leptonP4.Px() + nuP4.Py() * leptonP4.Py())) / (2 * leptonP4.E());
    double const a = 1. - std::pow(leptonP4.Pz() / leptonP4.E(), 2);
    double const b = -2 * (leptonP4.Pz() / leptonP4.E()) * lambda;
    double const c = nuP4.Pt() * nuP4.Pt() - lambda * lambda;
    
    if (a == 0.)
    {
        // This is actually a linear equation (although this should not happen in real life)
        assert(b != 0);
        nuP4.SetPz(-c / b);
        neutrinos.emplace_back(nuP4);
        
        return true;
    }
    
    
    double const discriminant = b * b - 4 * a * c;

    if (discriminant > 0.)
    {
        // The equation has two real-valued solutions. Both corresponding neutrino candidates are
        //recorded
        
        nuP4.SetPz((-b - std::sqrt(discriminant)) / (2 * a));
        neutrinos.emplace_back(nuP4);
        
        nuP4.SetPz((-b + std::sqrt(discriminant)) / (2 * a));
        neutrinos.emplace_back(nuP4);
    }
    else
    {
        // There are no real-valued solutions. Will minimally modify the value of MET, while
        //keeping its direction in the transverse plane, to obtain a zero discriminant
        
        double gamma_x = leptonP4.Px() / std::sqrt(1. + std::pow(nuP4.Py() / nuP4.Px(), 2));
        double gamma_y = leptonP4.Py() / std::sqrt(1. + std::pow(nuP4.Px() / nuP4.Py(), 2));
        
        if (nuP4.Px() < 0.)
            gamma_x = -gamma_x;
        
        if (nuP4.Py() < 0.)
            gamma_y = -gamma_y;
        
        
        // A quadratic equation for the adjusted value of MET:
        //  u * MET^2 + v * MET + w = 0.
        // Presumably, its solution sets the discriminant of the pz(nu) equation to zero
        double const u = std::pow(leptonP4.Pz() / leptonP4.E(), 2) +
          std::pow((gamma_x + gamma_y) / leptonP4.E(), 2) - 1.;
        double const v = std::pow(1. / leptonP4.E(), 2) * (gamma_x + gamma_y) *
          (std::pow(m_w, 2) - std::pow(leptonP4.M(), 2));
        double const w = std::pow((std::pow(m_w, 2) - std::pow(leptonP4.M(), 2)) /
          (2 * leptonP4.E()), 2);

        double const discriminantMET = v * v - 4 * u * w;
        assert(discriminantMET >= 0.);
        
        
        double adjustedMET;
        
        if (u == 0.)
        {
            assert(v != 0.);
            adjustedMET = -w / v;
            
            if (adjustedMET <= 0.)
            {
                // Give up reconstruction
                return true;
            }
        }
        else
        {
            double const met1 = (-v - std::sqrt(discriminantMET)) / (2 * u);
            double const met2 = (-v + std::sqrt(discriminantMET)) / (2 * u);
            
            if (met1 > 0. and met2 > 0.)
            {
                // Choose the solution which is closest to the measured MET
                if (std::fabs(nuP4.Pt() - met1) < std::fabs(nuP4.Pt() - met2))
                    adjustedMET = met1;
                else
                    adjustedMET = met2;
            }
            else
            {
                // Choose positive solution or give up reconstruction if there is none
                if (met1 > 0.)
                    adjustedMET = met1;
                else if (met2 > 0.)
                    adjustedMET = met2;
                else
                    return true;
            }
        }
        
        
        // Modify MET value
        nuP4.SetPtEtaPhiM(adjustedMET, 0., nuP4.Phi(), 0.);
        
        
        // By construction, with the adjusted value of MET the discriminant of the quadratic
        //equation for pz(nu) is zero. Find its solution.
        double const lambdaAdjusted = (m_w * m_w - leptonP4.M() * leptonP4.M() +
          2 * (nuP4.Px() * leptonP4.Px() + nuP4.Py() * leptonP4.Py())) / (2 * leptonP4.E());
        double const aAdjusted = 1. - std::pow(leptonP4.Pz() / leptonP4.E(), 2);
        double const bAdjusted = -2 * (leptonP4.Pz() / leptonP4.E()) * lambdaAdjusted;
        
        nuP4.SetPz(-bAdjusted / (2 * aAdjusted));
        neutrinos.emplace_back(nuP4);
    }
    
    
    
    // Always return true since this method does not perform event filtering
    return true;
}
