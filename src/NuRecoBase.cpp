#include <NuRecoBase.hpp>

#include <PECFwk/core/LeptonReader.hpp>
#include <PECFwk/core/JetMETReader.hpp>
#include <PECFwk/core/Processor.hpp>


NuRecoBase::NuRecoBase(std::string name /*= "NuReco"*/):
    AnalysisPlugin(name),
    leptonPluginName("Leptons"), leptonPlugin(nullptr),
    jetmetPluginName("JetMET"), jetmetPlugin(nullptr)
{}


NuRecoBase::NuRecoBase(NuRecoBase const &src) noexcept:
    AnalysisPlugin(src),
    leptonPluginName(src.leptonPluginName), leptonPlugin(nullptr),
    jetmetPluginName(src.jetmetPluginName), jetmetPlugin(nullptr)
{}


NuRecoBase::~NuRecoBase() noexcept
{}


void NuRecoBase::BeginRun(Dataset const &)
{
    // Save pointers to reader plugins
    leptonPlugin = dynamic_cast<LeptonReader const *>(
      GetMaster().GetPluginBefore(leptonPluginName, GetName()));
    
    jetmetPlugin = dynamic_cast<JetMETReader const *>(
      GetMaster().GetPluginBefore(jetmetPluginName, GetName()));
}


std::vector<Candidate> const &NuRecoBase::GetNeutrinos() const
{
    return neutrinos;
}
