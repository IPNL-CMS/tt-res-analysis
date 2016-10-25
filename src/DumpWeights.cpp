#include <DumpWeights.hpp>

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

#include <sstream>
#include <stdexcept>


using namespace std::string_literals;


DumpWeights::DumpWeights(std::string const &name, std::string const &weightCollectorName):
    AnalysisPlugin(name),
    fileServiceName("TFileService"), fileService(nullptr),
    triggerFilterName("TriggerFilter"), triggerFilter(nullptr),
    generatorPluginName("Generator"), generatorPlugin(nullptr),
    weightCollectorName(weightCollectorName), weightCollector(nullptr)
{}


DumpWeights::DumpWeights(std::string const &weightCollectorName):
    DumpWeights("DumpWeights", weightCollectorName)
{}


DumpWeights::DumpWeights(DumpWeights const &src):
    AnalysisPlugin(src),
    fileServiceName(src.fileServiceName), fileService(nullptr),
    triggerFilterName(src.triggerFilterName), triggerFilter(nullptr),
    generatorPluginName(src.generatorPluginName), generatorPlugin(nullptr),
    weightCollectorName(src.weightCollectorName), weightCollector(nullptr),
    systWeights(src.systWeights)
{}


void DumpWeights::BeginRun(Dataset const &dataset)
{
    if (not dataset.IsMC())
    {
        std::ostringstream message;
        message << "DumpWeights[\"" << GetName() << "\"]::BeginRun: The current dataset is data, "
          "but this plugin should only be used with simulation.";
        throw std::runtime_error(message.str());
    }
    
    
    // Save pointers to services and other plugins
    fileService = dynamic_cast<TFileService const *>(GetMaster().GetService(fileServiceName));
    triggerFilter = dynamic_cast<PECTriggerFilter const *>(GetDependencyPlugin(triggerFilterName));
    generatorPlugin =
      dynamic_cast<PECGeneratorReader const *>(GetDependencyPlugin(generatorPluginName));
    
    if (weightCollectorName != "")
        weightCollector =
          dynamic_cast<WeightCollector const *>(GetDependencyPlugin(weightCollectorName));
    
    
    // Adjust the size of the array to store alternative weights
    unsigned nSystWeights = 0;
    
    if (weightCollectorName != "")
        for (unsigned iPlugin = 0; iPlugin < weightCollector->GetNumPlugins(); ++iPlugin)
            nSystWeights += 2 * weightCollector->GetPlugin(iPlugin)->GetNumVariations();
    
    systWeights.resize(nSystWeights);
    
    
    // Create output tree
    tree = fileService->Create<TTree>("", "Weights", "Nominal and alternative weights");
    
    
    // Assign branch addresses
    ROOTLock::Lock();
    
    tree->Branch("weight", &weight);
    tree->Branch("systWeights", systWeights.data(),
      ("systWeights["s + std::to_string(systWeights.size()) + "]/F").c_str());
    
    ROOTLock::Unlock();
    
    
    // Common event weight in this dataset
    auto const &firstFile = dataset.GetFiles().front();
    weightDataset = firstFile.GetWeight();
}


Plugin *DumpWeights::Clone() const
{
    return new DumpWeights(*this);
}


bool DumpWeights::ProcessEvent()
{
    double const w = weightDataset * triggerFilter->GetWeight() *
      generatorPlugin->GetNominalWeight();
    
    
    if (weightCollector)
    {
        weight = w * weightCollector->GetWeight();
        unsigned curWeightIndex = 0;
        
        for (unsigned iPlugin = 0; iPlugin < weightCollector->GetNumPlugins(); ++iPlugin)
        {
            EventWeightPlugin const *plugin = weightCollector->GetPlugin(iPlugin);
            
            for (unsigned iVar = 0; iVar < plugin->GetNumVariations(); ++iVar)
            {
                systWeights.at(curWeightIndex) = w * weightCollector->GetWeightUp(iPlugin, iVar);
                systWeights.at(curWeightIndex + 1) =
                  w * weightCollector->GetWeightDown(iPlugin, iVar);
                curWeightIndex += 2;
            }
        }
    }
    else
        weight = w;
    
    tree->Fill();
    return true;
}
