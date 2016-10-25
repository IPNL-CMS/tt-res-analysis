#pragma once

#include <mensura/core/AnalysisPlugin.hpp>

#include <TTree.h>

#include <string>
#include <vector>


class PECGeneratorReader;
class PECTriggerFilter;
class TFileService;
class WeightCollector;


/**
 * \class DumpWeights
 * \brief A plugin to save event weights
 * 
 * Creates a TTree with nominal event weight and alternative weights that account for systematic
 * variations. The weights are read from a WeightCollector, and all provided systematic variations
 * are evaluated. In addition to the WeightCollector, this plugin always reads the nominal
 * generator-level weight and trigger weight (i.e. the integrated luminosity). They are
 * incorporated into all stored weights.
 * 
 * Nominal and alternative weights are written into two different branches. The alternative weights
 * are stored as an array, whose size is equal to the number of variations provided by the
 * WeightCollector and thus might depend on the dataset.
 * 
 * This plugin should only be used with simulation.
 */
class DumpWeights: public AnalysisPlugin
{
public:
    /**
     * \brief Constructor
     * 
     * The arguments are the name for the new plugin and the name of a WeightCollector plugin from
     * which event weights will be read (might be empty).
     */
    DumpWeights(std::string const &name, std::string const &weightCollectorName);
    
    /// A short-cut with default name "DumpWeights"
    DumpWeights(std::string const &weightCollectorName);
    
    /// Default move constructor
    DumpWeights(DumpWeights &&) = default;
    
    /// Assignment operator is deleted
    DumpWeights &operator=(DumpWeights const &) = delete;
    
private:
    /// Copy constructor
    DumpWeights(DumpWeights const &src);

public:
    /**
     * \brief Saves pointers to dependencies and sets up the output tree
     * 
     * Reimplemented from Plugin.
     */
    virtual void BeginRun(Dataset const &dataset) override;
    
    /**
     * \brief Creates a newly configured clone
     * 
     * Implemented from Plugin.
     */
    virtual Plugin *Clone() const override;

private:
    /**
     * \brief Saves event weights
     * 
     * Implemented from Plugin.
     */
    virtual bool ProcessEvent() override;

private:
    /// Name of TFileService
    std::string fileServiceName;
    
    /// Non-owning pointer to TFileService
    TFileService const *fileService;
    
    /// Name of trigger filter
    std::string triggerFilterName;
    
    /// Non-owning pointer to the trigger filter
    PECTriggerFilter const *triggerFilter;
    
    /// Name of plugin that provides generator-level weight
    std::string generatorPluginName;
    
    /// Non-owning pointer to plugin that provides generator-level weight
    PECGeneratorReader const *generatorPlugin;
    
    /// Name of a plugin that collects event weights
    std::string weightCollectorName;
    
    /// Non-owning pointer to a plugin that collects event weights
    WeightCollector const *weightCollector;
    
    /**
     * \brief Common event weight for the current dataset
     * 
     * Equals cross section divided by the number of events in the dataset.
     */
    double weightDataset;
    
    /// Non-owning pointer to output tree
    TTree *tree;
    
    // Output buffers
    Float_t weight;
    std::vector<Float_t> systWeights;
};
