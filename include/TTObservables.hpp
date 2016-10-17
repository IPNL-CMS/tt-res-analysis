#pragma once

#include <mensura/core/AnalysisPlugin.hpp>

#include <TTSemilepRecoBase.hpp>

#include <mensura/extensions/TFileService.hpp>

#include <TTree.h>


/**
 * \class TTObservables
 * \brief Saves some observables related to reconstructed top quarks
 * 
 * Relies on the presence of a reconstruction plugin with the default name "TTReco".
 */
class TTObservables: public AnalysisPlugin
{
public:
    /// Constructor
    TTObservables(std::string const name = "TTVars");
    
    /// Default copy constructor
    TTObservables(TTObservables const &) = default;
    
    /// Default move constructor
    TTObservables(TTObservables &&) = default;
    
public:
    /**
     * \brief Saves pointers to dependencies and sets up the output tree
     * 
     * Reimplemented from Plugin.
     */
    virtual void BeginRun(Dataset const &) override;
    
    /**
     * \brief Creates a newly configured clone
     * 
     * Implemented from Plugin.
     */
    virtual Plugin *Clone() const override;
    
    /// Specifies name of the plugin that performs tt reconstruction
    void SetRecoPluginName(std::string const &pluginName);
    
private:
    /**
     * \brief Computes observables and fills the tree
     * 
     * Implemented from Plugin.
     */
    virtual bool ProcessEvent() override;
    
private:
    /// Name of TFileService
    std::string fileServiceName;
    
    /// Non-owning pointer to TFileService
    TFileService const *fileService;
    
    /// Name of plugin that reconstructs event under the ttbar hypothesis
    std::string ttRecoPluginName;
    
    /// Non-owning pointer to plugin that reconstructs event under the ttbar hypothesis
    TTSemilepRecoBase const *ttRecoPlugin;
    
    /// Non-owning pointer to output tree
    TTree *tree;
    
    // Output buffers
    Float_t bfBestRank;
    UShort_t bfRecoStatus;
    Float_t bfMassTopLep, bfMassTopHad, bfMassWHad, bfPtTT;
    Float_t bfPtTopLep, bfPtTopHad;
    Float_t bfMassTT, bfEtaTT;
    Float_t bfDRTT;
};
