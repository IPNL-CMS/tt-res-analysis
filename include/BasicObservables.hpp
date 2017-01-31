#pragma once

#include <mensura/core/AnalysisPlugin.hpp>

#include <mensura/core/BTagger.hpp>

#include <TTree.h>

#include <string>


class BTagWPService;
class LeptonReader;
class JetMETReader;
class PileUpReader;
class TFileService;


/**
 * \class BasicObservables
 * \brief A plugin to store basic kinematical information
 * 
 * Saves TTree with a number of simple observables.
 */
class BasicObservables: public AnalysisPlugin
{
public:
    /// Constructor
    BasicObservables(BTagger const &bTagger);
    
    /// Default move constructor
    BasicObservables(BasicObservables &&) = default;
    
    /// Assignment operator is deleted
    BasicObservables &operator=(BasicObservables const &) = delete;
    
private:
    /**
     * \brief Copy constructor that produces a newly initialized clone
     * 
     * Behaviour of this copy constructor is appropriate only before processing of the first
     * dataset starts since attributes that are created in BeginRun are not copy. For this reason
     * the copy constructor must not be used in a generic case and made private to prevent this.
     */
    BasicObservables(BasicObservables const &src);

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

private:
    /**
     * \brief Computes representative observables for the current event and fills the tree
     * 
     * Implemented from Plugin.
     */
    virtual bool ProcessEvent() override;

private:
    /// Selected b-tagging algorithm and working point
    BTagger bTagger;
    
    /// Name of TFileService
    std::string fileServiceName;
    
    /// Non-owning pointer to TFileService
    TFileService const *fileService;
    
    /// Name of the service that provides b-tagging working points
    std::string bTagWPServiceName;
    
    /// Non-owning pointer to the service that provides b-tagging working points
    BTagWPService const *bTagWPService;
    
    /// Name of the plugin that produces leptons
    std::string leptonPluginName;
    
    /// Non-owning pointer to the plugin that produces leptons
    LeptonReader const *leptonPlugin;
    
    /// Name of the plugin that produces jets and MET
    std::string jetmetPluginName;
    
    /// Non-owning pointer to the plugin that produces jets and MET
    JetMETReader const *jetmetPlugin;
    
    /// Name of the plugin that reads information about pile-up
    std::string puPluginName;
    
    /// Non-owning pointer to the plugin that reads information about pile-up
    PileUpReader const *puPlugin;
    
    /// Non-owning pointer to output tree
    TTree *tree;
    
    // Output buffers
    Int_t nJet30, nJet20, nBJet30, nBJet20;
    Float_t Pt_Lep, Eta_Lep;
    Float_t Pt_J1, Eta_J1, Pt_J2, Eta_J2, Pt_J3, Pt_J4;
    Float_t Pt_BJ1;
    Float_t bTag_J1, bTag_J2;
    Float_t M_J1J2, DR_J1J2;
    Float_t Ht, St;
    Float_t MET, Phi_MET, MtW;
    Int_t nPV;
    Float_t Rho;
};
