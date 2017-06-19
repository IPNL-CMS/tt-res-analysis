#pragma once

#include <mensura/extensions/EventWeightPlugin.hpp>

#include <array>
#include <string>
#include <regex>


class GenParticleReader;


/**
 * \class TopPtWeight
 * \brief Implements empirical top pt reweighting
 * 
 * Computes nominal weight and two systematic variations for datasets whose source ID matches a
 * mask. The weights are normalized by their mean values before the event selection. Parameters for
 * the reweighting are hard-coded.
 */
class TopPtWeight: public EventWeightPlugin
{
public:
    /// Constructor from a given plugin name
    TopPtWeight(std::string const name = "TopPtWeight");
    
    /// Default copy constructor
    TopPtWeight(TopPtWeight const &) = default;
    
public:
    /**
     * \brief Saves pointers to dependencies and checks the current dataset is to be processed
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
    
    /**
     * \brief Selects datasets for which weights are to be evaluated
     * 
     * The plugin will only compute weights for datasets whose ID match at least one of the
     * provided masks.
     */
    void SelectDatasets(std::initializer_list<std::string> const &masks);

private:
    /**
     * \brief Formula to compute per-event weight without normalization by mean weights
     * 
     * Arguments are the generator-level top quarks' transverse momenta and variations for the two
     * systematic uncertainties, which take values 0, +1, or -1.
     */
    double ComputeTopPtWeight(double pt1, double pt2, int var1, int var2) const;
    
    /**
     * \brief Computes systematic weights
     * 
     * Reimplemented from Plugin.
     */
    virtual bool ProcessEvent() override;
    
private:
    /// Name of plugin that provides generator-level particles
    std::string genParticleReaderName;
    
    /// Non-owning pointer to plugin that provides generator-level particles
    GenParticleReader const *genParticleReader;
    
    /// Masks to choose in which datasets to compute weights
    std::vector<std::regex> datasetMasks;
    
    /// Flag showing whether weights should be computed for the current dataset
    bool processCurDataset;
    
    /// Parameters for the nominal weights
    std::array<double, 2> const nominalParams;
    
    /// Systematic shifs in the parameters for reweighting
    std::array<double, 2> const paramsVar1, paramsVar2;
    
    /**
     * \brief Mean weights before the event selection
     * 
     * Use to renormalize weights. Given in the order (var1, var2) = {(0, 0), (+1, 0), (-1, 0),
     * (0, +1), (0, -1)}.
     */
    std::array<double, 5> const meanWeights;
};
