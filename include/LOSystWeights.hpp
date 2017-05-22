#pragma once

#include <mensura/extensions/EventWeightPlugin.hpp>

#include <initializer_list>
#include <memory>
#include <regex>
#include <string>
#include <vector>
#include <utility>


class GeneratorReader;
namespace LHAPDF{class PDF;};


/**
 * \class LOSystWeights
 * \brief A plugin to compute systematic variation due to renormalization and factorization scales
 * 
 * This reweighting plugin computes weights to reproduce factor 2 variations in renormalization and
 * factorization scales. The procedure is only applicabe for leading-order generators. In addition,
 * the renormalization scale is assumed to be the same for all QCD vertices.
 * 
 * Five weights are evaluated for each processed event. Nominal weight is always set to unity. It
 * is followed by two variations for the renormalization scale, then two variations for the
 * factorization scale. Nominal scale and information about PDF initial state are accessed from a
 * GeneratorReader with a default name "Generator".
 * 
 * User can specify for which datasets weights need to be computed using method SelectDatasets. In
 * the remaining datasets only the nominal weight of unity will be reported.
 */
class LOSystWeights: public EventWeightPlugin
{
public:
    /**
     * \brief Constructs a new reweighting plugin with the given name
     * 
     * The other arguments are the number of QCD vertices and the name of the nominal PDF set.
     */
    LOSystWeights(std::string const &name, unsigned nQCDVert, std::string const &pdfSetName);
    
    /// A short-cut for the above version with a default name "LOSystWeights"
    LOSystWeights(unsigned nQCDVert, std::string const &pdfSetName);
    
    /// Copy constructor
    LOSystWeights(LOSystWeights const &src);
    
    /// Default move constructor
    LOSystWeights(LOSystWeights &&) = default;
    
    /// Assignment operator is deleted
    LOSystWeights &operator=(LOSystWeights const &) = delete;
    
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
    /// Computes alpha_s at the given scale (in GeV)
    static double AlphaS(double scale);
    
    /**
     * \brief Computes systematic weights
     * 
     * Reimplemented from Plugin.
     */
    virtual bool ProcessEvent() override;
    
private:
    /// Name of plugin that provides generator-level weights
    std::string generatorReaderName;
    
    /// Non-owning pointer to plugin that provides generator-level weights
    GeneratorReader const *generatorReader;
    
    /// Masks to choose in which datasets to compute weights
    std::vector<std::regex> datasetMasks;
    
    /// Flag showing whether weights should be computed for the current dataset
    bool processCurDataset;
    
    /**
     * \brief Energy scale is varied by the given factor
     * 
     * "Up" variation is given by scale * scaleVarFactor, "down" one is scale / scaleVarFactor.
     */
    double scaleVarFactor;
    
    /// Number of strong vertices used in reweighting for renomalization scale
    unsigned nQCDVert;
    
    /**
     * \brief Requested PDF set
     * 
     * It is shared among all clones of this plugin
     */
    std::shared_ptr<LHAPDF::PDF const> pdfSet;
};
