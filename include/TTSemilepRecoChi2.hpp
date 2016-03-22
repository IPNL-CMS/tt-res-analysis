#pragma once

#include <TTSemilepRecoBase.hpp>

#include <PECFwk/core/PhysicsObjects.hpp>

#include <vector>


class LeptonReader;
class NuRecoBase;


/**
 * \class TTSemilepRecoChi2
 * \brief A plugin to perform jet assignment in semileptonic ttbar events using customized chi^2
 * 
 * The plugin builds on top of TTSemilepRecoBase, implementing a chi^2 figure of merit to rank
 * alternative event interpretations. The chi^2 is constructed from summands provided by user with
 * the help of method AddChi2Term. The rank of an interpretation is defined as -chi^2.
 * 
 * Neutrino candidates are read from a dedicated reconstruction plugin with a default name
 * "NuReco". All candidates are considered for each way of jet assignment, and the one that gives
 * the smallest value of chi^2 is accepted.
 * 
 * The semileptonically decaying top quark is reconstructed using the leading charged lepton, which
 * is provided by a lepton trigger with a default name "Leptons".
 * 
 * If the event contains no charged leptons or no neutrino candidates have been reconstructed,
 * reconstruction is aborted. However, events are never rejected.
 */
class TTSemilepRecoChi2: public TTSemilepRecoBase
{
public:
    /// Supported types of summands in the chi^2
    enum class Expression
    {
        MassTopLep,  ///< Mass of semileptonically decaying top quark
        MassTopHad,  ///< Mass of hadronically decaying top quark
        MassWHad,    ///< Mass of the W boson from the hadronically decaying top quark
        PtTT         ///< Transverse momentum of the tt system
    };
    
private:
    /// An auxiliary structure to combine information about a single summand in the chi^2
    struct Chi2Term
    {
        /// Evaluates the chi^2 term
        double Eval(Lepton const &l, Candidate const &nu, Jet const &bTopLep,
          Jet const &bTopHad, Jet const &q1TopHad, Jet const &q2TopHad) const;
        
        /// Pointer to the function to evaluate the chi^2 term
        double (*expression)(Lepton const &l, Candidate const &nu, Jet const &bTopLep,
          Jet const &bTopHad, Jet const &q1TopHad, Jet const &q2TopHad);
        
        /// Mean value to be used in evaluation of chi^2
        double mean;
        
        /// Variance to be used in evaluation of chi^2
        double variance;
    };
    
public:
    /**
     * \brief Constructs a new plugin with the given name
     * 
     * User is encouraged to keep the default name.
     */
    TTSemilepRecoChi2(std::string name = "TTReco");
    
    /// Default move constructor
    TTSemilepRecoChi2(TTSemilepRecoChi2 &&) = default;
    
    /// Assignment operator is deleted
    TTSemilepRecoChi2 *operator=(TTSemilepRecoChi2 const &) = delete;
    
    /// Trivial destructor
    virtual ~TTSemilepRecoChi2() noexcept;
    
private:
    /**
     * \brief Copy constructor
     * 
     * Can only be called before processing of the first dataset has started.
     */
    TTSemilepRecoChi2(TTSemilepRecoChi2 const &src) noexcept;
    
public:
    /**
     * \brief Adds a new term to the chi^2 figure of merit
     * 
     * The functional form of the constraint is chosen from a predefined list by argument
     * expression, and remaining arguments provide the mean and the variance to calculate the
     * chi^2 term.
     */
    void AddChi2Term(Expression expression, double mean, double variance);
    
    /**
     * \brief Saves pointers to reader plugins
     * 
     * Reimplemented from TTSemilepRecoBase.
     */
    void BeginRun(Dataset const &) override;
    
    /**
     * \brief Creates a newly configured clone
     * 
     * Implemented from Plugin.
     */
    virtual Plugin *Clone() const override;
    
    /**
     * \brief Returns charged lepton from the t->blv decay
     * 
     * Returns the leading tight lepton. Throws an exception if the event contains no leptons.
     * 
     * Implemeted from TTSemilepRecoBase.
     */
    virtual Lepton const &GetLepton() const override;
    
    /**
     * \brief Returns reconstructed neutrino from the t->blv decay
     * 
     * This is one of neutrinos reconstructed by the dedicated plugin. Throws an exception if
     * the collection of neutrinos is empty or the jet assignment has been aborted.
     * 
     * Implemeted from TTSemilepRecoBase.
     */
    virtual Candidate const &GetNeutrino() const override;
    
private:
    /**
     * \brief Computes rank of the given event interpretation
     * 
     * The rank is defined as -chi^2. All neutrino solutions are considered, and the highest rank
     * (corresponding to smallest chi^2) is returned.
     * 
     * Implemented from TTSemilepRecoBase.
     */
    virtual double ComputeRank(Jet const &bTopLep, Jet const &bTopHad, Jet const &q1TopHad,
      Jet const &q2TopHad) override;
    
    /**
     * \brief Performs reconstruction of the current event
     * 
     * Calls PerformJetAssignment from the base class.
     * 
     * Reimplemented from TTSemilepRecoBase.
     */
    virtual bool ProcessEvent() override;
    
private:
    /// Name of the plugin that produces leptons
    std::string leptonPluginName;
    
    /// Non-owning pointer to the plugin that produces leptons
    LeptonReader const *leptonPlugin;
    
    /// Name of the plugin that reconstructs neutrinos
    std::string nuRecoPluginName;
    
    /// Non-owning pointer to the plugin that reconstructs neutrinos
    NuRecoBase const *nuRecoPlugin;
    
    /**
     * \brief Neutrino solution that is used in the best interpretation found so far
     * 
     * The pointer is non-owning and refers to the collection returned by the plugin for neutrino
     * reconstruction.
     */
    Candidate const *bestNu;
    
    /// Summands of the chi^2 figure of merit
    std::vector<Chi2Term> chi2Terms;
    
    /**
     * \brief Minimal value of chi^2 found in the current event so far
     * 
     * Needed to find the best neutrino solution. Reset to +infinity at the start of reconstruction
     * of each event.
     */
    double minChi2;
};
