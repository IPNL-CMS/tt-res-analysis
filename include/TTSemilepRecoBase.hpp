#pragma once

#include <mensura/core/AnalysisPlugin.hpp>

#include <mensura/core/PhysicsObjects.hpp>

#include <limits>
#include <string>
#include <vector>


class JetMETReader;


/**
 * \class TTSemilepRecoBase
 * \brief An abstract base class to perform jet assignment in a semileptonic ttbar event
 * 
 * This is a base class for a plugin to identify reconstructed jets corresponding to the four
 * quarks in the final state of tt -> blv bqq. The plugin implements looping over all possible ways
 * to assign four reconstructed jets to the four quarks, or interpretations. A derived class must
 * implement method to compute the rank of an interpretation. The plugin accepts interpretation
 * that gets the highest rank.
 * 
 * Reconstruction of the neutrino is delegated to the derived class. If multiple candidates can
 * be reconstructed in a single event, it must choose the most suitable one. It provides
 * reconstructed neutrino and also selected charged lepton by implementing pure virtual methods
 * GetLepton and GetNeutrino.
 * 
 * No reconstruction is performed if an event contains less than four jets satisfying the
 * selection. However, this plugin never rejects events.
 * 
 * This plugin relies on a jet reader with the default name "JetMET".
 */
class TTSemilepRecoBase: public AnalysisPlugin
{
public:
    /// Jets to be identified in the final state of a ttbar system
    enum class DecayJet
    {
        bTopLep,   ///< Jet from semileptonically decaying top quark
        bTopHad,   ///< Jet from hadronization of b quark from hadronically decaying top quark
        q1TopHad,  ///< Leading light-flavour jet from hadronically decaying top quark
        q2TopHad   ///< Subleading light-flavour jet from hadronically decaying top quark
    };
    
public:
    /**
     * \brief Constructs a new plugin with the given name
     * 
     * User is encouraged to keep the default name.
     */
    TTSemilepRecoBase(std::string name = "TTReco");
    
    /// Default move constructor
    TTSemilepRecoBase(TTSemilepRecoBase &&) = default;
    
    /// Assignment operator is deleted
    TTSemilepRecoBase &operator=(TTSemilepRecoBase const &) = delete;
    
    /// Trivial destructor
    virtual ~TTSemilepRecoBase() noexcept;
    
protected:
    /**
     * \brief Copy constructor
     * 
     * Can only be called before processing of the first dataset has started.
     */
    TTSemilepRecoBase(TTSemilepRecoBase const &src) noexcept;
    
public:
    /**
     * \brief Saves pointer to jet reader
     * 
     * User may override this method, but in this case the overriding method must initialize the
     * pointer jetmetPlugin (possibly by calling method of the base class).
     * 
     * Reimplemented from Plugin.
     */
    void BeginRun(Dataset const &) override;
    
    /**
     * \brief Returns jet corresponding to the given quark in the final state tt -> blv bqq
     * 
     * The behaviour is undefined if reconstruction has failed.
     */
    Jet const &GetJet(DecayJet type) const;
    
    /// A pure virtual method to return charged lepton from the t->blv decay
    virtual Lepton const &GetLepton() const = 0;
    
    /// A pure virual method to return reconstructed neutrino from the t->blv decay
    virtual Candidate const &GetNeutrino() const = 0;
    
    /**
     * \brief Returns rank of the accepted interpretation of the current event
     * 
     * If reconstruction has failed, returns -infinity.
     */
    double GetRank() const;
    
    /**
     * \brief Reports if reconstruction of the current event has been successful
     * 
     * Returns false if reconstruction of the current event has been aborted. This happens if the
     * event contains less than four jets passing the selection or if reconstruction failure is
     * declared by a derived class.
     */
    bool GetRecoStatus() const;
    
    /// Compute and return four-momentum of reconstructed leptonically decaying top quark
    TLorentzVector GetTopLepP4() const;
    
    /// Compute and return four-momentum of reconstructed hadronically decaying top quark
    TLorentzVector GetTopHadP4() const;
    
    /**
     * \brief Sets jet selection
     * 
     * Only jets satisfying this selection are tried as decay products of the top quarks.
     */
    void SetJetSelection(double minPt, double maxAbsEta = std::numeric_limits<double>::infinity());
    
protected:
    /**
     * \brief Performs jet assignment in the current event
     * 
     * Considers all possible ways to choose four reconsted jets and match them to decay products
     * of a pair of top quarks. Jets must satisfy the selection on pt and |eta|. For each possible
     * interpretation its rank is computed with method ComputeRank, and the interpretation with the
     * highest rank is accepted.
     * 
     * If the event contain less than for jets satisfying the selection, reconstruction is not
     * performed. In this case the highest rank is set to -infinity.
     * 
     * Provided collection of jets must exist and not be reallocated until the end of processing of
     * the current event. This plugin will save pointers to identified jets.
     */
    void PerformJetAssignment(std::vector<Jet> const &jets);
    
    /// Declares that reconstuction of the current event has been aborted
    void SetRecoFailure();
    
private:
    /// A pure virtual method to calculate rank of a given interpretation of the current event
    virtual double ComputeRank(Jet const &bTopLep, Jet const &bTopHad, Jet const &q1TopHad,
      Jet const &q2TopHad) = 0;
    
    /**
     * \brief Performs reconstruction of the current event by calling PerformJetAssignment
     * 
     * Implemented from Plugin.
     */
    virtual bool ProcessEvent() override;
    
protected:
    /// Name of the plugin that produces MET
    std::string jetmetPluginName;
    
    /// Non-owning pointer to the plugin that produces MET
    JetMETReader const *jetmetPlugin;
    
private:
    /// Selection on jet transverse momentum
    double minPt;
    
    /// Selection on absolute value of jet pseudorapidity
    double maxAbsEta;
    
    /**
     * \brief Indices of jets that pass the selection on pt and |eta|
     * 
     * This vector is only used in the method PerformJetAssignment but placed here to avoid
     * reallocation of memory for each event.
     */
    std::vector<unsigned> selectedJetIndices;
    
    /// Flag showing whether reconstruction of the current event has been successful
    bool recoSuccess;
    
    /**
     * \brief Rank of the best interpretation constructed so far
     * 
     * When reconstruction of a new event starts, this variable is reset to -infinity. After all
     * interpretations in an event have been considered, it contains the rank of the best
     * interpretation.
     */
    double highestRank;
    
    /**
     * \brief Non-owning pointers to jets identified as decay products of the top quarks
     * 
     * Pointers refer to jets in the collection provided by the jet reader.
     */
    Jet const *bTopLep, *bTopHad, *q1TopHad, *q2TopHad;
};
