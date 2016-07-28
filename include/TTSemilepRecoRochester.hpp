#pragma once

#include <TTSemilepRecoBase.hpp>

#include <mensura/core/PhysicsObjects.hpp>

#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>

#include <memory>


class LeptonReader;


/**
 * \class TTSemilepRecoRochester
 * \brief Reconstructs semileptonic tt decays following the Rochester algorithm
 * 
 * Reconstruction of the neutrino is delegated to class NuRecoRochester, which exploits constraints
 * from masses of the top quark and the W boson. Each event interpretation is ranked using a
 * (logarithm of the) product of two likelihood functions. One describes distribution of the
 * figure of merit used in the neutrino reconstruction, the other one represents the joint
 * distribution of the two reconstructed masses.
 * 
 * All possible jet assignments are considered, without any restriction due to b-tagging.
 * 
 * A version of this algorithm was used in TOP-16-008 (AN-16-020).
 */
class TTSemilepRecoRochester: public TTSemilepRecoBase
{
public:
    /**
     * \brief Creates reconstruction plugin with the given name
     * 
     * User is encouraged to keep the default name. In order to complete configuration of the
     * plugin, likelihood functions for the reconstruction must be provided using method
     * SetLikelihood.
     */
    TTSemilepRecoRochester(std::string name = "TTReco");
    
private:
    /// Copy constructor
    TTSemilepRecoRochester(TTSemilepRecoRochester const &src);
    
public:
    /**
     * \brief Saves pointers to reader plugins
     * 
     * Reimplemented from TTSemilepRecoBase.
     */
    virtual void BeginRun(Dataset const &dataset) override;
    
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
     * Throws an exception if the event reconstruction has failed for the current event.
     * 
     * Implemeted from TTSemilepRecoBase.
     */
    virtual Candidate const &GetNeutrino() const override;
    
    /**
     * \brief Provides likelihood function for reconstruction
     * 
     * The path to a ROOT file is resolved using FileInPath. Throws exceptions if the file is not
     * found or it does not contain one of the required histograms. The histograms do not need to
     * be normalized.
     */
    void SetLikelihood(std::string const &path,
      std::string const histNeutrinoName = "nusolver_chi2_right",
      std::string const histMassName = "mWhad_vs_mtophad_right");
    
private:
    /**
     * \brief Computes rank of the given event interpretation
     * 
     * The rank is defined as the logarithm of the likelihood computed with the figure of merit for
     * the neutrino reconstruction and masses of hadronically decaying t quark and W boson. A
     * simple caching is implemented for the neutrino reconstruction.
     * 
     * Implemented from TTSemilepRecoBase.
     */
    virtual double ComputeRank(Jet const &bTopLep, Jet const &bTopHad, Jet const &q1TopHad,
      Jet const &q2TopHad) override;
    
    /**
     * \brief Performs reconstruction of the current event
     * 
     * Calls PerformJetAssignment from the base class. Identifies reason for reconstruction
     * failure.
     * 
     * Reimplemented from TTSemilepRecoBase.
     */
    virtual bool ProcessEvent() override;
    
private:
    /// Name of the plugin that produces leptons
    std::string leptonPluginName;
    
    /// Non-owning pointer to the plugin that produces leptons
    LeptonReader const *leptonPlugin;
    
    /**
     * \brief Lepton selected for reconstruction of the tt system
     * 
     * The pointer is null if the current event contains no leptons.
     */
    Lepton const *lepton;
    
    /// Short-cut to access MET in the current event
    Candidate const *met;
    
    /**
     * \brief Histogram describing likelihood of neutrino solutions
     * 
     * The histograms is shared among all clones of this.
     */
    std::shared_ptr<TH1> likelihoodNeutrino;
    
    /**
     * \brief Histogram describing likelihood of reconstructed masses
     * 
     * Arguments are masses of hadronically decaying W boson and hadronically decaying top quark.
     * The histograms is shared among all clones of this.
     */
    std::shared_ptr<TH2> likelihoodMass;
    
    /// Current best neutrino candidate
    Candidate neutrino;
    
    /**
     * \brief Pointer to b-quark jet from t -> blv in the last considered interpretation
     * 
     * It allows to implement caching of neutrino reconstruction. The pointer must be reset to null
     * at the start of processing of each new event.
     */
    Jet const *cachedBTopLep;
    
    /// Cached four-momentum of reconstructed neutrino
    TLorentzVector cachedP4Nu;
    
    /// Cached log-likelihood corresponding to the neutrino
    double cachedLogLikelihoodNu;
    
    /// Flags that help to deduce reason of failed reconstruction
    bool neutrinoReconstructed, neutrinoLikelihoodInRange, massLikelihoodInRange;
};
