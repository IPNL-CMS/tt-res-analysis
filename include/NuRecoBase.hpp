#pragma once

#include <PECFwk/core/AnalysisPlugin.hpp>

#include <PECFwk/core/PhysicsObjects.hpp>

#include <string>
#include <vector>


class LeptonReader;
class JetMETReader;


/**
 * \class NuRecoBase
 * \brief An abstract base class to perform reconstruction of full neutrino momentum
 * 
 * This is a base class for a plugin to reconstruct momentum of a neutrino in events with a
 * leptonically decaying W boson. For the sake of convenience, it sets up pointers to relevant
 * reader plugins (with default names "Leptons" and "JetMET"). User must implement method
 * ProcessEvent performing actual reconstruction of the neutrino and save the candidates in the
 * dedicated collection.
 */
class NuRecoBase: public AnalysisPlugin
{
public:
    /**
     * \brief Constructs a reconstruction plugin with the given name
     * 
     * User is encouraged to keep the default name.
     */
    NuRecoBase(std::string name = "NuReco");
    
    /// Copy constructor
    NuRecoBase(NuRecoBase const &src) noexcept;
    
    /// Default move constructor
    NuRecoBase(NuRecoBase &&) = default;
    
    /// Trivial destructor
    virtual ~NuRecoBase() noexcept;
    
public:
    /**
     * \brief Saves pointers to reader plugins
     * 
     * Reimplemented from Plugin.
     */
    void BeginRun(Dataset const &) override;
    
    /**
     * \brief Returns neutrinos reconstructed in the current event
     * 
     * The returned collection might be empty in case no neutrinos are reconstructed. It might
     * also contain more entries than expected from physics in case several alternative candidates
     * are reconstructed.
     */
    std::vector<Candidate> const &GetNeutrinos() const;
    
protected:
    /// Name of the plugin that produces leptons
    std::string leptonPluginName;
    
    /// Non-owning pointer to the plugin that produces leptons
    LeptonReader const *leptonPlugin;
    
    /// Name of the plugin that produces MET
    std::string jetmetPluginName;
    
    /// Non-owning pointer to the plugin that produces MET
    JetMETReader const *jetmetPlugin;
    
    /// Collection of neutrinos reconstructed in the current event
    std::vector<Candidate> neutrinos;
};
