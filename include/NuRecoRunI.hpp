#pragma once

#include <NuRecoBase.hpp>


/**
 * \class NuRecoRunI
 * \brief A plugin to reconstruct neutrino in W->lv events
 * 
 * This plugin reconstructs neutrino from the leading lepton and MET. It exploits the W mass
 * constraint solve for the longitudinal component of neutrino momentum. If two solutions are
 * found, both corresponding neutrino candidates are built. When underlying quadratic equation has
 * no real-valued solutions, the value of MET is modified, while keeping its direction in the
 * transverse plane, until the discriminant turns zero. The neutrino candidate is built using the
 * resulting z component of the momentum and the modified MET.
 * 
 * This method was used in LHC Run 1 in AN-12-488 and AN-14-097. It is different from what is
 * usually used in analyses targeting single top quark production.
 */
class NuRecoRunI: public NuRecoBase
{
public:
    /**
     * \brief Constructs a new plugin with the given name
     * 
     * User is encouraged to keep the default name.
     */
    NuRecoRunI(std::string name = "NuReco");
    
    /// Default copy constructor
    NuRecoRunI(NuRecoRunI const &) = default;
    
    /// Default move constructor
    NuRecoRunI(NuRecoRunI &&) = default;
    
    /// Trivial destructor
    virtual ~NuRecoRunI() noexcept;
    
public:
    /**
     * \brief Creates a newly configured clone
     * 
     * Implemented from Plugin.
     */
    virtual Plugin *Clone() const override;
    
private:
    /**
     * \brief Reconstructs neutrino in the current event
     * 
     * Always returns true even if reconstruction fails.
     * Implemented from Plugin.
     */
    virtual bool ProcessEvent() override;
};
