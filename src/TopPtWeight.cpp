#include <TopPtWeight.hpp>

#include <mensura/core/GenParticleReader.hpp>
#include <mensura/core/Dataset.hpp>

#include <TLorentzVector.h>

#include <cmath>
#include <sstream>
#include <stdexcept>


TopPtWeight::TopPtWeight(std::string const name /*= "TopPtWeight"*/):
    EventWeightPlugin(name),
    genParticleReaderName("GenParticles"), genParticleReader(nullptr),
    datasetMasks({std::regex(".*")}), processCurDataset(false),
    nominalParams{6.15024e-02, -5.17833e-04},
    paramsVar1{0.03243, -1.404e-4}, paramsVar2{-4.353e-07, -1.005e-4},
    meanWeights{0.9985, 1.0142, 0.9832, 0.9865, 1.0107}
{}


void TopPtWeight::BeginRun(Dataset const &dataset)
{
    std::string const &datasetID = dataset.GetSourceDatasetID();
    
    
    // Check if weights need to be computed for the current dataset
    processCurDataset = false;
    
    for (auto const &mask: datasetMasks)
    {
        if (std::regex_match(datasetID, mask))
        {
            processCurDataset = true;
            break;
        }
    }
    
    
    if (processCurDataset)
    {
        // Save pointer to plugin that provides generator-level weights
        genParticleReader =
          dynamic_cast<GenParticleReader const *>(GetDependencyPlugin(genParticleReaderName));
        
        weights.resize(5);
        weights[0] = 1.;
    }
    else
    {
        weights.resize(1);
        weights[0] = 1.;
    }
}


Plugin *TopPtWeight::Clone() const
{
    return new TopPtWeight(*this);
}


void TopPtWeight::SelectDatasets(std::initializer_list<std::string> const &masks)
{
    datasetMasks.clear();
    
    for (auto const &mask: masks)
        datasetMasks.emplace_back(mask);
}


double TopPtWeight::ComputeTopPtWeight(double pt1, double pt2, int var1, int var2) const
{
    double const p0 = nominalParams[0] + var1 * paramsVar1[0] + var2 * paramsVar2[0];
    double const p1 = nominalParams[1] + var1 * paramsVar1[1] + var2 * paramsVar2[1];
    double rawWeight = std::sqrt(std::exp(p0 + p1 * pt1) * std::exp(p0 + p1 * pt2));
    
    return rawWeight;
}


bool TopPtWeight::ProcessEvent()
{
    // Do nothing if no reweighting is needed for the current dataset
    if (not processCurDataset)
        return true;
    
    
    // Find top quarks. Compute their pt from sum of momenta of their daughters to get pt of the
    //last top quarks in the event history.
    unsigned nTopFound = 0;
    std::array<TLorentzVector, 2> topP4;
    
    for (auto const &p: genParticleReader->GetParticles())
    {
        if (std::abs(p.GetPdgId()) != 6)
            continue;
        
        
        ++nTopFound;
        
        if (nTopFound == 3)
        {
            std::ostringstream message;
            message << "TopPtWeight[\"" << GetName() << "\"]::ProcessEvent: Found more than "
              "two top quarks in an event.";
            throw std::runtime_error(message.str());
        }
        
        
        unsigned nDaughters = 0;
        
        for (auto const &d: p.GetDaughters())
        {
            if (std::abs(d->GetPdgId()) != 24)
            {
                topP4[nTopFound - 1] += d->P4();
                ++nDaughters;
            }
            else
            {
                for (auto const &dW: d->GetDaughters())  // There are no chains like W->W->...
                {
                    topP4[nTopFound - 1] += dW->P4();
                    ++nDaughters;
                }
            }
        }
        
        if (nDaughters != 3)
        {
            std::ostringstream message;
            message << "TopPtWeight[\"" << GetName() << "\"]::ProcessEvent: Found a top quark "
              "with " << nDaughters << " daughters.";
            throw std::runtime_error(message.str());
        }
    }
    
    if (nTopFound < 2)
    {
        std::ostringstream message;
        message << "TopPtWeight[\"" << GetName() << "\"]::ProcessEvent: Found " << nTopFound <<
          " < 2 top quarks in an event.";
        throw std::runtime_error(message.str());
    }
    
    
    // Compute event weights taking into account normalization by mean weights
    double const pt1 = topP4[0].Pt(), pt2 = topP4[1].Pt();
    
    weights.at(0) = ComputeTopPtWeight(pt1, pt2, 0, 0) / meanWeights[0];
    weights.at(1) = ComputeTopPtWeight(pt1, pt2, 1, 0) / meanWeights[1];
    weights.at(2) = ComputeTopPtWeight(pt1, pt2, -1, 0) / meanWeights[2];
    weights.at(3) = ComputeTopPtWeight(pt1, pt2, 0, 1) / meanWeights[3];
    weights.at(4) = ComputeTopPtWeight(pt1, pt2, 0, -1) / meanWeights[4];
    
    
    /// The plugin does not perform any filtering
    return true;
}
