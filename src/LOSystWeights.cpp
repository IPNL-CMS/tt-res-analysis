#include <LOSystWeights.hpp>

#include <mensura/core/Dataset.hpp>
#include <mensura/core/FileInPath.hpp>
#include <mensura/core/GeneratorReader.hpp>

#include <LHAPDF/LHAPDF.h>

#include <cmath>


LOSystWeights::LOSystWeights(std::string const &name, unsigned nQCDVert_,
  std::string const &pdfSetName):
    EventWeightPlugin(name),
    generatorReaderName("Generator"), generatorReader(nullptr),
    datasetMasks({std::regex(".*")}), processCurDataset(false),
    scaleVarFactor(2.),
    nQCDVert(nQCDVert_),
    pdfSet(LHAPDF::mkPDF(pdfSetName, 0))
{}


LOSystWeights::LOSystWeights(unsigned nQCDVert, std::string const &pdfSetName):
    LOSystWeights("LOSystWeights", nQCDVert, pdfSetName)
{}


LOSystWeights::LOSystWeights(LOSystWeights const &src):
    EventWeightPlugin(src),
    generatorReaderName(src.generatorReaderName), generatorReader(nullptr),
    datasetMasks(src.datasetMasks), processCurDataset(src.processCurDataset),
    scaleVarFactor(src.scaleVarFactor),
    nQCDVert(src.nQCDVert),
    pdfSet(src.pdfSet)
{}


void LOSystWeights::BeginRun(Dataset const &dataset)
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
        generatorReader =
          dynamic_cast<GeneratorReader const *>(GetDependencyPlugin(generatorReaderName));
        
        weights.resize(5);
        weights[0] = 1.;
    }
    else
    {
        weights.resize(1);
        weights[0] = 1.;
    }
}


Plugin *LOSystWeights::Clone() const
{
    return new LOSystWeights(*this);
}


void LOSystWeights::SelectDatasets(std::initializer_list<std::string> const &masks)
{
    datasetMasks.clear();
    
    for (auto const &mask: masks)
        datasetMasks.emplace_back(mask);
}


double LOSystWeights::AlphaS(double scale)
{
    double const alpha0 = 0.1184;
    double const mZ = 91.1876;
    unsigned const nf = 4;
    
    double const b0 = (33 - 2 * nf) / (12 * M_PI);
    return alpha0 / (1 + alpha0 * b0 * 2 * std::log(scale / mZ));
}


bool LOSystWeights::ProcessEvent()
{
    // Do nothing if no weights need to be computed for the current dataset
    if (not processCurDataset)
        return true;
    
    
    double const scale = generatorReader->GetScale();
    unsigned iVar;
    
    
    // Variation of renormalization scale
    iVar = 0;
    double const alphaSNominal = AlphaS(scale);
    weights.at(1 + 2 * iVar) = std::pow(AlphaS(scale * scaleVarFactor) / alphaSNominal, 2);
    weights.at(2 + 2 * iVar) = std::pow(AlphaS(scale / scaleVarFactor) / alphaSNominal, 2);
    
    
    // Variation of factorization scale
    iVar = 1;
    int const id1 = generatorReader->GetPdfPart().first;
    int const id2 = generatorReader->GetPdfPart().second;
    double const &x1 = generatorReader->GetPdfX().first, x2 = generatorReader->GetPdfX().second;
    
    double const pdfNominal = pdfSet->xfxQ(id1, x1, scale) * pdfSet->xfxQ(id2, x2, scale);
    weights.at(1 + 2 * iVar) = pdfSet->xfxQ(id1, x1, scale * scaleVarFactor) * \
      pdfSet->xfxQ(id2, x2, scale * scaleVarFactor) / pdfNominal;
    weights.at(2 + 2 * iVar) = pdfSet->xfxQ(id1, x1, scale / scaleVarFactor) * \
      pdfSet->xfxQ(id2, x2, scale / scaleVarFactor) / pdfNominal;
    
    
    return true;
}
