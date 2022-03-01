//============================================================================
// Name    : HLA-PRG-LA.cpp
// Author    :
// Version   :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>

#include <map>
#include <assert.h>
#include <string>
#include <vector>
#include <exception>
#include <stdexcept>
#include <chrono>
#include <cstdio>

#include "mapper/processBAM.h"
#include "mapper/reads/PRGContigBAMAlignment.h"
#include "mapper/reads/verboseSeedChain.h"
#include "mapper/aligner/extensionAligner.h"
#include "mapper/bwa/BWAmapper.h"
#include "mapper/bowtie2/Bowtie2mapper.h"

#include "simulator/trueReadLevels.h"

#include "Graph/Graph.h"
#include "Graph/graphSimulator/simpleGraphSimulator.h"
#include "simulator/simulator.h"
#include "Graph/GraphAndEdgeIndex.h"

#include "Utilities.h"
#include "pathFinder.h"

#include "hla/HLATyper.h"

#include "linearALTs/linearALTs.h"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/set.hpp>

int main(int argc, char *argv[]) {

    std::vector<std::string> ARG (argv + 1, argv + argc + !argc);
    std::map<std::string, std::string> arguments;

    /*
       arguments["action"] = "PRGmapping";
       arguments["action"] = "testChainExtension";
       arguments["action"] = "testAlignments2Chains";
       arguments["action"] = "testPRGMapping";
       arguments["action"] = "testPRGMappingUnpaired";
     */

    // arguments["action"] = "testRealBAM";

    // arguments["action"] = "prepareGraph";
    // arguments["action"] = "TestHLATyping";

    // arguments["action"] = "testCheckPresence";
    // arguments["action"] = "checkKIRgraph";

    // arguments["action"] = "KIR";

    // arguments["PRG_graph_dir"] = "/gpfs1/well/gsk_hla/HLA-PRG-LA/graphs/PRG_MHC_GRCh38_withIMGT";

    for(unsigned int i = 0; i < ARG.size(); i++)
    {
        if((ARG.at(i).length() > 2) && (ARG.at(i).substr(0, 2) == "--"))
        {
            std::string argname = ARG.at(i).substr(2);
            std::string argvalue = ARG.at(i+1);
            arguments[argname] = argvalue;
        }
    }


    /*

       std::map<std::string, std::map<std::string, std::pair<std::string, std::string>>> trueHLA_test;
       std::map<std::string, std::map<std::string, std::pair<std::set<std::string>, std::set<std::string>>>> inferredHLA_test;

       hla::HLATyper::read_inferred_types("NA12892", inferredHLA_test, "/Net/birch/data/dilthey/MHC-PRG/tmp/hla/C_Platinum_NA12892/R1_bestguess.txt");
       hla::HLATyper::read_true_types(trueHLA_test, "/Net/birch/data/dilthey/1000GHLA/G1000_GSK_combined.txt.manuallyAmended");
       hla::HLATyper::evaluate_HLA_types(trueHLA_test, inferredHLA_test);
       assert(2 == 5);
     */

    // some overlap tests
    assert(!Utilities::intervalsOverlap(1, 10, 11, 20));
    assert(!Utilities::intervalsOverlap(5, 11, 1, 4));
    assert(Utilities::intervalsOverlap(5, 11, 8, 11));
    assert(Utilities::intervalsOverlap(8, 11, 1, 9));
    assert(Utilities::intervalsOverlap(8, 11, 1, 9));
    assert(Utilities::intervalsOverlap(8, 11, 9, 10));
    assert(Utilities::intervalsOverlap(9, 10, 8, 11));
    assert(Utilities::intervalsOverlap(1, 10, 2, 3));
    assert(Utilities::intervalsOverlap(2, 3, 1, 10));

    if(arguments.count("action") == 0)
    {
        std::cerr << "\n\nMissing --action parameter. Please don't try calling me directly; use HLA-LA.pl instead (see documentation on GitHub).\n" << std::endl;
        throw std::runtime_error("Missing arguments -- see above.");
    }

    std::set<std::string> noBinariesRequired = {"prepareGraph", "testBinary"};

    if((noBinariesRequired.count(arguments.at("action"))  == 0) && (!arguments.count("bwa_bin")))
    {
        throw std::runtime_error("Please specify arguments --bwa_bin");
    }
    if((noBinariesRequired.count(arguments.at("action"))  == 0) && (!arguments.count("samtools_bin")))
    {
        throw std::runtime_error("Please specify arguments --samtools_bin");
    }

    if(!std::system(NULL))
    {
        std::cerr << "\n\nMissing shell - std::system(NULL) has returned a 0 value.\n" << std::endl;
        throw std::runtime_error("Missing shell");
    }

    pathFinder pF(arguments);
    assert(arguments.count("action"));

    if(arguments.at("action") == "HLA")
    {

        // ../bin/HLA-PRG-LA --action HLA --sampleID NA12878 --BAM /gpfs1/well/gsk_hla/temp_mapping_2/NA12878_PRG/merged.bam --outputDirectory /gpfs1/well/gsk_hla/HLA-PRG-LA/working/NA12878_PLATINUM --trueHLA /Net/birch/data/dilthey/1000GHLA/G1000_GSK_combined.txt.manuallyAmended
        // ../bin/HLA-PRG-LA --action HLA --sampleID NA12892_REDUCED --BAM /gpfs1/well/gsk_hla/PRG_Remapping/BAMs/PLATINUM_NA12892/merged.bam --outputDirectory /gpfs1/well/gsk_hla/HLA-PRG-LA/working/NA12892_PLATINUM_RED --trueHLA /Net/birch/data/dilthey/1000GHLA/G1000_GSK_combined.txt.manuallyAmended

        unsigned int maxThreads = 1;

        std::cout << "Running script...\n";

        assert(arguments.count("sampleID"));
        assert(arguments.count("BAM") || (arguments.count("FASTQ1") && arguments.count("FASTQ2")));
        assert(arguments.count("outputDirectory"));
        assert(arguments.count("PRG_graph_dir"));

        std::string sampleID = arguments.at("sampleID");
        std::string outputDirectory = arguments.at("outputDirectory");
        std::string PRG_graph_dir = arguments.at("PRG_graph_dir");
        if(arguments.count("FASTQ1"))
        {
            assert(arguments.count("mapAgainstCompleteGenome"));
            assert(!arguments.count("BAM"));
        }
        if(arguments.count("maxThreads"))
        {
            maxThreads = Utilities::StrtoI(arguments.at("maxThreads"));
            std::cout << "Set maxThreads to " << maxThreads << "\n" << std::flush;
        }

        std::string file_true_HLA_types = (arguments.count("trueHLA") ? arguments.at("trueHLA") : "");

        if(!Utilities::directoryExists(outputDirectory))
        {
            Utilities::makeDir(outputDirectory);
        }
        // t

        std::map<std::string, std::map<std::string, std::pair<std::set<std::string>, std::set<std::string> > > > inferredHLA;


        std::cout << "BAMprocessor...\n";
        mapper::processBAM BAMprocessor (PRG_graph_dir, maxThreads);
        std::cout << "BAMprocessor complete...\n";
        std::pair<double, double> IS_estimate;
        std::string BAM_remapped = outputDirectory + "/remapped_with_a.bam";
        std::string PRGonlyReferenceGenomePath = PRG_graph_dir + "/mapping_PRGonly/referenceGenome.fa";
        std::string extendedReferenceGenomePath;
        if(Utilities::fileExists(PRG_graph_dir + "/extendedReferenceGenomePath.txt"))
        {
            extendedReferenceGenomePath  = Utilities::getFirstLine(PRG_graph_dir + "/extendedReferenceGenomePath.txt");
        }
        else
        {
            extendedReferenceGenomePath  = PRG_graph_dir + "/extendedReferenceGenome/extendedReferenceGenome.fa";
            assert(Utilities::fileExists(extendedReferenceGenomePath));
        }
        mapper::bwa::BWAmapper bwaMapper(pF, maxThreads);
        bool remap_with_a = true;
        if(arguments.count("remap_with_a"))
        {
            remap_with_a = Utilities::StrtoB(arguments.at("remap_with_a"));
        }

        bool remapped_against_extended_reference_genome;

        std::string longReads;

        std::cout << "Begin analysis...\n";
        assert(arguments.count("FASTQ1"));
        assert(arguments.count("FASTQ2"));
        assert(arguments.count("mapAgainstCompleteGenome"));
        assert(arguments.count("longReads"));

        longReads = arguments.at("longReads");
        assert((longReads == "0") || (longReads == "ont2d") || (longReads == "pacbio"));
        if(longReads == "0")
        {
            longReads = "";
        }

        if(longReads.length())
        {
            assert(arguments.count("FASTQU"));
        }

        bool mapAgainstCompleteGenome = Utilities::StrtoB(arguments.at("mapAgainstCompleteGenome"));

        std::string referenceGenomeForMapping = mapAgainstCompleteGenome ? extendedReferenceGenomePath : PRGonlyReferenceGenomePath;
        if(longReads.length() != 0)
        {
            bwaMapper.mapLong(referenceGenomeForMapping, arguments.at("FASTQU"), BAM_remapped, remap_with_a, longReads);
        }
        else
        {
            bwaMapper.map(referenceGenomeForMapping, arguments.at("FASTQ1"), arguments.at("FASTQ2"), BAM_remapped, remap_with_a);
        }

        std::cout << Utilities::timestamp() << "Remapping done.\n" << std::flush;
        assert(Utilities::fileExists(BAM_remapped));
        assert(Utilities::fileExists(BAM_remapped+".bai"));

        if(longReads.length() == 0)
        {
            IS_estimate = BAMprocessor.estimateInsertSize(BAM_remapped, mapAgainstCompleteGenome);
        }

        remapped_against_extended_reference_genome = mapAgainstCompleteGenome;

        std::cout << "Getting graph data...\n";
        Graph* g = BAMprocessor.getGraph();
        hla::HLATyper HLAtyper(g, PRG_graph_dir, "");

        // BAMprocessor.alignReads(BAM_remapped, 0, IS_estimate.first, IS_estimate.second, outputDirectory, false, &HLAtyper);
        BAMprocessor.alignReads_and_inferHLA(BAM_remapped, 0, IS_estimate.first, IS_estimate.second, outputDirectory, remapped_against_extended_reference_genome, &HLAtyper, 1, longReads);

        std::string expected_HLA_type_inference_output = outputDirectory + "/hla/R1_bestguess.txt";
        assert(Utilities::fileExists(expected_HLA_type_inference_output));
        hla::HLATyper::read_inferred_types(sampleID, inferredHLA, expected_HLA_type_inference_output);

        if(file_true_HLA_types.length())
        {
            std::map<std::string, std::map<std::string, std::pair<std::string, std::string> > > trueHLA;
            hla::HLATyper::read_true_types(trueHLA, file_true_HLA_types);
            hla::HLATyper::evaluate_HLA_types(trueHLA, inferredHLA);
        }
    }
    else
    {
        throw std::runtime_error("Invalid --action: " + arguments.at("action"));
    }

    return 0;
}

