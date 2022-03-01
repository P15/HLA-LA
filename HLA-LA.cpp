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

    // what follows is the new multi-step mapping
    if(arguments.at("action") == "multiHLA")
    {
        // ../bin/HLA-PRG-LA --action HLA --sampleID NA12878 --BAM /gpfs1/well/gsk_hla/temp_mapping_2/NA12878_PRG/merged.bam --outputDirectory /gpfs1/well/gsk_hla/HLA-PRG-LA/working/NA12878_PLATINUM --trueHLA /Net/birch/data/dilthey/1000GHLA/G1000_GSK_combined.txt.manuallyAmended
        // ../bin/HLA-PRG-LA --action HLA --sampleID NA12892_REDUCED --BAM /gpfs1/well/gsk_hla/PRG_Remapping/BAMs/PLATINUM_NA12892/merged.bam --outputDirectory /gpfs1/well/gsk_hla/HLA-PRG-LA/working/NA12892_PLATINUM_RED --trueHLA /Net/birch/data/dilthey/1000GHLA/G1000_GSK_combined.txt.manuallyAmended
        assert(arguments.count("sampleID"));
        assert(arguments.count("BAM") || (arguments.count("FASTQ1") && arguments.count("FASTQ2")));
        assert(arguments.count("outputDirectory"));
        assert(arguments.count("PRG_graph_dir"));
        std::string sampleID = arguments.at("sampleID");
        std::string BAM = arguments.at("BAM");
        std::string outputDirectory = arguments.at("outputDirectory");
        std::string PRG_graph_dir = arguments.at("PRG_graph_dir");
        if(arguments.count("FASTQ1"))
        {
            assert(arguments.count("mapAgainstCompleteGenome"));
            assert(! arguments.count("BAM"));
        }

        std::string file_true_HLA_types = (arguments.count("trueHLA") ? arguments.at("trueHLA") : "");
        if(! Utilities::directoryExists(outputDirectory))
        {
            Utilities::makeDir(outputDirectory);
        }
        //
        std::map<std::string, std::map<std::string, std::pair<std::set<std::string>, std::set<std::string>>>> inferredHLA;
        mapper::processBAM BAMprocessor (PRG_graph_dir);
        std::pair<double, double> IS_estimate;
        if(arguments.count("BAM"))
        {
            assert(Utilities::fileExists(BAM));
            IS_estimate = BAMprocessor.estimateInsertSize(BAM, true);
            bool alwaysStartFromScratch = true;
            if(alwaysStartFromScratch)
            {
                Utilities::make_or_clearDirectory(outputDirectory);
            }
            else
            {
                if(!Utilities::directoryExists(outputDirectory))
                {
                    Utilities::makeDir(outputDirectory);
                }
            }
            std::vector<std::string> sampledReferenceGenomes = Utilities::getAllLines(PRG_graph_dir + "/sampledReferenceGenomes.txt");
            // extract reads
            std::string fastq_extracted_forRemapping_1;
            std::string fastq_extracted_forRemapping_2;
            std::string dir_extractedReads_for_remapping = outputDirectory + "/extractedReads_forRemapping";
            if(alwaysStartFromScratch || (!Utilities::directoryExists(dir_extractedReads_for_remapping)))
            {
                Utilities::make_or_clearDirectory(dir_extractedReads_for_remapping);
                std::string PRG_sequences_file = PRG_graph_dir + "/sequences.txt";
                std::map<std::string, std::pair<int, int>> regions_for_extraction;
                {
                    // read regions for extraction
                    std::ifstream PRG_covered_regions_stream;
                    PRG_covered_regions_stream.open(PRG_sequences_file.c_str());
                    assert(PRG_covered_regions_stream.is_open());
                    assert(PRG_covered_regions_stream.good());
                    std::string line;
                    std::getline(PRG_covered_regions_stream, line);
                    std::string headerLine = line;
                    Utilities::eraseNL(headerLine);
                    std::vector<std::string> headerFields = Utilities::split(headerLine, "\t");
                    assert(headerFields.at(0) == "SequenceID");
                    assert(headerFields.at(1) == "Name");
                    assert(headerFields.at(2) == "FASTAID");
                    assert(headerFields.at(3) == "Chr");
                    assert(headerFields.at(4) == "Start_1based");
                    assert(headerFields.at(5) == "Stop_1based");
                    while(PRG_covered_regions_stream.good())
                    {
                        std::getline(PRG_covered_regions_stream, line);
                        Utilities::eraseNL(line);
                        if(line.length() == 0)
                        {
                            continue;
                        }
                        std::vector<std::string> fields = Utilities::split(line, "\t");
                        assert(fields.size() == headerFields.size());
                        std::string BAMid;
                        int startIndex_0based;
                        int stopIndex_0based;
                        std::string FastaID = fields.at(2);
                        std::string Chr = fields.at(3);
                        std::string start_str = fields.at(4);
                        std::string stop_str = fields.at(5);
                        if(Chr != "")
                        {
                            BAMid = Chr;
                            assert(start_str.length());
                            assert(stop_str.length());
                            startIndex_0based = Utilities::StrtoI(start_str);
                            stopIndex_0based = Utilities::StrtoI(stop_str);
                            assert(startIndex_0based >= 0);
                            assert(stopIndex_0based >= 0);
                            assert(startIndex_0based <= stopIndex_0based);
                        }
                        else
                        {
                            assert(FastaID.length());
                            BAMid = FastaID;
                            assert(start_str.length() == 0);
                            assert(stop_str.length() == 0);
                            startIndex_0based = -1;
                            stopIndex_0based = -1;
                        }
                        assert(BAMid.length());
                        assert(((startIndex_0based == -1) && (stopIndex_0based == -1)) || ((startIndex_0based != -1) && (stopIndex_0based != -1) && (startIndex_0based <= stopIndex_0based)));
                        regions_for_extraction[BAMid] = std::make_pair(startIndex_0based, stopIndex_0based);
                    }
                }
                linearALTs::linearALTs::extractReadsFromBAM(dir_extractedReads_for_remapping, BAM, regions_for_extraction);
                fastq_extracted_forRemapping_1 = dir_extractedReads_for_remapping + "/R_1.fq";
                fastq_extracted_forRemapping_2 = dir_extractedReads_for_remapping + "/R_2.fq";
                assert(Utilities::fileExists(fastq_extracted_forRemapping_1));
                assert(Utilities::fileExists(fastq_extracted_forRemapping_2));
            }
            else
            {
                fastq_extracted_forRemapping_1 = dir_extractedReads_for_remapping + "/R_1.fq";
                fastq_extracted_forRemapping_2 = dir_extractedReads_for_remapping + "/R_2.fq";
                assert(Utilities::fileExists(fastq_extracted_forRemapping_1));
                assert(Utilities::fileExists(fastq_extracted_forRemapping_2));
            }
            std::vector<std::string> remapped_BAMs;
            for(unsigned int genomeI = 0; genomeI < sampledReferenceGenomes.size(); genomeI++)
            {
                std::string sampledReferenceGenome = sampledReferenceGenomes.at(genomeI);
                std::string BAM_remapped = outputDirectory + "/sampled_" + Utilities::ItoStr(genomeI) + ".bam";
                std::cout << Utilities::timestamp() << "Carry out remapping step for sampled genome " << genomeI << " - input " << BAM << ", output " << BAM_remapped << "\n" << std::flush;
                if(alwaysStartFromScratch || (! Utilities::fileExists(BAM_remapped)))
                {
                    mapper::bwa::BWAmapper bwaMapper(pF);
                    bwaMapper.map(sampledReferenceGenome, fastq_extracted_forRemapping_1, fastq_extracted_forRemapping_2, BAM_remapped, true);
                }
                else
                {
                    assert(Utilities::fileExists(BAM_remapped));
                }
                std::cout << Utilities::timestamp() << "Remapping done.\n" << std::flush;
                assert(Utilities::fileExists(BAM_remapped+".bai"));
                remapped_BAMs.push_back(BAM_remapped);
            }
        }
        Graph* g = BAMprocessor.getGraph();
        hla::HLATyper HLAtyper(g, PRG_graph_dir, "");
        BAMprocessor.alignReadsMulti(remapped_BAMs, 0, IS_estimate.first, IS_estimate.second, outputDirectory, true, &HLAtyper);
        std::string expected_HLA_type_inference_output = outputDirectory + "/hla/R1_bestguess.txt";
        assert(Utilities::fileExists(expected_HLA_type_inference_output));
        hla::HLATyper::read_inferred_types(sampleID, inferredHLA, expected_HLA_type_inference_output);
        if(file_true_HLA_types.length())
        {
            std::map<std::string, std::map<std::string, std::pair<std::string, std::string>>> trueHLA;
            hla::HLATyper::read_true_types(trueHLA, file_true_HLA_types);
            hla::HLATyper::evaluate_HLA_types(trueHLA, inferredHLA);
        }
    }
    else if(arguments.at("action") == "HLA")
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
