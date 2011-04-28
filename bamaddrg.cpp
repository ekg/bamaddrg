#include <map>
#include <vector>
#include <string>
#include <getopt.h>
#include <iostream>
#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include "api/SamReadGroup.h"

using namespace BamTools;
using namespace std;

void printUsage(int argc, char** argv) {

    cerr << "usage: " << argv[0] << " [-b FILE [-s NAME [-r GROUP]]]" << endl
         << endl
         << "options:" << endl
         << "    -h, --help         this dialog" << endl
         << "    -b, --bam FILE     use this BAM as input" << endl
         << "    -s, --sample NAME  optionally apply this sample name to the preceeding BAM file" << endl
         << "    -r, --read-group GROUP  optionally apply this read group to the preceeding BAM file" << endl
         << endl
         << "Merges the alignments in the supplied BAM files, using the supplied sample names" << endl
         << "and read groups to specifically add read group (RG) tags to each alignment.  The" << endl
         << "output is uncompressed, and is suitable for input into downstream alignment systems" << endl
         << "which require RG tag information." << endl
         << endl
         << "Sample names and read groups may be specified by supplying a sample name or read group" << endl
         << "argument after each listed BAM file." << endl
         << endl
         << "When no sample names are supplied, the names of the BAM files are used as the sample" << endl
         << "names and read groups.  When no read groups are supplied, the sample names are used" << endl
         << "as read groups." << endl
         << endl
         << "author: Erik Garrison <erik.garrison@bc.edu>" << endl;

}

int main(int argc, char** argv) {

    vector<string> inputFilenames;
    vector<string> sampleNames;
    vector<string> readGroups;
    
    string currFileName;
    string currReadGroup;
    string currSampleName;

    // parse command-line options
    int c;

    while (true) {

        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"bam",  no_argument, 0, 'b'},
            {"read-group", no_argument, 0, 'r'},
            {"sample", no_argument, 0, 's'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hb:s:r:",
                         long_options, &option_index);

        if (c == -1)
            break;
 
        switch (c) {

            case '?':
                printUsage(argc, argv);
                return 0;
                break;

            case 'h':
                printUsage(argc, argv);
                return 0;
                break;

            case 'b':
                if (!currFileName.empty()) {
                    inputFilenames.push_back(currFileName);
                    if (currSampleName.empty()) {
                        currSampleName = currFileName;
                    }
                    sampleNames.push_back(currSampleName);
                    // use the sample name by default
                    readGroups.push_back(currReadGroup.empty() ? currSampleName : currReadGroup);
                    currFileName.clear();
                    currSampleName.clear();
                    currReadGroup.clear();
                }
                currFileName = optarg;
                break;

            case 's':
                currSampleName = optarg;
                break;

            case 'r':
                currReadGroup = optarg;
                break;

            default:
                return 1;
                break;
        }
    }

    // catch last iteration
    if (!currFileName.empty()) {
        inputFilenames.push_back(currFileName);
        if (currSampleName.empty()) {
            currSampleName = currFileName;
        }
        sampleNames.push_back(currSampleName);
        // use the sample name by default
        readGroups.push_back(currReadGroup.empty() ? currSampleName : currReadGroup);
        currFileName.clear();
        currSampleName.clear();
        currReadGroup.clear();
    }

    if (inputFilenames.empty()) {
        cerr << "no input files specified" << endl;
        return 1;
    }

    map<string, string> filenameToReadGroup;
    map<string, string> readGroupToSampleName;
    vector<SamReadGroup> newReadGroups;

    vector<string>::iterator b = inputFilenames.begin();
    vector<string>::iterator r = readGroups.begin();
    vector<string>::iterator s = sampleNames.begin();
    for (; b != inputFilenames.end(); ++b, ++r, ++s) {
        cerr << *b << " " << *s << " " << *r << endl;
        filenameToReadGroup[*b] = *r; 
        readGroupToSampleName[*r] = *s;
        SamReadGroup samRG(*r);
        samRG.Sample = *s;
        newReadGroups.push_back(samRG);
    }

    BamMultiReader reader;
    if (!reader.Open(inputFilenames)) {
        cerr << "could not open input BAM files" << endl;
        return 1;
    }

    const RefVector references = reader.GetReferenceData();

    // add read groups
    SamHeader header = reader.GetHeader();
    for (vector<SamReadGroup>::iterator r = newReadGroups.begin(); r != newReadGroups.end(); ++r) {
        header.ReadGroups.Add(*r);
    }

    BamWriter writer;
    if (!writer.Open("stdout", header.ToString(), references)) {
        cerr << "could not open BAM output stream" << endl;
        return 1;
    }

    BamAlignment al;
    while (reader.GetNextAlignment(al)) {
        if (!al.EditTag("RG", "Z", filenameToReadGroup[al.Filename])) {
            cerr << "could not add or edit RG tag on alignment " << al.Name << endl;
            return 1;
        }
        writer.SaveAlignment(al);
    }

    reader.Close();
    writer.Close();

    return 0;

}
