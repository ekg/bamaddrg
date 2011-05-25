#include <map>
#include <vector>
#include <string>
#include <getopt.h>
#include <stdlib.h>
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
         << "    -u, --uncompressed write uncompressed BAM output" << endl
         << "    -s, --sample NAME  optionally apply this sample name to the preceeding BAM file" << endl
         << "    -d, --delete NAME  removes this sample name and all associated RGs from the header" << endl
         << "    -r, --read-group GROUP  optionally apply this read group to the preceeding BAM file" << endl
         << "    -R, --region REGION  limit alignments to those in this region (chr:start..end)" << endl
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

void setRegion(BamMultiReader& reader, string& regionStr) {

    // parse the region string
    if (!regionStr.empty()) {

        map<string, int> refLength;
        map<string, int> refID;

        int id = 0;
        RefVector references = reader.GetReferenceData();
        for (RefVector::iterator r = references.begin(); r != references.end(); ++r) {
            refLength[r->RefName] = r->RefLength;
            refID[r->RefName] = id++;
        }

        // parse the region string
        string startSeq;
        int startPos;
        int stopPos;

        size_t foundFirstColon = regionStr.find(":");

        // we only have a single string, use the whole sequence as the target
        if (foundFirstColon == string::npos) {
            startSeq = regionStr;
            startPos = 0;
            stopPos = -1;
        } else {
            startSeq = regionStr.substr(0, foundFirstColon);
            size_t foundRangeDots = regionStr.find("..", foundFirstColon);
            if (foundRangeDots == string::npos) {
                startPos = atoi(regionStr.substr(foundFirstColon + 1).c_str());
                // differ from bamtools in this regard, in that we process only
                // the specified position if a range isn't given
                stopPos = startPos + 1;
            } else {
                startPos = atoi(regionStr.substr(foundFirstColon + 1, foundRangeDots - foundRangeDots - 1).c_str());
                // if we have range dots specified, but no second number, read to the end of sequence
                if (foundRangeDots + 2 != regionStr.size()) {
                    stopPos = atoi(regionStr.substr(foundRangeDots + 2).c_str()); // end-exclusive, bed-format
                } else {
                    stopPos = refLength[startSeq];
                }
            }
        }

        if (stopPos == -1) {
            stopPos = refLength[startSeq];
        }

        int startSeqRefID = refID[startSeq];

        if (!reader.LocateIndexes()) {
            cerr << "region specified, but could not open load BAM index" << endl;
            exit(1);
        } else {
            reader.SetRegion(startSeqRefID, startPos, startSeqRefID, stopPos);
        }

    }

}

int main(int argc, char** argv) {

    vector<string> inputFilenames;
    vector<string> sampleNames;
    vector<string> readGroups;

    map<string, int> samplesToDelete;
    
    string currFileName;
    string currReadGroup;
    string currSampleName;
    string regionStr;

    bool writeUncompressed = false;

    // parse command-line options
    int c;

    while (true) {

        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"bam",  required_argument, 0, 'b'},
            {"uncompressed",  no_argument, 0, 'u'},
            {"read-group", required_argument, 0, 'r'},
            {"delete", required_argument, 0, 'd'},
            {"sample", required_argument, 0, 's'},
            {"region", required_argument, 0, 'R'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hb:d:s:r:R:",
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

            case 'u':
                writeUncompressed = true;
                break;

            case 'd':
                samplesToDelete[optarg];
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

            case 'R':
                regionStr = optarg;
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
    vector<SamReadGroup> readGroupsToDelete;

    vector<string>::iterator b = inputFilenames.begin();
    vector<string>::iterator r = readGroups.begin();
    vector<string>::iterator s = sampleNames.begin();
    for (; b != inputFilenames.end(); ++b, ++r, ++s) {
        //cerr << *b << " " << *s << " " << *r << endl;
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

    setRegion(reader, regionStr);

    const RefVector references = reader.GetReferenceData();

    // add read groups
    SamHeader header = reader.GetHeader();
    for (vector<SamReadGroup>::iterator r = newReadGroups.begin(); r != newReadGroups.end(); ++r) {
        header.ReadGroups.Add(*r);
    }

    vector<SamReadGroup> rgsToDelete;
    for (SamReadGroupIterator g = header.ReadGroups.Begin(); g != header.ReadGroups.End(); ++g) {
        if (samplesToDelete.find(g->Sample) != samplesToDelete.end()) {
            rgsToDelete.push_back(*g);
        }
    }

    for (vector<SamReadGroup>::iterator r = rgsToDelete.begin(); r != rgsToDelete.end(); ++r) {
        header.ReadGroups.Remove(*r);
    }

    BamWriter writer;

    if (writeUncompressed) {
        writer.SetCompressionMode(BamWriter::Uncompressed);
    }

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
