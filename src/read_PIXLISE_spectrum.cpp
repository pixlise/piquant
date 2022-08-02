// Copyright (c) 2018-2022 California Institute of Technology (“Caltech”) and
// University of Washington. U.S. Government sponsorship acknowledged.
// All rights reserved.
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// * Neither the name of Caltech nor its operating division, the Jet Propulsion
//   Laboratory, nor the names of its contributors may be used to endorse or
//   promote products derived from this software without specific prior written
//   permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#include "read_PIXLISE_spectrum.h"

using namespace std;

// We can exclude this if it's not needed, so compilation is simpler/faster by not requiring code
// generation of protobuf serialiser, protobuf library headers/linking
#if defined EXCLUDE_PIXLISE_READER

int read_PIXLISE_spectrum(std::ostream &termOutFile,
        const std::string &spectrumPathName,
        const std::string &spectrumSelector,
        std::vector <XraySpectrum> &spectra,
        std::vector <float> &conditionsArray,
        std::string &optic_file
    )
{
    termOutFile << "This build of PIQUANT does not support reading PIXLISE binary files" << endl;
    return -1;
}

#else

#include <iostream>
#include <fstream>
#include <map>
#include "data-formats/experiment.pb.h"
#include "upper_trim.h"
#include "read_EMSA_PIXL.h"


// Defined in read_EMSA_PIXL
int parse_EMSA_description( const int index, const string &s );
// Defined in read_spectrum_file
void print_spectrum_summary(const std::vector <XraySpectrum> &spectra, std::ostream &termOutFile);

// Utility functions to make this file more readable...
string getMetaByLabel(const Experiment_Location_DetectorSpectrum &detector, const vector<string> &meta_labels, const string &label);
void getSpectrumUncompressed(const Experiment_Location_DetectorSpectrum &detector, vector<float> &out_spectrum_values);
int processMetadataValue(int pmc, size_t numSpectrumChannels, const string &label, Experiment_MetaDataType metaType, const Experiment_Location_MetaDataItem &meta,
    // We read into the following variables:
        XraySpectrum &this_spectrum,
        Spec_Aux_Info &spec_info_hold,
        bool &kev_units,
        bool &livetime_XIA,
        float &ev_ch,
        float &ev_start,
        std::vector <float> &conditionsArray,
        std::string &optic_file,
    // Logging...
        std::ostream &termOutFile
        );

struct SelectorParams
{
    SelectorParams() : _matchCount(0)
    {
    }

    SelectorParams(const std::string &readType, const std::string &detectorId) : _matchCount(0), _readType(readType), _detectorId(detectorId)
    {
    }

    string toString(int pmc) const
    {
        ostringstream o;
        o << pmc << "|" << _readType << "|" << _detectorId;
        return o.str();
    }

    int _matchCount;
    std::string _readType;
    std::string _detectorId;
};

class SpectrumMatcher
{
public:
    SpectrumMatcher() {}

    int getMatchCount(const string &readType, const string &detectorId, std::vector<SelectorParams>::const_iterator &out_it)
    {
        for(auto it = _selectors.begin(); it != _selectors.end(); it++)
        {
            if(readType == it->_readType && detectorId == it->_detectorId)
            {
                // Increment the match count & return
                it->_matchCount++;
                out_it = it;
                return it->_matchCount;
            }
        }

        // No matches!
        return 0;
    }

    bool makeSpectrumFileNameValues(std::string &out_ReadType, std::string &out_Detector) const
    {
        // If A and B vary, return Combined. If readType varies, return Mixed
        if(_selectors.empty())
        {
            return false;
        }

        auto first = _selectors.begin();
        string readType = first->_readType;
        string detectorId = first->_detectorId;

        for(auto it = _selectors.begin(); it != _selectors.end(); it++)
        {
            if(readType != "Mixed" && it->_readType != readType)
            {
                readType = "Mixed";
            }

            if(detectorId != "Combined" && it->_detectorId != detectorId)
            {
                detectorId = "Combined";
            }
        }

        out_ReadType = readType;
        out_Detector = detectorId;
        return true;
    }

    void addSelectorParams(const SelectorParams &params)
    {
        _selectors.push_back(params);
    }

    bool allMatched() const
    {
        size_t matchCount = 0;
        for(auto it = _selectors.begin(); it != _selectors.end(); it++)
        {
            matchCount += it->_matchCount;
        }

        return matchCount == _selectors.size();
    }

protected:
    std::vector<SelectorParams> _selectors;
};

typedef std::map<int, SpectrumMatcher> MatcherMap;
bool makeSelectorMatcher(const std::string &spectrumSelector, MatcherMap &out_matcher, std::string &out_matcherOverallFileName, std::ostream &termOutFile);

int read_PIXLISE_spectrum(std::ostream &termOutFile,
        const std::string &spectrumPathName,
        const std::string &spectrumSelector,
        std::vector <XraySpectrum> &spectra,
        std::vector <float> &conditionsArray,
        std::string &optic_file
    )
{
    std::string selectorPreview = spectrumSelector;
    if(selectorPreview.length() > 50)
    {
        selectorPreview = selectorPreview.substr(0, 50);
        selectorPreview += "...";
    }
    termOutFile << "Reading spectrum from file: " << spectrumPathName << " with selector: " << selectorPreview << endl;

    // Open the file
    std::ifstream fin(spectrumPathName.c_str(), std::ios::binary);
    if(!fin)
    {
        termOutFile << "Failed to read PIXLISE binary file: " << spectrumPathName << endl;
        return -1;
    }

    // Read with protobuf deserialisation code
    Experiment exp;
    if(!exp.ParseFromIstream(&fin))
    {
        termOutFile << "Failed to parse PIXLISE binary file: " << spectrumPathName << endl;
        return -1;
    }

    vector<string> meta_labels;
    for(int c = 0; c < exp.meta_labels_size(); c++)
    {
        auto label = exp.meta_labels(c);
        meta_labels.push_back(label);
    }
    vector<Experiment_MetaDataType> meta_types;
    for(int c = 0; c < exp.meta_labels_size(); c++)
    {
        meta_types.push_back(exp.meta_types(c));
    }

    MatcherMap matcherMap;
    std::string overallFileNameColumn;
    if(!makeSelectorMatcher(spectrumSelector, matcherMap, overallFileNameColumn, termOutFile))
    {
        // Logs already written if error...
        return -1;
    }

    termOutFile << "Parsed " << spectrumPathName << ", created selector matcher: " << selectorPreview << " with " << matcherMap.size() << " entries. Overall file name column: " << overallFileNameColumn << endl;

    // Get the spectra for the specified PMC
    spectra.clear();
    size_t pmcsMatched = 0;

    for(int c = 0; c < exp.locations_size(); c++)
    {
        auto loc = exp.locations(c);
        auto locPMC = atoi(loc.id().c_str());

        // Check if its one of the PMCs we were asked to look for
        auto matcherIt = matcherMap.find(locPMC);

        if(matcherIt != matcherMap.end())
        {
            // We've found a PMC we're interested in, now find all spectra specified
            // We have multiple detectors, each will have its own spectra that it read, so for a given PMC we return at least 1 spectrum...
            // This is mostly copied from read_EMSA_PIXL, as we're basically trying to supply the same data, but from a different source
            for(int det = 0; det < loc.detectors_size(); det++)
            {
                Experiment_Location_DetectorSpectrum detector = loc.detectors(det);

                auto readType = getMetaByLabel(detector, meta_labels, "READTYPE");
                auto detectorId = getMetaByLabel(detector, meta_labels, "DETECTOR_ID");

                // If we don't find the above, complain
                if(detectorId.empty() || readType.empty())
                {
                    termOutFile << "PIXLISE binary file: " << spectrumPathName << " pmc: " << locPMC << " was missing READTYPE and/or DETECTOR_ID" << endl;
                }

                // If we don't match, skip
                std::vector<SelectorParams>::const_iterator selectorMatchedIt;
                auto matchCount = matcherIt->second.getMatchCount(readType, detectorId, selectorMatchedIt);
                if(matchCount == 0)
                {
                    // Not one we're interested in
                    continue;
                }
                if(matchCount > 1)
                {
                    // Weird, we matched twice, shouldn't be!
                    termOutFile << "PIXLISE binary file: " << spectrumPathName << " pmc: " << locPMC << " readtype: " << readType << ", detectorId: " << detectorId << " was matched multiple times!" << endl;
                    return -1;
                }

                // Only consider it if it has beam location info
                if(!loc.has_beam())
                {
                    termOutFile << "PIXLISE binary file: " << spectrumPathName << " pmc: " << locPMC << " readtype: " << readType << ", detectorId: " << detectorId << " had no beam location!" << endl;
                    return -1;
                }

                auto beam = loc.beam();

                // Otherwise, we've matched this one, so add it to the spectrum list we're returning
                XraySpectrum this_spectrum;

                // Read out the spectrum values as an array of floats
                Spec_Aux_Info spec_info_hold;
                bool kev_units = false;
                bool livetime_XIA = false;
                float ev_ch = 0;
                float ev_start = 0;

                vector<float> spectrum_values;
                getSpectrumUncompressed(detector, spectrum_values);

                // Now process all the meta tags, if they exist, using the lookup
                for(int metaIdx = 0; metaIdx < detector.meta_size(); metaIdx++)
                {
                    // Get meta value
                    auto meta = detector.meta(metaIdx);

                    // Get the index to look up label/type with
                    auto label_idx = meta.label_idx();

                    // Get label and type of this metadata value
                    auto label = meta_labels[label_idx];
                    auto metaType = meta_types[label_idx];

                    // Get data for this metadata tag
                    if(processMetadataValue(locPMC, spectrum_values.size(), label, metaType, meta,
                        this_spectrum,
                        spec_info_hold,
                        kev_units,
                        livetime_XIA,
                        ev_ch,
                        ev_start,
                        conditionsArray,
                        optic_file,
                        termOutFile) != 0)
                    {
                        return -1;
                    }
                }

                // Set specific fields...
                spec_info_hold.pmc = locPMC;

                // Some comes from beam location:
                spec_info_hold.x = beam.x();
                spec_info_hold.y = beam.y();
                spec_info_hold.z = beam.z();

                spec_info_hold.i = beam.image_i();
                spec_info_hold.j = beam.image_j();

                // geom_corr is optional, the only way to tell it's set is if it's non-zero
                if(beam.geom_corr() != 0)
                {
                    conditionsArray[GEOMETRY_INDEX] = beam.geom_corr();
                }

                // Now that we've read out the meta values and spectrum counts for this PMC/detector combo, lets set up the spectrum
                if(kev_units)
                {
                    ev_ch *= 1000;
                    ev_start *= 1000;
                }
                this_spectrum.calibration(ev_start, ev_ch);
                this_spectrum.aux_info_replace(spec_info_hold);

                if(livetime_XIA)
                {
                    //  See writeup in file "JPL-XIA_PIXL_FPGA_Specification_v2.06.pdf", page 9
                    this_spectrum.header_info_change().live_time_DSPC = this_spectrum.live_time();
                    const Spec_Header_Info &header = this_spectrum.header_info();
                    if(header.triggers > 0)
                    {
                        float total_counts_in = header.events + header.overflows + header.underflows;
                        this_spectrum.live_time( header.live_time_DSPC * total_counts_in / header.triggers );
                    }
                    else if(header.live_time_DSPC != 0)
                    {
                        termOutFile << "Unexpected livetime_XIA/live_time situation found for selector: " << selectorMatchedIt->toString(locPMC) << endl;
                        return -1;
                    }
                }

                this_spectrum.meas(spectrum_values);

                // Set the file name (using the source file we saved in binary
                // file, so this should look like an MSA read to the rest of piquant!)
                string file_name = overallFileNameColumn;
                this_spectrum.file_name(file_name);

                spectra.push_back(this_spectrum);

                termOutFile << "Read spectrum for selector: " << selectorMatchedIt->toString(locPMC) << " from: \"" << spectrumPathName << "\"" << endl;
            }

            // At this point we've scanned everything for this PMC (which was matched). Check that all specifiers were matched
            if(!matcherIt->second.allMatched())
            {
                termOutFile << "Failed to match all selectors: " << spectrumSelector << " for PMC: " << locPMC << " in dataset file: \"" << spectrumPathName << "\"" << endl;
                return -1;
            }

            // If all PMCs have been matched, we can stop here
            pmcsMatched++;
            if(pmcsMatched >= matcherMap.size())
            {
                termOutFile << "Read " << spectra.size() << " spectra specified by: " << selectorPreview << " from: \"" << spectrumPathName << "\" successfully" << endl;
                print_spectrum_summary(spectra, termOutFile);
                return 0;
            }
        }
    }

    // Didn't match something anyway...
    termOutFile << "Failed to match all selectors: " << spectrumSelector << " in dataset file: \"" << spectrumPathName << "\"" << endl;
    return -1;
}

// TODO: should be a unit tested utility function!
void splitString(const std::string &str, char delim, vector<string> &out_result)
{
    out_result.clear();

    istringstream inStr(str);
    string s;
    while(getline(inStr, s, delim))
    {
        out_result.push_back(s);
    }
}

bool makeSelectorMatcher(const std::string &spectrumSelectorRaw, MatcherMap &out_matcher, std::string &out_matcherOverallFileName, std::ostream &outLog)
{
    out_matcher.clear();

    // New optional feature, string can start with a tag, separated from the rest by a :
    // So if we have a : we read a tag out
    // Specifically the tag is useful for something like an ROI ID
    std::string tag;
    std::string spectrumSelector = spectrumSelectorRaw;

    auto tagSepPos = spectrumSelector.find(":");
    if(tagSepPos != string::npos)
    {
        tag = spectrumSelector.substr(0, tagSepPos);
        spectrumSelector = spectrumSelector.substr(tagSepPos+1, string::npos);
    }

    // We're expecting strings of the form:
    // 15|Normal|A,15|Normal|B
    vector<string> specs;
    splitString(spectrumSelector, ',', specs);

    // Run through each one and break them further
    for(auto it = specs.begin(); it != specs.end(); it++)
    {
        vector<string> bits;
        splitString(*it, '|', bits);

        // Should have 3 parts: PMC, READTYPE, DETECTOR_ID
        if(bits.size() != 3)
        {
            outLog << "Failed to parse pmcs file line: " << spectrumSelector << endl;
            return false;
        }

        // Make sure 1st is a number
        auto pmc = atoi(bits[0].c_str());
        if(pmc <= 0)
        {
            outLog << "Failed to parse PMC=" << bits[0] << " on pmcs file line: " << spectrumSelector << endl;
            return false;
        }

        // Second should be a read type
        if(bits[1] != "Normal" && bits[1] != "Dwell" && bits[1] != "BulkSum" && bits[1] != "MaxValue")
        {
            outLog << "Failed to parse READTYPE=" << bits[1] << " on pmcs file line: " << spectrumSelector << endl;
            return false;
        }

        // Finally, expecting detector ID A or B
        if(bits[2] != "A" && bits[2] != "B")
        {
            outLog << "Failed to parse DETECTOR_ID=" << bits[2] << " on pmcs file line: " << spectrumSelector << endl;
            return false;
        }

        auto existingIt = out_matcher.find(pmc);
        if(existingIt == out_matcher.end())
        {
            // Doesn't exist yet, add an empty vector here
            out_matcher[pmc] = SpectrumMatcher();
        }

        out_matcher[pmc].addSelectorParams(SelectorParams(bits[1], bits[2]));
        //outLog << "Added matcher for PMC: " << pmc << " with: " << *it << endl;
    }

    if(out_matcher.empty())
    {
        outLog << "Failed to find any matching information on pmcs file line: " << spectrumSelector << endl;
        return false;
    }

    // Run through all matchers, for the case where there are multiple PMCs, we want to scan all matchers for READTYPE and DETECTOR_ID
    // to form a file name column entry for anything matched by this matcher map.
    // In other words, if match string is: 41|Normal|A,41|Normal|B,42|Normal|A,42|Normal|B, we want
    // to form an overall file name of Normal_Combined, but if 41|Normal|A,42|Normal|A,43|Normal|A we want
    // Normal_A.
    // If there is only one matcher, this is taken care of already
    std::string overallReadType, overallDetectorId; 

    for(auto it = out_matcher.begin(); it != out_matcher.end(); it++)
    {
        std::string readType, detectorId;
        if(it->second.makeSpectrumFileNameValues(readType, detectorId))
        {
            if(overallReadType.empty())
            {
                overallReadType = readType;
            }
            else if(readType != overallReadType)
            {
                overallReadType = "Mixed";
            }

            if(overallDetectorId.empty())
            {
                overallDetectorId = detectorId;
            }
            else if(detectorId != overallDetectorId)
            {
                overallDetectorId = "Combined";
            }
        }
    }

    out_matcherOverallFileName = overallReadType+"_"+overallDetectorId;

    // If there is a tag, add it
    if(!tag.empty())
    {
        out_matcherOverallFileName += "_" + tag;
    }

    return true;
}

string getMetaByLabel(const Experiment_Location_DetectorSpectrum &detector, const vector<string> &meta_labels, const string &label)
{
    for(int metaIdx = 0; metaIdx < detector.meta_size(); metaIdx++)
    {
        auto meta = detector.meta(metaIdx);
        if(meta_labels[meta.label_idx()] == label)
        {
            return meta.svalue();
        }
    }

    return "";
}


void getSpectrumUncompressed(const Experiment_Location_DetectorSpectrum &detector, vector<float> &out_spectrum_values)
{
    // We have to read in the values, but keep in mind they're compressed...
    // at the moment the only way we store them is zero-run-length encoding
    bool lastWas0 = false;
    for(int ch = 0; ch < detector.spectrum_size(); ch++)
    {
        auto val = detector.spectrum(ch);
        if(val == 0)
        {
            // Next value will tell us how many
            lastWas0 = true;
        }
        else
        {
            if(lastWas0)
            {
                // this is telling us how many 0's there are
                for(int z = 0; z < val; z++)
                {
                    out_spectrum_values.push_back(0);
                }
                lastWas0 = false;
            }
            else
            {
                // it's just a value, store it
                out_spectrum_values.push_back(val);
            }
        }
    }
}

int processMetadataValue(int pmc, size_t numSpectrumChannels, const string &label, Experiment_MetaDataType metaType, const Experiment_Location_MetaDataItem &meta,
        // We read into the following variables:
        XraySpectrum &this_spectrum,
        Spec_Aux_Info &spec_info_hold,
        bool &kev_units,
        bool &livetime_XIA,
        float &ev_ch,
        float &ev_start,
        std::vector <float> &conditionsArray,
        std::string &optic_file,
    // Logging...
        std::ostream &termOutFile
        )
{
    bool handled = true;

    // Some labels we ignore, some we store data for... See read_EMSA_spectrum for what we need to store, this is just
    // another data format that intends to supply the same values as the MSA files did.
    if(metaType == Experiment_MetaDataType_MT_STRING)
    {
        auto value = meta.svalue();

        // A bunch of these wants to interpret the string as a float...
        float fValue = 0;
        bool isFloat = false;
        {
            istringstream valStream(value);
            valStream >> fValue;
            isFloat = !(!valStream);
        }
        int iValue = 0;
        bool isInt = false;
        {
            istringstream valStream(value);
            valStream >> iValue;
            isInt = !(!valStream);
        }

        if(label == "FORMAT")
        {
            if(value != "EMSA/MAS spectral data file")
            {
                termOutFile << "Unexpected data format found for PMC: " << pmc << ": " << value << endl;
                return -1;
            }
        }
        else if(label == "VERSION")
        {
            if(value != "TC202v2.0 PIXL")
            {
                termOutFile << "Unexpected data version found for PMC: " << pmc << ": " << value << endl;
                return -1;
            }
        }
        else if(label == "SIGNALTYPE")
        {
            if(value != "XRF")
            {
                termOutFile << "Unexpected signal type found for PMC: " << pmc << ": " << value << endl;
                return -1;
            }
        }
        else if(label == "DATATYPE")
        {
            if(value != "Y" && value != "YY")
            {
                termOutFile << "Unexpected data type found for PMC: " << pmc << ": " << value << endl;
                return -1;
            }
        }
        else if(label == "COMMENT")
        {
            spec_info_hold.comments.push_back(value);
        }
        else if(label == "TITLE")
        {
            spec_info_hold.titles.push_back(value);
        }
        else if(label == "DATE")
        {
            spec_info_hold.date = value;
        }
        else if(label == "TIME")
        {
            spec_info_hold.time = value;
        }
        else if(label == "OWNER")
        {
            spec_info_hold.owner = value;
        }
        else if(label == "NPOINTS")
        {
            // Make sure it matches the spectrum point count we read
            if(!isInt || numSpectrumChannels != iValue)
            {
                termOutFile << "Unexpected NPOINTS found for PMC: " << pmc << ": " << value << ", expected: " << numSpectrumChannels << endl;
                return -1;
            }
        }
        else if(label == "NCOLUMNS")
        {
            /*
            // Make sure it matches the spectrum point count we read
            if(!isInt || loc.detectors_size() != iValue)
            {
                termOutFile << "Unexpected NCOLUMNS found for PMC: " << pmc << ": " << value << ", expected: " << loc.detectors_size() << endl;
                return -1;
            }
            */
            // We're expecting 1 column, because we read 1 spectrum per detector
            // As a work-around temporarily we allow 1 or 2 because we still haven't resolved wether the idea is piquant
            // loads the 2 detectors and ONLY outputs 1 row per PMC vs 1 per detector and PMC
            if(!isInt || (iValue != 1 && iValue != 2))
            {
                termOutFile << "Unexpected NCOLUMNS found for PMC: " << pmc << ": " << value << ", expected: 1" << endl;
                return -1;
            }
        }
        else if(label == "XUNITS")
        {
            value = upper_trim(value);
            if(value == "EV")
            {
                kev_units = false;
            }
            else if(value == "KEV")
            {
                kev_units = true;
            }
            else
            {
                termOutFile << "Unexpected x-units found for PMC: " << pmc << ": " << value << endl;
                return -1;
            }
        }
        else if(label == "YUNITS")
        {
            if(value != "COUNTS")
            {
                termOutFile << "Unexpected y-units found for PMC: " << pmc << ": " << value << endl;
                return -1;
            }
        }
        // IGNORED in read_EMSA_PIXL: XLABEL
        // IGNORED in read_EMSA_PIXL: YLABEL
        else if(label == "OPTICFILE")
        {
            // If it's a number, act differently...
            if(!isInt)
            {
                optic_file = value;
                conditionsArray[TEST_OPTIC_TYPE_INDEX] = 4;
            }
            else
            {
                conditionsArray[TEST_OPTIC_TYPE_INDEX] = iValue;
            }
        }

        // These were multi-column values in MSA files, but the dataset converter already takes this into acount and
        // reads them into 1 value per detector spectrum so we can just parse them direct as single values
        // NOTE: these are all read as floats, so check for that!
        else if(!isFloat && (
                (label == "TRIGGERS") ||
                (label == "EVENTS") ||
                (label == "OVERFLOWS") ||
                (label == "UNDERFLOWS") ||
                (label == "BASE_EVENTS") ||
                (label == "RESETS") ||
                (label == "OVER_ADCMAX")
            )
        )
        {
            termOutFile << "Unexpected " << label << " value found for PMC: " << pmc << ": " << value << endl;
            return -1;
        }
        else if(label == "TRIGGERS")
        {
            //triggers_line = line_number;
            livetime_XIA = true;
            this_spectrum.header_info_change().triggers = fValue;
        }
        else if(label == "EVENTS")
        {
            this_spectrum.header_info_change().events = fValue;
        }
        else if(label == "OVERFLOWS")
        {
            this_spectrum.header_info_change().overflows = fValue;
        }
        else if(label == "UNDERFLOWS")
        {
            this_spectrum.header_info_change().underflows = fValue;
        }
        else if(label == "BASE_EVENTS")
        {
            this_spectrum.header_info_change().baseline_samples = fValue;
        }
        else if(label == "RESETS")
        {
            this_spectrum.header_info_change().preamp_resets = fValue;
        }
        else if(label == "OVER_ADCMAX")
        {
            this_spectrum.header_info_change().saturates = fValue;
        }
        // Single values...
        else if(label == "DETECTOR_ID")
        {
            spec_info_hold.det_ID = value;
        }
        else if(label == "IPOSITION")
        {
            // Read from beam location instead...
            //spec_info_hold.i = fValue;
        }
        else if(label == "JPOSITION")
        {
            // Read from beam location instead...
            //spec_info_hold.j = fValue;
        }
        else
        {
            handled = false;
        }
    }
    else if(metaType == Experiment_MetaDataType_MT_INT)
    {
        auto ivalue = meta.ivalue();
/*      Took this out, we set it directly anyway...
        if(label == "PMC")
        {
            spec_info_hold.pmc = value;
        }
        else*/ if(label == "RTT")
        {
            spec_info_hold.rtt = ivalue;
        }
        else
        {
            handled = false;
        }
    }
    else if(metaType == Experiment_MetaDataType_MT_FLOAT)
    {
        auto fvalue = meta.fvalue();

        if(label == "XPERCHAN")
        {
            ev_ch = fvalue;
        }
        else if(label == "OFFSET")
        {
            ev_start = fvalue;
        }
        else if(label == "LIVETIME")
        {
            this_spectrum.live_time( fvalue );
        }
        else if(label == "REALTIME")
        {
            this_spectrum.real_time( fvalue );
        }
        else if(label == "XPOSITION")
        {
            // Read from beam location instead...
            //spec_info_hold.x = fvalue;
        }
        else if(label == "YPOSITION")
        {
            // Read from beam location instead...
            //spec_info_hold.y = fvalue;
        }
        else if(label == "ZPOSITION")
        {
            // Read from beam location instead...
            //spec_info_hold.z = fvalue;
        }
        else
        {
            handled = false;
        }
    }
    else
    {
        handled = false;
    }

    if(!handled)
    {
        // Any that weren't handled above, we try to handle here
        for(int i = 0; i<XRF_PARAMETER_LAST; i++)
        {
            if(i == TEST_OPTIC_TYPE_INDEX)
            {
                continue;  //  handled above since it might be a file name
            }

            string upper_test = upper_trim(get_EMSA_keyword(i));
            // There can be 1 or 2 # at the start...
            if(upper_test[0] == '#')
            {
                upper_test = upper_test.substr(1);
            }
            if(upper_test[0] == '#')
            {
                upper_test = upper_test.substr(1);
            }

            if(label == upper_test)
            {
                auto value = meta.svalue();
                if(metaType != Experiment_MetaDataType_MT_STRING)
                {
                    termOutFile << "Expected string for " << label << " when reading conditionsArray for PMC: " << pmc << ": " << value << endl;
                    return -1;
                }

                istringstream vstr(value);
                vstr >> conditionsArray[i];
                if(!vstr)
                {
                    conditionsArray[i] = parse_EMSA_description(i, value);
                    if(conditionsArray[i] < 0)
                    {
                        termOutFile << "Failed to read " << label << " as conditionsArray for PMC: " << pmc << ": " << value << endl;
                        return -1;
                    }
                }

                if(i == TUBE_CURRENT_INDEX)
                {
                    conditionsArray[i] /= 1000;   //  convert from microAmps to milliAmps
                }

                // Found it, no longer need to process...
                break;
            }
        }
    }

    return 0;
}

#endif
