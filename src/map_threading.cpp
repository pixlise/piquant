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

#include "map_threading.h"


#include <sstream>
#include <mutex>
#include <thread>
#include <chrono>

#include "read_spectrum_file.h"
#include "fpSetupConditions.h"
#include "quantCombineSpectra.h"
#include "quantWriteMap.h"
#include "quantWriteResults.h"
#include "quantUnknown.h"
#include "read_EMSA_PIXL.h"
#include "setup_spectrum_parameters.h"
#include "time_code.h"
#include "read_PIXLISE_spectrum.h"


class SpectrumMapJob
{
public:
    SpectrumMapJob(
        const string &map_spec_file,

        const XRFconditionsInput &condStruct_config,

        const ARGUMENT_LIST &arguments,
        const bool oxidesOutput,

        const XraySpectrum &configSpectrum,

        const vector <ElementListEntry> &element_list,

        int jobId,

        int sequence_number,

        const string &pmcSpecifier  // If not empty, we assume the map_spec_file is a PIXLISE binary
                                    // file and need to process the given PMC in there. Otherwise
                                    // we process as before. Can be a number or number,A or number,B
                                    // to specify what detector to read. This is a temporary measure
                                    // while we still operate on MSAs that have 1 column (FM data)
    ) :
        _map_spec_file(map_spec_file),
        _condStruct_config(condStruct_config),
        _arguments(arguments),
        _oxidesOutput(oxidesOutput),
        _configSpectrum(configSpectrum),
        _element_list(element_list),
        _jobId(jobId),
        _sequence_number(sequence_number),
        _pmcSpecifier(pmcSpecifier)
    {
        _logger.setf( ios::fixed, ios::floatfield );
        _logger.precision(2);
    }

    int getJobId() const { return _jobId; }
    bool getError() const { return _error; }
    int getResultCode() const { return _result_code; }
    const ostringstream &getMapOutput() const { return _map_row; }
    const ostringstream &getResultString() const { return _logger; }
    const string &getSpectrumFile() const { return _map_spec_file; }
    string getRunTimeSec() const { return _runtimeSec; }

    void run()
    {
        // So we can wrap the run with some timing code...
        auto tm = time_code("mapSpectrum", false);

        runInternal();

        ostringstream str;
        str.precision(3);
        str << tm.elapsedSince(false);
        _runtimeSec = str.str();
    }

private:

    void runInternal()
    {
        // Spectrum we've read
        vector <XraySpectrum> spectrum_vec;

        XraySpectrum singleSpectrum = _configSpectrum;

        int result = 0;
        _error = false;

        XRFconditionsInput condStruct_Map;

        if(!_pmcSpecifier.empty() && _map_spec_file.length() > 4 && _map_spec_file.substr(_map_spec_file.length()-4) == ".bin")
        {
            // We're reading a PIXLISE binary file, and processing the spectra for a given PMC in there
            result = read_PIXLISE_spectrum(_logger, _map_spec_file, _pmcSpecifier, spectrum_vec, condStruct_Map.conditionsVector, condStruct_Map.optic_file_name );
            if ( result != 0 )
            {
                _logger << "read_PIXLISE_spectrum failed, result = " << result << " file " << endl;
                _error = true;
                _result_code = -1;
                return;
            }
        }
        else
        {
            result = read_spectrum_file( _logger, _map_spec_file, spectrum_vec, condStruct_Map );
            if ( result != 0 )
            {
                _logger << "read_spectrum_file failed, result = " << result << " file " << endl;
                _error = true;
                _result_code = -1;
                return;
            }
        }

        //  Set up energy calibration, background parameters, and measurement conditions
        setup_spectrum_parameters( _arguments, _configSpectrum.calibration(), spectrum_vec,
                _condStruct_config, condStruct_Map, _logger );

//_logger << "Read " << spectrum_vec.size() << " spectra, displaying [0]:" << endl;
//_logger << spectrum_vec[0].toString() << endl;

        if( spectrum_vec.size() <= 0 ) {
            _logger << "No spectra in file " << _map_spec_file << endl;
            _error = true;
            _result_code = -1;
            return;
        } else {
            //  Combine the spectrum information from several detectors (or the selected detector) into the variable where they will be used
            //      NB: quantCombineSpectra modifies the spectra in the input list to match them to a single energy axis
            //          for proper plotting
            result = quantCombineSpectra( spectrum_vec, singleSpectrum, _arguments.detector_select );
            if( result < 0 ) {
                _error = true;
                _result_code = -1;
                return;
            };
        }
        singleSpectrum.seq_number( _sequence_number );

        if( ! singleSpectrum.calibration().good() ) {
            _logger << "Bad energy calibration, can't quantify spectrum." << endl;
            _error = true;
            _result_code = -1;
            return;
        }
        if( singleSpectrum.live_time() <= 0 ) {
            _logger << "*** Error - live time is bad, can't quantify spectrum." << endl;
            _error = true; //  Plot can be vs channels, all others are not possible without calibration
            _result_code = 0;
            return;
        }
//tm.split("setup_spectrum");
        //  Set up new instrument measurement conditions
        XRFconditions mapConditions;

        result = fpSetupConditions ( condStruct_Map, mapConditions );
        if( result < 0 ) {
            _logger << "fpSetupConditions failed, result " << result;
            _logger << "   error in parameter with keyword " << get_EMSA_keyword( -(result+100) ) << endl;
            _error = true;
            _result_code = -500 + result;
            return;
        };
        _logger << endl;

//tm.split("setup_conditions");
        XrayMaterial unknown;
        result = quantUnknown( unknown, _element_list, mapConditions, singleSpectrum, _arguments.calibration_file, _logger );
        if ( result < 0 ) {
            _logger << "quantUnknown failed, result = " << result << "   file " << _arguments.spectrum_file << endl;
            _error = true;
        };
//tm.split("quantUnknown");
    //  Write full results to output file and put results in element list for map and calibration files
        float element_sum = 0;
        //  Normalize result if argument is not zero
        if( _arguments.normalization > 0 ) unknown.normalize( _arguments.normalization / 100 );
        result = quantWriteResults( unknown, mapConditions.detector, _element_list,
                            singleSpectrum, _oxidesOutput, _logger, element_sum );
        if ( result != 0 )
        {
            _logger << "quantWriteResults failed, result = " << result << endl;
            _error = true;
            _result_code = -1;
            return;
        }
//tm.split("quantWriteResults");

        // Save map row output locally
        quantWriteMapRow(_map_row,
            _arguments.quant_map_outputs, // TIMTIME: This was upper_trim'd twice, once when saved into ARGUMENT_LIST and once when calling quantWriteMap
            _element_list,
            mapConditions.detector,
            singleSpectrum, element_sum);

//tm.split("quantWriteMapRow");
        _result_code = 0;
    }

private:
    const string _map_spec_file;

    const XRFconditionsInput _condStruct_config;

    const ARGUMENT_LIST _arguments;
    const bool _oxidesOutput;

    const XraySpectrum _configSpectrum;

    vector <ElementListEntry> _element_list;

    const int _jobId;

    const int _sequence_number;

    const string _pmcSpecifier;

// Outputs
    ostringstream _logger;
    ostringstream _map_row;
    int _result_code;
    bool _error;
    string _runtimeSec;
};


class MTSpectrumMapJobList
{
public:
    void add(SpectrumMapJob *job)
    {
        const std::lock_guard<std::mutex> lock(_mutex);
        _jobs.push_back(job);
    }

    bool empty()
    {
        bool empty = false;
        {
            const std::lock_guard<std::mutex> lock(_mutex);
            empty = _jobs.empty();
        }
        return empty;
    }

    size_t size()
    {
        return _jobs.size();
    }

    SpectrumMapJob *remove()
    {
        SpectrumMapJob *result = 0;

        {
            const std::lock_guard<std::mutex> lock(_mutex);

            if(_jobs.size() > 0)
            {
                result = _jobs[_jobs.size()-1];
                _jobs.pop_back();
            }
        }

        return result;
    }

    void orderByMapFileName(const vector<string> &order)
    {
        const std::lock_guard<std::mutex> lock(_mutex);

        if(order.size() != _jobs.size())
        {
            cout << "ERROR: orderByMapFileName: ordering list size is different to job list size" << endl;
            return;
        }

        vector<SpectrumMapJob *> orderedJobs;
        for(auto itOrd = order.rbegin(); itOrd != order.rend(); itOrd++)
        {
            for(auto it = _jobs.begin(); it != _jobs.end(); it++)
            {
                if(*itOrd == (*it)->getSpectrumFile())
                {
                    orderedJobs.push_back(*it);
                    _jobs.erase(it);
                    break;
                }
            }
        }

        _jobs = orderedJobs;
    }

private:
    vector<SpectrumMapJob *> _jobs;
    std::mutex _mutex;
};

MTSpectrumMapJobList _mapJobQ;
MTSpectrumMapJobList _mapOutputQ;
vector<string> _mapFileOrder;

bool _mapJobRunning = false;

void setMapJobRunning(bool mapJobRunning)
{
    _mapJobRunning = mapJobRunning;
}

void outputMapFile(ostream &logger, const ARGUMENT_LIST &arguments, const vector <ElementListEntry> &element_list, const bool oxidesOutput)
{
    // We've run through, if we have any outputs, save to the output map file
    if(_mapOutputQ.empty())
    {
        logger << "No map data to output!" << endl;
    }
    else
    {
        // Write the header
        std::ofstream fout(arguments.map_file);
        std::ofstream logout(arguments.map_file+"_log.txt");
        if(!fout)
        {
            logger << "Failed open: " << arguments.map_file << " for writing." << endl;
        }
        else
        {
            // TIMTIME: What title should we put here?
            quantWriteMapHeader(fout, "Insert Title Here", arguments.quant_map_outputs, element_list, oxidesOutput);

            // Order it so we output lines in the same order we read the spectra in
            _mapOutputQ.orderByMapFileName(_mapFileOrder);

            SpectrumMapJob *job = _mapOutputQ.remove();
            while(job)
            {
                if(job->getError())
                {
                    logger << "Map row for: " << job->getJobId() << " had ERROR! Result code: " << job->getResultCode() << endl;
                }

                logout << "=================================================================" << endl;
                logout << "= " << job->getSpectrumFile() << " error=" << (job->getError() ? "true" : "false") << " result=" << job->getResultCode() << " runtime: " << job->getRunTimeSec() << "sec" << endl;
                logout << "=================================================================" << endl;

                logout << job->getResultString().str() << endl << endl;

                if(!job->getError())
                {
                    fout << job->getMapOutput().str();
                }

                // Done with this!
                delete job;
                job = _mapOutputQ.remove();
            }

            logger << "Map file written to " << arguments.map_file << endl;
            logger << "          map quantitative output options ";
            if( arguments.quant_map_outputs.length() <= 0 ) logger << "default (percents only)" << endl;
            else logger << arguments.quant_map_outputs << endl;
        }
    }
}


// mapSpectrum(quantUnknown) took: 4.1769 sec
//             quantUnknown(start) took: 0.0013 sec
// mapSpectrum(quantWriteMap) took: 0.0038 sec

void queueMapSpectrum(const std::string &map_spec_file,

    const XRFconditionsInput &condStruct_config,

    const ARGUMENT_LIST &arguments,
    const bool oxidesOutput,

    const XraySpectrum &configSpectrum,
    int n_map_spectra,

    vector <ElementListEntry> &element_list,

    int sequence_number,
    const string &pmcSpecifier)
{
    _mapFileOrder.push_back(map_spec_file);

    auto job = new SpectrumMapJob(
        map_spec_file,

        condStruct_config,

        arguments,
        oxidesOutput,

        configSpectrum,

        element_list,

        _mapJobQ.size()+1,

        sequence_number,
        pmcSpecifier
        );

    cout << "Queued: \"" << map_spec_file << "\", pmc spec: \"" << pmcSpecifier << "\"" << endl;
    _mapJobQ.add(job);
}

#define DBG_THREAD 1

void processMapJob()
{
#ifdef DBG_THREAD
    auto id = std::this_thread::get_id();
    cout << id << " processMapJob start" << endl;
#endif

    while(_mapJobRunning || !_mapJobQ.empty())
    {
        // Get it
        auto job = _mapJobQ.remove();

        if(!job)
        {
#ifdef DBG_THREAD
            cout << id << " Waiting for map job!" << endl;
#endif

            // Nothing to do, wait around
            std::this_thread::sleep_for(std::chrono::milliseconds(50));
        }
        else
        {
#ifdef DBG_THREAD
            cout << id << " Dequeued map job: " << job->getJobId() << endl;
#endif

            // Process
            job->run();

#ifdef DBG_THREAD
            cout << id << " Job ran: " << job->getJobId() << endl;
#endif

            // Save output results
            _mapOutputQ.add(job);

#ifdef DBG_THREAD
            cout << id << " Output saved for job: " << job->getJobId() << endl;
#endif
        }
    }

#ifdef DBG_THREAD
    cout << id << " processMapJob end" << endl;
#endif
}
