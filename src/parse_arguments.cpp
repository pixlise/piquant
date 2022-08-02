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

//
//  parse_arguments.cpp
//  PIQUANT
//
//  Created by W. T. Elam on 1/14/2017.
//  Copyright (c) 2017 APL/UW. All rights reserved.
//

#include <iostream>
#include <sstream>
#include "parse_arguments.h"
#include "parse_records.h"
#include "upper_trim.h"
#include "XRFconstants.h"

//  Written Jan. 15, 2017
//      Parse argument list for PIQUANT Subprocess
//      Write helpful information to cout if errors
//  Modified May 9, 2017 to match documentation
//      Change order of sub-command help list
//      Update other help output to reflect actual operation
//  Modified Sept. 30, 2017
//      Return invalid arguments error messages in arguments.invalid_arguments instead of writing them to cout
//          (cout is not redirected yet)
//  Modified Nov. 8, 2017
//      Add ems sub-command to convert output of SEND_SDD_DATA command to EDR (csv) format
//  Modified Jan. 3, 2017
//      Move tab, single quote, double quote, blank, comma, and underscore definitions to XRFconstants.h
//  Modified Mar. 2, 2018
//      Added -q option to specify outputs to map file
//  Modified July 27, 2018
//      Change background arguments to vector <float> and handle as many as follow -b
//  Modified May 15, 2019
//      Clear out defaults for option arguments if option letter is found (so defaults can be set at allocation)
//      remove upper_trim for quant map options (-q)
//  Modified May 16, 2019
//      Add -f and -g options, capability to turn off adjustments to energy calibration and detector resolution in fits
//      Add -s option to convolve Compton scatter components with a Gaussian (detector resolution)  brute force, very expensive in compute time
//  Modified May 21, 2019
//      Implement Evaluate action, with correct argument list
//  Modified July 3, 2019
//      Change default for convolution of Compton scatter to true (for use with calculated instead of SNIP background)
//  Modified Oct. 22, 2020  Add optic response sub-command
//  Modified Oct. 26, 2020  Add actual sub-command found to invalid sub-command message
//  Modified Nov. 4, 2020   Select standard from standards file using -s option (by number or name)  [change -s for Compton convolve to -v]
//                          Add -w option to write evaluation file during calibration for plotting and debugging
//  Modified Jan. 13, 2021  Change -w option to -v, use -w for evaluate minimum weight
//  Modified Feb. 2, 2021   Add -a for crossover of SNIP bkg fit (fit below crossover and not above)
//  Modified May 10, 2021   Add -bh and -bx background options, eliminate -a option
//  Modified May 14, 2021   Move shelf factor and slope to XrayDetector and control via -T option
//  Modified June 27, 2021  Add command line option to normalize element sum to 100% (or any value)
//  Modified July 9, 2021   Add command line option to change Fe oxide ratio (-Fe)

using namespace std;

//  Utility function to figure out command and write list if not found
int parse_command( const string &sub_command, PIQUANT_SUBCOMMAND &cmd );

int parse_arguments( const int argc, const char * argv[],
            PIQUANT_SUBCOMMAND &cmd, ARGUMENT_LIST &arguments ) {

    int result = 0;
    arguments.eV_ch = 0;
    if( argc > 1 ) {
        string sub_command( argv[1] );
        result = parse_command( sub_command, cmd );
        if( result < 0 ) return result;
    } else {
        string sub_command;
        result = parse_command( sub_command, cmd );
        if( result < 0 ) return result;
    }

    //  Break up the arguments into the file list according to the sub-command
    //  enum PIQUANT_SUBCOMMAND { ENERGY_CAL, PLOT, PRIMARY, CALCULATE, CALIBRATE, QUANTIFY, EVALUATE, MAP };
    int term_file_index = 1;
    switch( cmd ) {
        case ENERGY_CAL:
            term_file_index = 4;
            if( argc < term_file_index ) {
                cout << endl;
                cout << "Not enough arguments for energy calibration." << endl;
                cout << "   Spectrum file" << endl;
                cout << "   Element list for one or two largest peaks (K lines only at present)" << endl;
                cout << "     (comma or space separated, probably needs to be in quotes, no tabs)" << endl;
                cout << endl;
                return -2001;
            } else {
                string spec_file( argv[2] );
                arguments.spectrum_file = spec_file;
                string elements( argv[3] );
                arguments.element_list = elements;
            }
            break;
        case PLOT:
            term_file_index = 4;
            if( argc < term_file_index ) {
                cout << endl;
                cout << "Not enough arguments for plot." << endl;
                cout << "   Spectrum file (or CSV file)" << endl;
                cout << "   Plot file (required, this is the plot output, CSV format)" << endl;
                cout << endl;
                return -2002;
            } else {
                string spec_file( argv[2] );
                arguments.spectrum_file = spec_file;
                string plot_file( argv[3] );
                arguments.plot_file = plot_file;
            }
            break;
        case PRIMARY:
            term_file_index = 4;
            if( argc < term_file_index ) {
                cout << endl;
                cout << "Not enough arguments for primary spectrum calculation." << endl;
                cout << "   Configuration file" << endl;
                cout << "   Plot file (required, calculated spectrum output, CSV format)" << endl;
                cout << endl;
                return -2003;
            } else {
                string config( argv[2] );
                arguments.configuration_file = config;
                string plot_file( argv[3] );
                arguments.plot_file = plot_file;
            }
            break;
        case CALCULATE:
            term_file_index = 5;
            if( argc < term_file_index ) {
                cout << endl;
                cout << "Not enough arguments for spectrum calculation." << endl;
                cout << "   Configuration file" << endl;
                cout << "   Standards file (only the first standard is processed)" << endl;
                cout << "   Plot file (required, calculated spectrum output, CSV format)" << endl;
                cout << endl;
                return -2004;
            } else {
                string config( argv[2] );
                arguments.configuration_file = config;
                string stds( argv[3] );
                arguments.standards_file = stds;
                string plot_file( argv[4] );
                arguments.plot_file = plot_file;
            }
            break;
        case CALIBRATE:
            term_file_index = 6;
            if( argc < term_file_index ) {
                cout << endl;
                cout << "Not enough arguments for quantitative calibration." << endl;
                cout << "   Configuration file" << endl;
                cout << "   Standards file" << endl;
                cout << "   Calibration file (overwritten)" << endl;
                cout << "   Element fit control list (optional but must be present if plot file or any options)" << endl;
                cout << endl;
                return -2005;
            } else {
                string config( argv[2] );
                arguments.configuration_file = config;
                string stds( argv[3] );
                arguments.standards_file = stds;
                string cal( argv[4] );
                arguments.calibration_file = cal;
                string elements( argv[5] );
                arguments.element_list = elements;
           }
            break;
        case QUANTIFY:
            term_file_index = 6;
            if( argc < term_file_index ) {
                cout << endl;
                cout << "Not enough arguments for quantification." << endl;
                cout << "   Configuration file" << endl;
                cout << "   Calibration file" << endl;
                cout << "   Spectrum file" << endl;
                cout << "   Element list for quantification (required, see user manual)" << endl;
                cout << endl;
                return -2006;
            } else {
                string config( argv[2] );
                arguments.configuration_file = config;
                string cal( argv[3] );
                arguments.calibration_file = cal;
                string spec_file( argv[4] );
                arguments.spectrum_file = spec_file;
                string elements( argv[5] );
                arguments.element_list = elements;
                //  See if there is a plot file in the argument list
                if( argc > term_file_index ) {
                    string plot_file( argv[term_file_index] );
                    arguments.plot_file = plot_file;
                    term_file_index++;
                }
            }
            break;
        case EVALUATE:
            term_file_index = 7;
            if( argc < term_file_index ) {
                cout << endl;
                cout << "Not enough arguments for evaluate." << endl;
                cout << "   Configuration file" << endl;
                cout << "   Standards file (each standard in this file is processed as an unknown)" << endl;
                cout << "   Calibration file (used for quantification of each standard)" << endl;
                cout << "   Element list for quantification (added to element list derived for each standard)" << endl;
                cout << "   Map file (overwritten, contains results of quantifying each standard in standards file)" << endl;
                cout << endl;
                return -2007;
            } else {
                string config( argv[2] );
                arguments.configuration_file = config;
                string stds( argv[3] );
                arguments.standards_file = stds;
                string cal( argv[4] );
                arguments.calibration_file = cal;
                string elements( argv[5] );
                arguments.element_list = elements;
                string map_file( argv[6] );
                arguments.map_file = map_file;
            }
            break;
        case MAP:
            term_file_index = 7;
            if( argc < term_file_index ) {
                cout << endl;
                cout << "Not enough arguments for mapping." << endl;
                cout << "   Configuration file" << endl;
                cout << "   Calibration file" << endl;
                cout << "   Spectrum file" << endl;
                cout << "   Element list for quantification (required, see user manual)" << endl;
                cout << "   Map file (overwritten)" << endl;
                cout << endl;
                return -2008;
            } else {
                string config( argv[2] );
                arguments.configuration_file = config;
                string cal( argv[3] );
                arguments.calibration_file = cal;
                string spec_file( argv[4] );
                arguments.spectrum_file = spec_file;
                string elements( argv[5] );
                arguments.element_list = elements;
                string map_file( argv[6] );
                arguments.map_file = map_file;
            }
            break;
        case COMPARE:
            term_file_index = 6;
            if( argc < term_file_index ) {
                cout << endl;
                cout << "Not enough arguments for comparing measured to calculated." << endl;
                cout << "   Configuration file" << endl;
                cout << "   Standards file" << endl;
                cout << "   Spectrum file" << endl;
                cout << "   Plot file (required, this is the plot output, CSV format)" << endl;
                cout << endl;
                return -2009;
            } else {
                string config( argv[2] );
                arguments.configuration_file = config;
                string stds( argv[3] );
                arguments.standards_file = stds;
                string spec_file( argv[4] );
                arguments.spectrum_file = spec_file;
                string plot_file( argv[5] );
                arguments.plot_file = plot_file;
           }
        case FIT_ONE_STANDARD:
            term_file_index = 6;
            if( argc < term_file_index ) {
                cout << endl;
                cout << "Not enough arguments for fitting one standard." << endl;
                cout << "   Configuration file" << endl;
                cout << "   Standards file" << endl;
                cout << "   Element list" << endl;
                cout << "   Plot file (required, this is the plot output, CSV format)" << endl;
                cout << endl;
                return -2010;
            } else {
                string config( argv[2] );
                arguments.configuration_file = config;
                string stds( argv[3] );
                arguments.standards_file = stds;
                string elements( argv[4] );
                arguments.element_list = elements;
                string plot_file( argv[5] );
                arguments.plot_file = plot_file;
           }
            break;
        case BULK_SUM_MAX:
            term_file_index = 7;
            if( argc < term_file_index ) {
                cout << endl;
                cout << "Not enough arguments for sum." << endl;
                cout << "   Configuration file" << endl;
                cout << "   Calibration file" << endl;
                cout << "   Spectrum file (or CSV file)" << endl;
                cout << "   Element list for quantification (required, see user manual)" << endl;
                cout << "   Plot file (required, calculated spectrum output, CSV format)" << endl;
                cout << endl;
                return -2011;
            } else {
                string config( argv[2] );
                arguments.configuration_file = config;
                string cal( argv[3] );
                arguments.calibration_file = cal;
                string spec_file( argv[4] );
                arguments.spectrum_file = spec_file;
                string elements( argv[5] );
                arguments.element_list = elements;
                string plot_file( argv[6] );
                arguments.plot_file = plot_file;
            }
            break;
        case EM_SDD_DATA:
            term_file_index = 4;
            if( argc < term_file_index ) {
                cout << endl;
                cout << "Not enough arguments for ems sub-command." << endl;
                cout << "   Input file of SEND_SDD_DATA SDF contents (required, CSV format)" << endl;
                cout << "   Output EDR file (required, CSV format)" << endl;
                cout << endl;
                return -2012;
            } else {
                string spec_file( argv[2] );
                arguments.spectrum_file = spec_file;
                string map_file( argv[3] );
                arguments.map_file = map_file;
            }
            break;
        case PRINT_VERSION:
            // don't care about any other params
            return 0;
        case OPTIC_RESPONSE:
            term_file_index = 7;
            if( argc < term_file_index ) {
                cout << endl;
                cout << "Not enough arguments for computing optic response." << endl;
                cout << "   Configuration file" << endl;
                cout << "   Standards file" << endl;
                cout << "   Spectrum file" << endl;
                cout << "   Element list for optic absorption edges and ignored elements in fit (optional but empty string required)" << endl;
                cout << "   Plot file (required, this is the plot output, CSV format)" << endl;
                cout << "See terminal output for optic response curve." << endl;
                cout << endl;
                return -2013;
            } else {
                string config( argv[2] );
                arguments.configuration_file = config;
                string stds( argv[3] );
                arguments.standards_file = stds;
                string spec_file( argv[4] );
                arguments.spectrum_file = spec_file;
                string elements( argv[5] );
                arguments.element_list = elements;
                string plot_file( argv[6] );
                arguments.plot_file = plot_file;
           }
            break;
        default:
            return -2020;
            break;
    }
    int arg_index;
    for( arg_index=term_file_index; arg_index<argc; arg_index++ ) {
        string temp( argv[arg_index] );
        if( temp.length() <= 0 ) continue;
        if( temp.substr( 0, 1 ) == "-" ) {  //  Then its an option
            vector <string> records;
            int result = parse_records( COMMA_CHARACTER, temp, records );
            if( result >= 0 && records.size() > 0 ) {
                if( records[0] == "-e" ) {  //  Energy calibration
                    float temp_eV_start, temp_eV_ch = 0;
                    if( records.size() < 3 ) {
                        result = -1;
                    } else {
                        istringstream temp_stream1( records[1] );
                        temp_stream1 >> temp_eV_start;
                        if( ! temp_stream1 ) result = -1;
                        istringstream temp_stream2( records[2] );
                        temp_stream2 >> temp_eV_ch;
                        if( ! temp_stream2 ) result = -1;
                    }
                    if( result >= 0 && temp_eV_ch > 0 ) {
                        arguments.eV_start = temp_eV_start;
                        arguments.eV_ch = temp_eV_ch;
                    } else {
                        arguments.invalid_arguments += "Invalid energy calibration in argument list: " + temp;
                        return -2024;
                    }
                } else if( records[0] == "-b" || records[0] == "-bh" || records[0] == "-bx" ) {  //  Background controls
                    vector <float> bkg_temp_params;
                    result = 0;
                    string bad_bkg_parameters;
                    //  Parse each in order, any can be left out
                    int i_bkg_arg;
                    for( i_bkg_arg=1; i_bkg_arg<records.size(); i_bkg_arg++ ) {
                        if( records[i_bkg_arg].length() <= 0 ) {
                            bkg_temp_params.push_back( 0 );
                        } else {
                            istringstream temp_stream1( records[i_bkg_arg] );
                            float temp_bkg_value = 0;
                            temp_stream1 >> temp_bkg_value;
                            if( ! temp_stream1 ) {
                                result = -1;
                                bad_bkg_parameters += " " + records[i_bkg_arg];
                            } else {
                                bkg_temp_params.push_back( temp_bkg_value );
                            }
                        }
                    }
                    if( result < 0 ) {
                        arguments.invalid_arguments += "Invalid background parameter in argument list:" + bad_bkg_parameters;
                        return -2027;
                    } else {
                        if( records[0] == "-b" ) arguments.bkg_args = bkg_temp_params;
                        else if( records[0] == "-bh" ) arguments.bh_args = bkg_temp_params;
                        else if( records[0] == "-bx" ) arguments.bx_args = bkg_temp_params;
                    }
                } else if( records[0] == "-T" ) {  //  Control for detector shelf adjustment factor and slope vvs energy
                    vector <float> det_temp_params;
                    result = 0;
                    string bad_det_parameters;
                    //  Parse each in order, any can be left out
                    int i_det_arg;
                    for( i_det_arg=1; i_det_arg<records.size(); i_det_arg++ ) {
                        if( records[i_det_arg].length() <= 0 ) {
                            det_temp_params.push_back( 0 );
                        } else {
                            istringstream temp_stream1( records[i_det_arg] );
                            float temp_det_value = 0;
                            temp_stream1 >> temp_det_value;
                            if( ! temp_stream1 ) {
                                result = -1;
                                bad_det_parameters += " " + records[i_det_arg];
                            } else {
                                det_temp_params.push_back( temp_det_value );
                            }
                        }
                    }
                    if( result < 0 ) {
                        arguments.invalid_arguments += "Invalid detector shelf parameter in argument list:" + bad_det_parameters;
                        return -2027;
                    } else {
                        arguments.detector_shelf_parameters = det_temp_params;
                    }
                } else if( records[0] == "-d" ) {  //  Choose which of multiple detectors to include
                    int temp_det_sel = -1;
                    result = -1;
                    if( records.size() > 1 ) {
                        istringstream temp_stream1( records[1] );
                        temp_stream1 >> temp_det_sel;
                        if( ! temp_stream1 ) {
                            result = -1;
                        } else {
                            arguments.detector_select = temp_det_sel;
                            result = 1;
                        }
                    }
                    if( result < 0 ) {
                        arguments.invalid_arguments += "Invalid detector selection in argument list: " + temp;
                        return -2025;
                    }
                } else if( records[0] == "-m" ) {  //  Maximum number of spectrum files to read for map
                    int temp_max_map = -1;
                    result = -1;
                    if( records.size() > 1 ) {
                        istringstream temp_stream1( records[1] );
                        temp_stream1 >> temp_max_map;
                        if( ! temp_stream1 ) {
                            result = -1;
                        } else {
                            arguments.max_map_arg = temp_max_map;
                            result = 1;
                        }
                    }
                    if( result < 0 ) {
                        arguments.invalid_arguments += "Invalid maximum number of map spectra in argument list: " + temp;
                        return -2026;
                    }
                } else if( records[0] == "-q" ) {  //  Specify outputs to map file
                    result = -1;
                    if( records.size() > 1 ) {
                        arguments.quant_map_outputs = records[1];
                        result = 1;
                    }
                    if( result < 0 ) {
                        arguments.invalid_arguments += "No output selection for map files: " + temp;
                        return -2026;
                    }
                } else if( records[0] == "-f" ) {  //  Turn off adjustments to energy calibration in fits
                    arguments.fit_adjust_energy = false;
                } else if( records[0] == "-g" ) {  //  Turn off adjustments to detector resolution in fits
                    arguments.fit_adjust_width = false;
                } else if( records[0] == "-v" ) {  //  Turn on convolution of Compton components with detector resolution
                    arguments.convolve_Compton = true;
                } else if( records[0] == "-c" ) {  //  Treat some elements as carbonates instead of oxides
                    arguments.carbonates = true;
                } else if(records[0] == "-t") {
                    int tmp = -1;
                    result = -1;
                    if( records.size() > 1 ) {
                        istringstream temp_stream1( records[1] );
                        temp_stream1 >> tmp;
                        if( ! temp_stream1 ) {
                            result = -1;
                        } else {
                            arguments.map_threads = tmp;
                            result = 1;
                        }
                    }
                    if( result < 0 ) {
                        arguments.invalid_arguments += "Invalid thread count in argument list: " + temp;
                        return -2027;
                    }
                } else if( records[0] == "-s" ) {  //  Select standard from input file by number or name
                    int temp_std = -1;
                    result = -1;
                    if( records.size() > 1 ) {
                        istringstream temp_stream1( records[1] );
                        temp_stream1 >> temp_std;
                        if( ! temp_stream1 ) {
                            //  Treat the option as a string
                            arguments.standard_name = records[1];
                            arguments.standard_selected = true;
                            result = 1;
                        } else {
                            arguments.standard_selection = temp_std;
                            arguments.standard_selected = true;
                            result = 1;
                        }
                    }
                    if( result < 0 ) {
                        arguments.invalid_arguments += "Invalid standard selection in argument list: " + temp;
                        return -2028;
                    }
                } else if( records[0] == "-w" ) {  //  Minimum weight in stds file for inclusion in evaluate output
                    float temp_wgt = 0;
                    result = 0;
                    if( records.size() < 1 ) {
                        result = -1;
                    } else {
                        istringstream temp_stream1( records[1] );
                        temp_stream1 >> temp_wgt;
                        if( ! temp_stream1 ) result = -1;
                     }
                    if( result >= 0 ) {
                        arguments.min_wgt_eval = temp_wgt;
                    } else {
                        arguments.invalid_arguments += "Evaluation weight missing or invalid in argument list";
                        return -2030;
                    }
                } else if( records[0] == "-u" ) {  //  Output evaluation file during Calibration or plot file during Evaluate
                    result = -1;
                    if( records.size() > 1 ) {
                        arguments.cal_eval_file = records[1];
                        result = 1;
                    }
                    if( result < 0 ) {
                        arguments.invalid_arguments += "File name missing for -u option";
                        return -2029;
                    }
                } else if(records[0] == "-n") {
                    int tmp = -1;
                    result = -1;
                    if( records.size() > 1 ) {
                        istringstream temp_stream1( records[1] );
                        temp_stream1 >> tmp;
                        if( ! temp_stream1 ) {
                            result = -1;
                        } else {
                            arguments.normalization = tmp;
                            result = 1;
                        }
                    }
                    if( result < 0 ) {
                        arguments.invalid_arguments += "Invalid normalization in argument list: " + temp;
                        return -2030;
                    }
                } else if( records[0] == "-Fe" ) {  //  Iron default oxide ratio
                    float temp_oxide_ratio = -1;
                    result = 0;
                    if( records.size() < 1 ) {
                        result = -1;
                    } else {
                        istringstream temp_stream1( records[1] );
                        temp_stream1 >> temp_oxide_ratio;
                        if( ! temp_stream1 ) result = -1;
                        if( temp_oxide_ratio < 0 ) result = -1;
                     }
                    if( result >= 0 ) {
                        arguments.iron_oxide_ratio = temp_oxide_ratio;
                    } else {
                        arguments.invalid_arguments += "Iron oxide ratio missing or invalid in argument list";
                        return -2031;
                    }
                } else {
                    arguments.invalid_arguments += "Invalid option in argument list: " + temp;
                    return -2023;
                }
            } else {
                arguments.invalid_arguments += "Invalid option format in argument list: " + temp;
                return -2022;
            }
        } else if( arg_index == term_file_index ) { //  Then its the terminal output file
            arguments.terminal_text_file = temp;
        } else {
            arguments.invalid_arguments += "Too many arguments that do not start with a minus sign. " + temp;
            return -2021;
        }
    }

    return 0;
}


//  Figure out command and write list of not found
int parse_command( const string &sub_command, PIQUANT_SUBCOMMAND &cmd ) {

    //  make upper case
    string cmd_uc = upper_trim( sub_command );

    //  Interpret sub-command
    //  enum PIQUANT_subcommand ( ENERGY_CAL, PLOT, PRIMARY, CALCULATE, CALIBRATE, QUANTIFY, EVALUATE, MAP );
    if( cmd_uc.substr(0,4) == "ENE" ) {
        cmd = ENERGY_CAL;
    } else if( cmd_uc.substr(0,3) == "PLO" ) {
        cmd = PLOT;
    } else if( cmd_uc.substr(0,3) == "PRI" ) {
        cmd = PRIMARY;
    } else if( cmd_uc.substr(0,4) == "CALC" ) {
        cmd = CALCULATE;
    } else if( cmd_uc.substr(0,4) == "CALI"
        || ( cmd_uc.length() == 3 && cmd_uc.substr(0,3) == "CAL" ) ) {
        cmd = CALIBRATE;
    } else if( cmd_uc.substr(0,3) == "QUA" ) {
        cmd = QUANTIFY;
    } else if( cmd_uc.substr(0,3) == "EVA" ) {
        cmd = EVALUATE;
    } else if( cmd_uc.substr(0,3) == "MAP" ) {
        cmd = MAP;
    } else if( cmd_uc.substr(0,3) == "COM" ) {
        cmd = COMPARE;
    } else if( cmd_uc.substr(0,3) == "FIT" ) {
        cmd = FIT_ONE_STANDARD;
    } else if( cmd_uc.substr(0,3) == "SUM" ) {
        cmd = BULK_SUM_MAX;
    } else if( cmd_uc.substr(0,3) == "EMS" ) {
        cmd = EM_SDD_DATA;
    } else if( cmd_uc.substr(0,3) == "VER" ) {
        cmd = PRINT_VERSION;
    } else if( cmd_uc.substr(0,3) == "OPT" ) {
        cmd = OPTIC_RESPONSE;
    } else {
        cout << endl;
        cout << "Invalid sub-command; " << cmd_uc << ", possibilities are (only the first 3 letters are checked):" << endl; // What about CALI vs CAL vs CALC?
        cout << "   energy_calibrate - use one or two elements to associate with largest peaks and find energy calibration" << endl;
        cout << "   plot             - plot the spectrum, or a CSV file with appropriate format" << endl;
        cout << "   primary_spectrum - calculate the primary spectrum, with and without optic and filter" << endl;
        cout << "   calculate        - calculate a spectrum of the first standard in the list" << endl;
        cout << "   compare          - compare a measured spectrum to its calculated spectrum" << endl;
        cout << "   optic            - compute an optic response curve using a measured spectrum and its known composition" << endl;
        cout << "   calibrate        - perform quantitative calibration using standards and write the calibration file (can be just cal)" << endl;
        cout << "   quantify         - use the calibration file and element list to fit and quantify a spectrum" << endl;
        cout << "   evaluate         - quantify each standard as an unknown and check against the known values" << endl;
        cout << "   map              - quantify a set of spectra and write a map file" << endl;
        cout << "   sum              - calculate sum and maximum value spectra from a set of spectra" << endl;
        cout << "   ems              - convert output of SEND_SDD_DATA command (SDF contents in csv file) to EDR (csv) format" << endl;
        cout << "   version          - print piquant version" << endl;
        cout << endl;
        return -2000;
    }
    return 0;
}
