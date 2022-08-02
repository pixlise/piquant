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

#include <algorithm>
#include "quantWriteCalibrationTXT.h"
#include "XRFconstants.h"
#include "XRFcontrols.h"
#include "upper_trim.h"


using namespace std;

//      Write and read calibration files with elements and calibration factors
//      in the old text format with just the calibration factors

//  Written Apr. 12, 2017
//  Modified June 27, 2017 to initialize ECFs to zero in output vector from quantReadCalibrationTXT
//  Modified Dec. 6, 2017 to fix comments to match new terminology
//  Modified July 31, 2018
//      Re-write to match quantWriteCalibrationCSV (use element list entries loaded in quantWriteResults)
//      Use weights from standards input if CSV standards file was read
//      Use acceptance criteria based on min % and max rsd if old TXT standards file was read
//      Add header with date and time


int quantWriteCalibrationTXT( vector <StandardInformation> &standards, const string &calibrationFileName,
                        std::string date_and_time ) {

    //  Use the information in vector <StandardInformation> standards
    //  to find the element calibration factors for the input list of elements
    //  The calibration factor is the weighted average of coefficients
    //  Code and algorithm taken from unkQuant.cpp of March 24, 2015
    //  Via loadStandardQuant.cpp (calculateCalibration function) of April 1, 2015

    //  Check to be sure some standards have been loaded
    if( standards.size() <= 0 ) return -410;
    //  First generate a complete element list
    vector <Element> stdElements;
    int is;
    for( is=0; is<standards.size(); is++ ) {
        const vector <ElementListEntry> &std_element_list = standards[is].element_list;
        int ise;
        for ( ise=0; ise<std_element_list.size(); ise++ ) {
            //  Use acceptance criteria based on min % and max rsd if old TXT standards file was read
            if( !standards[is].user_weights ) {
                //  Skip low-Z elements since we can't calibrate them
                if( std_element_list[ise].element.Z() < calibration_minimum_z ) continue;
            }
            bool found = false;
            int ie;
            for ( ie=0; ie<stdElements.size(); ie++ ) {
                if ( stdElements[ie] == std_element_list[ise].element ) {
                    found = true;
                    continue;
                };
			};
            if ( ! found ) {
                stdElements.push_back( std_element_list[ise].element );
            };
		};
	};
    //  Sort the element list in increasing order by atomic number
    sort( stdElements.begin(), stdElements.end() );

    //  Now compute averages over the list of standards for each element
    const int nElements = stdElements.size();
    vector <float> avgFitCoeff( nElements, 0 );
    vector <float> avgWeight( nElements, 0 );
    vector <float> avgRSD( nElements, 0 );
    vector <int> stdsCount( nElements, 0 );
    float overallFitCoeff = 0;
    float overallWeight = 0;
    int ie;
    //  Calculate average fit coefficients for each element
    for ( ie=0; ie<nElements; ie++ ) {
        vector <float> indivFitCoeff;
        vector <float> indivWeight;
        cout << "Calculating ECF for element " << stdElements[ie].symbol() << endl;
        int is;
        for( is=0; is<standards.size(); is++ ) {
            vector <ElementListEntry> &std_element_list = standards[is].element_list;
            int i;
            for ( i=0; i<std_element_list.size(); i++ ) {
                float coeff = std_element_list[i].coefficient;
                //float fitErrorPct = std_element_list[i].rel_err_coeff;
                float weight = std_element_list[i].weight;
                if ( std_element_list[i].element == stdElements[ie] ) {
                    cout << "      Standard: ";
                    if( standards[is].names.size() > 0 ) cout << standards[is].names[0];
                    stdsCount[ie] ++;
                    if( coeff <= 0 ) {
                        cout << "  coefficient is zero or negative" << endl;
                        continue;
                    }
                    if( !standards[is].user_weights ) {
                    //  Use acceptance criteria based on min % and max rsd if old TXT standards file was read
                        if( std_element_list[i].element.Z() < calibration_minimum_z ) weight = 0;
                        if ( std_element_list[i].percent < calibration_minimum_fraction * 100
                            || std_element_list[i].rel_err_coeff > calibration_maximum_rsd * 100 ) weight = 0;
                        if( weight <= 0 ) cout << "  failed minimum fraction or maximum RSD, ";
                    }
                    if( weight <= 0 ) {
                        cout << "  not included, weight is zero" << endl;
                        continue;
                    }
                    avgFitCoeff[ie] += coeff * weight;
                    avgWeight[ie] += weight;
                    indivFitCoeff.push_back( coeff );
                    indivWeight.push_back( weight );
                    cout << "   % " << std_element_list[i].percent << "  wgt " << weight << "   coeff " << coeff << endl;
                }
            };
        };
        if ( avgWeight[ie] > 0 ) {
            overallFitCoeff += avgFitCoeff[ie];
            overallWeight += avgWeight[ie];
            avgFitCoeff[ie] /= avgWeight[ie];
            cout << "      Final ECF for " << stdElements[ie].symbol() << "   " << avgFitCoeff[ie];
            cout << "    " << stdsCount[ie] << " standards, total weight " << avgWeight[ie];
            cout << endl;
        } else {
            cout << "      Total weight was zero for " << stdElements[ie].symbol() << endl;
            continue;
        };
        //			calculate the relative standard deviation over all of the samples
        if ( indivFitCoeff.size() > 0 ) {
            int i;
            for ( i=0; i<indivFitCoeff.size(); i++ ) {
                float diff = indivFitCoeff[i] - avgFitCoeff[ie];
                avgRSD[ie] += indivWeight[i] * diff * diff;
            };
            avgRSD[ie] = sqrt ( avgRSD[ie] / avgWeight[ie] ) / avgFitCoeff[ie];
        };
//        cout << "avgRSD " << stdElements[ie].symbol() << "   " << avgRSD[ie] << endl;
    };
    //		use the overall average if no fit coefficient available from standards
    if ( overallWeight > 0 ) {
        overallFitCoeff /= overallWeight;
    } else {
        overallFitCoeff = 1;
    };
    for ( ie=0; ie<nElements; ie++ ) if ( avgWeight[ie] <= 0 || avgFitCoeff[ie] <= 0 ) {
        avgFitCoeff[ie] = overallFitCoeff;
        cout << "Using average ECF for element " << stdElements[ie].symbol() << "   " << avgFitCoeff[ie] << endl;
    }

    //  Now write the element list and calibration factors to the calibration file

    ofstream outFile ( calibrationFileName.c_str() );
    if( ! outFile ) {
        cout << "Couldn't open calibration file, file name " << calibrationFileName << endl;
        return -1;
    }
	outFile.setf( ios::fixed, ios::floatfield );
	outFile.precision(4);
    outFile << "0 PIQUANT, Text Calibration File,     written, " << date_and_time << endl;
    outFile << "0 Lines that start with zero are skipped." << endl;
    //	Now write element list, with fake DLL version for this calibration at the end
    const int element_width = 10;
    int versionDLL = 299;
    float avg_Compton_cps = 0;
    float avg_Rayleigh_cps = 0;

    //  Write the number of elements and the element atomic numbers
    outFile << stdElements.size();
    for ( ie=0; ie<stdElements.size(); ie++ ) {
        outFile << setw(element_width) << stdElements[ie].Z();
    }
    //  Add the DLL version (a legacy placeholder) to the end of the element list
    outFile << "  " << versionDLL;
    outFile << endl;
    //	Write the calibration factor for each element in the list
    for ( ie=0; ie<stdElements.size(); ie++ ) {
        outFile << "  " << avgFitCoeff[ie];
    }
    //	Add Compton and Rayleigh info (counts per second per milliAmp)
    //  For now these are zeros as placeholders, to be replaced once C and R fits are working
    outFile << "  " << avg_Compton_cps << "  " << avg_Rayleigh_cps;
    outFile << endl;
    outFile.close();

    return 0;

};


int quantReadCalibrationTXT( const std::string &calibrationFileName,
        std::vector <Element> &cal_element_list, std::vector <float> &cal_factor_list, std::ostream &logger ) {
	//  Read the old-style text calibration file and return the element calibration factors
    float avg_Compton_read = 1;
    float avg_Rayleigh_read = 1;
    int cal_version = 0;
    int ne_in = 0;
    //		open calibration file
	ifstream CalFile(calibrationFileName.c_str(), ios::in);
	if ( !CalFile ) {
		logger << "Can't read calibration file, using unity calibration factors" << endl;
		return -1;
	} else {
	    logger << "Reading calibration file from " << calibrationFileName << endl;
	};
    while ( ! ( ! CalFile ) ) {
        CalFile >> ne_in;
		if ( !CalFile ) break;
        //			skip this line if no entries (can be comment)
		if ( ne_in <= 0 ) {
			string skip;
			getline( CalFile, skip );
			continue;
        } else {
            int i;
            for ( i=0; i<ne_in; i++ ) {
                int element_Z;
                CalFile >> element_Z;
                if( Element::check_Z( element_Z ) ) {
                    Element elTrial( element_Z );
                    cal_element_list.push_back( elTrial );
                } else {
                    logger << "Invalid element Z " << element_Z << endl;
                };
            }
            CalFile >> cal_version;
            cal_factor_list.resize( ne_in, 0 );
            for ( i=0; i<ne_in; i++ ) {
                float element_factor;
                CalFile >> element_factor;
                cal_factor_list[i] = element_factor;
            }
            CalFile >> avg_Compton_read;
            CalFile >> avg_Rayleigh_read;
            break;
        }
    }
    CalFile.close();
    return cal_element_list.size();
};

