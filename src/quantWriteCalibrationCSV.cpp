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

#include <iostream>
#include <fstream>
#include <iomanip>
#include "quantWriteCalibrationCSV.h"
#include "XRFconstants.h"
#include "XRFcontrols.h"


using namespace std;

//      Write calibration file in new CSV format
//          that matches the new standards input format

//  Written Dec. 9, 2017
//  Modified Jan. 26, 2018
//      Add net peak intensity to element list and calibration file output
//      Don't write irrelevant info for Ignore, Exclude, and Matrix element list entries
//  Modified July 31, 2018
//      Use weights from standards input if CSV standards file was read
//      Use acceptance criteria based on min % and max rsd if old TXT standards file was read
//  Modified Nov. 4, 2019
//      Add atomic number to calibration file as an extra column (to aid in sorting by element)



int quantWriteCalibrationCSV( std::vector <StandardInformation> &standards,
        std::string &calibrationFileName, std::string date_and_time ) {

    //  Use the information in vector <StandardInformation> standards
    //  to find the element calibration factors for the input list of elements
    //  for each standard and put them in the element list
    //  The information for each standard is written separately so
    //  ECFs can be calculated based on different standards for each
    //  unknown if desired (not implemented yet).
    //  See setupStandardsCSV.cpp for format description

    //  Check to be sure some standards have been loaded
    if( standards.size() <= 0 ) return -410;

    //  Use the information in the element list, which was loaded during quantWriteResults
    //      with the information from the correct spectrum fit component
    //  Only the ECF info needs to be loaded here
    int is;
    for( is=0; is<standards.size(); is++ ) {
        int ie;
        for( ie=0; ie<standards[is].element_list.size(); ie++ ) {
            //  Set EFC to zero to avoid old values being written
            standards[is].element_list[ie].ecf = 0;
            standards[is].element_list[ie].ecf_sigma = 0;
            //  Skip if qualifier is I, X, or M
            if( standards[is].element_list[ie].qualifier == IGNORE
                    || standards[is].element_list[ie].qualifier == EXCLUDE
                    || standards[is].element_list[ie].qualifier == MATRIX ) continue;
            float coeff = standards[is].element_list[ie].coefficient;
            if( coeff <= 0 || coeff == COEFFICIENT_NO_COMPONENT ) continue;    //  No coefficient for this component, skip this element list entry
            if( ! standards[is].user_weights ) {
                //  Use acceptance criteria based on min % and max rsd if old TXT standards file was read
               if ( standards[is].element_list[ie].percent < calibration_minimum_fraction * 100
                    || standards[is].element_list[ie].rel_err_coeff > calibration_maximum_rsd * 100 ) {
                        standards[is].element_list[ie].weight = 0;
                }
            }
            standards[is].element_list[ie].ecf = coeff;
            standards[is].element_list[ie].ecf_sigma = standards[is].element_list[ie].rel_err_coeff;
        }
    }

    //  Now write everything to the calibration file

    ofstream outFile ( calibrationFileName.c_str(), ios::out );
    if( ! outFile ) {
        cout << "Couldn't open calibration file, file name " << calibrationFileName << endl;
        return -1;
    }
	outFile.setf( ios::fixed, ios::floatfield );
	outFile.precision(4);
	//  Write headers, first one used to identify file type
    outFile << "PIQUANT, Calibration File,     written, " << date_and_time << endl;
    outFile << "Element, Emission line, Fit qualifier, Type, Percent, ";
    outFile << "Uncertainty,  Oxide ratio, Weight, ECF, ECF Sigma, Intensity, Atomic number" << endl;

    //  Loop over list of standards and write an entry for each one
    for( is=0; is<standards.size(); is++ ) {
        //  Write out any comments that preceded the standard keyword
        int in;
        for( in=0; in<standards[is].preceding_comments.size(); in++ ) {
            outFile << "COMMENT";
            outFile << ", " << standards[is].preceding_comments[in];
            outFile << endl;
        }
        //  Standard keyword and name(s) of standard
        outFile << "STANDARD";
        for( in=0; in<standards[is].names.size(); in++ ) {
            outFile << ", " << standards[is].names[in];
        }
        outFile << endl;
        //  Write all of the comments associated with this standard
        for( in=0; in<standards[is].comments.size(); in++ ) {
            outFile << "COMMENT";
            outFile << ", " << standards[is].comments[in];
            outFile << endl;
        }
        //  Write keywords for carbonates, thickness, and density
        if( standards[is].carbonates ) outFile << "Carbonates" << endl;
        if( standards[is].mat.thickness() > 0 ) {
            outFile << "Thickness";
            outFile << ", " << standards[is].mat.thickness();
            outFile << endl;
        }
        if( standards[is].mat.density() > 0 ) {
            outFile << "Density";
            outFile << ", " << standards[is].mat.density();
            outFile << endl;
        }
        //  Write a line for each entry in the element list and all associated information
        //  At this point all of the relevant information has been copied into this list
        int ie;
        for( ie=0; ie<standards[is].element_list.size(); ie++ ) {
            outFile << standards[is].element_list[ie].element.symbol();
            outFile << ", ";
            switch( standards[is].element_list[ie].quant_level ) {
                case K_LEVEL:   outFile << "K";    break;
                case L_LEVEL:   outFile << "L";    break;
                case M_LEVEL:   outFile << "M";    break;
                case N_LEVEL:   outFile << "N";    break;
                default:                            break;
            };
            outFile << ", ";
            switch( standards[is].element_list[ie].qualifier ) {
                case IGNORE:    outFile << "I";    break;
                case FORCE:     outFile << "F";    break;
                case EXCLUDE:   outFile << "X";    break;
                case MATRIX:    outFile << "M";    break;
                default:                            break;
            };
            if( standards[is].element_list[ie].qualifier == IGNORE ) {
                outFile << endl;
                continue;
            }
            outFile << ", ";
            switch( standards[is].element_list[ie].type ) {
                case NO_COMPONENT:      outFile << "Nothing";    break;
                case ELEMENT:           outFile << "el";    break;
                case COMPTON:           outFile << "inc";    break;   //  Compton scatter = incoherent scatter
                case RAYLEIGH:          outFile << "coh";    break;   //  Rayleigh scatter = incoherent scatter
                case SNIP_BKG:        outFile << "SNIP_bkg";    break;    //  Get rid of element symbol and edge level, do not apply to this component
                case CONTINUUM:        outFile << "Cont_bkg";    break;    //  Get rid of element symbol and edge level, do not apply to this component
                case PRIMARY_LINES:     outFile << "pri";    break;    //  Primary spectrum from anode element lines
                case PRIMARY_CONTINUUM: outFile << "continuum";    break;    //  Get rid of element symbol and edge level, do not apply to this component
                case La:                outFile << "La";    break;
                case Lb1:               outFile << "Lb1";    break;
                case OPTIC_TRANS:       outFile << "Optic";    break;
                case DETECTOR_CE:       outFile << "ComptonEscape";    break;
            };
            outFile.precision(4);
            outFile << ", ";
            if( standards[is].element_list[ie].percent >= 0 ) outFile << standards[is].element_list[ie].percent << "%";
            outFile.precision(1);
            outFile << ", " << standards[is].element_list[ie].uncertainty << "%";
            outFile.precision(1);
            outFile << ", " << standards[is].mat.stoichiometry( standards[is].element_list[ie].element ).formula_ratio;
            if( standards[is].element_list[ie].qualifier == EXCLUDE
                    || standards[is].element_list[ie].qualifier == MATRIX ) {
                outFile << endl;
                continue;
            }
            outFile.precision(2);
            outFile << ",   " << standards[is].element_list[ie].weight;
            outFile.precision(4);
            outFile << ", " << standards[is].element_list[ie].ecf;
            outFile.precision(1);
            outFile << ", " << standards[is].element_list[ie].ecf_sigma << "%";
            //  Add the peak intensity to the element list
            outFile << ", " << standards[is].element_list[ie].intensity;
            //  Add the atomic number to help in sorting
            outFile << ", " << standards[is].element_list[ie].element.Z();
            outFile << endl;
        }
        //  Spectrum keyword and name of spectrum file (end of entry)
        outFile << "SPECTRUM";
        outFile << ", " << standards[is].spectrumFileName;
        outFile << endl;
        //  Separate entries with a blank line
        outFile << endl;
    }

    outFile.close();

    return 0;

};

