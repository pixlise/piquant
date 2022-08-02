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

#include <exception>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include "upper_trim.h"
#include "parse_arguments.h"
#include "parse_records.h"
#include "XRFconstants.h"
#include "XRFcontrols.h"
#include "XRFutilities.h"
#include "setupStandardsCSV.h"

//  This function reads both standards input and calibration files
//  The standards input file contains information needed to process a standard
//  The calibration file contains the information about the standard and
//      the results of processing the standard spectrum
//  The files are in comma-separated value format with the following entries:
//      Each line either starts with one of the keywords or is an element line
//      The keywords (not case sensitive) are comment, standard, and spectrum
//          The comment indicates comments that will be copied to the calibration file
//          The standard keyword clears all info and is followed by any names for the standard
//          The spectrum keyword is followed by the name of a spectrum file.  It causes
//              the spectrum to be processed and a calibration file entry produced.
//              More than one spectrum keyword can be included for any standard
//      If the line does not start with a keyword, then it is an element entry
//          An element entry has several fields in a required order to allow easy parsing and archiving
//          The entries are separated by commas and are described in order:
//              Element symbol (case sensitive)
//              Element emission line description (currently K, L, M or N)
//              Element processing qualifier (optional, but must have its own comma-separated field, i.e ,, or , ,)
//                  I - Fit the element emission line (or all lines if none specified) in the spectrum but ignore the results
//                  X - Exclude the element from the spectrum and the standard composition
//                  F - Forces element to be included in the spectrum fit and composition (for future use)
//                  M - The element is part of the matrix (for future use)
//              Spectrum component type, if no component is given Element is assumed
//                  Element emission line
//                  Coherent (Rayleigh) scatter line
//                  Incoherent (Compton) scatter line
//                  Background
//              Composition - fraction or percent of this element in the standard
//                  Can optionally be followed by %, f, or ppm.  If none specified, percent assumed
//              Uncertainty in composition, relative, in the same units as the composition (percent of fraction, not ppm)
//              Oxide ratio of this element (negative means use default oxide ratio for this element)
//              Weight to be given to this standard when computing the element correction factor for this element
//              (The following information is for calibration files only and is ignored for standards input files)
//              Element calibration factor resulting from fitting the element emission line to the spectrum
//              Relative error of the calibration factor from the fitting function as a percent
//

//  Written Nov. 1, 2017        (based on setupStandardsTXT)
//  Nov. 3, 2017    Finished initial testing
//  Modified Dec. 6, 2017
//      Move add element list entry to new function in parse_element_list
//      Remove handling of oxide ratio -1 (XrayMaterial already handles it properly)
//      Add ecf and ecf_sigma so it can handle calibration files
//  Modified Dec. 8, 2017
//      Added above comments and spectrum component information
//      Don't use underscore to separate element, emission line, qualifier, and component type,
//          put them into separate comma-separated fields.
//      Don't include element list entries with IGNORE qualifier or COMPTON or RAYLEIGH
//          component type into element list for standard composition
//  Modified Dec. 10, 2017
//      Skip over headers when reading calibration file
//      Pick up comments that are not between standard and spectrum keywords
//  Modified Dec. 23, 2017
//      Get rid of termOutFile, use cout (as redirected if appropriate)
//  Modified Jan. 3, 2017
//      Fix qualifier M, string was "X" changed to "M" in translation if block
//      Move tab, single quote, double quote, blank, comma, and underscore definitions to XRFconstants.h
//  Modified Mar. 7, 2018
//      Use utility function to extract file path for spectrum file
//      Use utility to check file extensions
//  Modified July 31, 2018
//      Indicate that weights were read in and specified by the user (not true for TXT standards input files)
//  Modified Sep. 21, 2018
//      Change handling of suffix on numbers to work on Mac (strip suffix before reading with stringstream)
//  Modified May 6, 2018
//      strip_suffix was moved to XRFutilities
//  Modified June 6, 2019
//      Handle EXCLUDE and MATRIX elements correctly (for standards and evaluate)
//      Implement negative oxide ratio entry means use default oxide ratio (from XrayMaterial.default_oxide_ratio)
//  Modified Dec. 16, 2019
//      Add element qualifier OUTPUT ("O") to force the element to be included in the evaluate list (with zeros if not in any standard in this run)
//      Elements in input list must be ignored here if they have that qualifier
//      Re-arrange how elements are added to only add if qualifier matches (NO_QUALIFIER or FORCE or MATRIX to add to std)
//  Modified Nov. 9, 2020
//      Process standards without SPECTRUM keyword, even if there is an end of file just after the input lines
//  Modified Nov. 24, 2020
//      Put element percent into element list as given percent
//  Modified Dec. 15, 2020
//      Added check for "CO3" in element list to indicate some elements to be included as carbonates instead of oxides
//  Modified Feb. 1, 202`
//      Fix bug reading standards files with only commas on blank lines (between standards)
//  Modified Apr. 8, 2021   Ignore blank standard names (zero length or a single blank), fix multiple SPECTRUM keywords for the same standard


using namespace std;

int setupStandardsCSV( const std::string &standardsInputFileName,
                      std::vector <StandardInformation> &standards_out, const float miniumum_amount ) {

    standards_out.clear();

    //  Get path from input file name (argument 1)
    string standardsPathName;
    string standardsFileOnly;
    /*bool path_found =*/ extract_path( standardsInputFileName, standardsPathName, standardsFileOnly );

//		open file of standard compositions and spectrum file names
	string StdCalEntryName = standardsInputFileName;
	cout << "Reading standard information from file " << StdCalEntryName << endl;
    if( miniumum_amount != 0 ) cout << "Ignoring elements less than " << 10000 * miniumum_amount << " ppm" << endl;
	ifstream StdCalEntryFile(StdCalEntryName.c_str(), ios::in);
	if ( !StdCalEntryFile ) {
		cout << "Cannot open file " << StdCalEntryName << endl;
		return -1;
	};
	bool error = false;
    StandardInformation standard_entry;
    //  Indicate that weights were read in as specified by the user (not true for TXT standards input files)
    standard_entry.user_weights = true;
    standard_entry.carbonates = false;  //  Can be set by CARBONATES keyword for each standard
    float thickness_save = 0;
    float density_save = 0;
    int line_number = 0;
    int error_number = 0;
    //  Hold on to any comments that are read before the first standard keyword
    //      or between spectrum ans standard keywords
    vector <string> stray_comments;
    bool between_standards = true;
//		read and process lines in standard compositions file
	while ( ! ( !StdCalEntryFile ) ) {
        bool entry_error = false;
        string input_str;
        getline( StdCalEntryFile, input_str );
        line_number++;
		//  Skip empty line (unless this is the end of the file, in which case we may need to process the last standard)
        if( ! ( !StdCalEntryFile ) && input_str.size() <= 0 ) continue;
		//  Get rid of trailing CR if file is Windows line endings on Linux or Mac
		if( (int) (input_str.data()[input_str.length()-1]) == 13 ) input_str.erase( input_str.length()-1,1);
//		cout << "setupStandardsCSV   in " << line_number << "  " << input_str.size() << " : " << input_str << ":  ! " << !StdCalEntryFile << endl;
        //  Parse line into comma separated fields
        vector <string> records;
        int result = parse_records( COMMA_CHARACTER, input_str, records );
        if( result < 0 ) {
            cout << "*** Error parsing comma separated entries on line " << line_number << ". ***" << endl;
            entry_error = true;
            error = true;
            continue;
        }
        //  Ignore line is first entry is zero length and between standards (line with only commas inserted by Excel)
        if( records.size() == 0 || ( between_standards && records[0].length() == 0 ) ) continue;

        //  Check for keywords first, if not assume element entry
        //  Note that STANDARD keyword must be handled after the processing at the end of this block if

        string keyword;
        if( records.size() > 0 ) keyword = upper_trim( records[0] );
		//  If this is a calibration file, skip both header lines
		if( keyword == "PIQUANT" ) continue;
		else if( keyword == "ELEMENT" ) continue;
        //  CARBONATES keyword
		else if( keyword == "CARBONATES" ) {
            standard_entry.carbonates = true;
            continue;
        //  FRACTIONS keyword
		} else if( keyword == "FRACTIONS" ) {
            if( records.size() > 1 && ( upper_trim( records[1] ) =="FORMULA" || upper_trim( records[1] ) =="OXIDE" ) ) {
                standard_entry.input_fractions_are_formula = true;
            } else standard_entry.input_fractions_are_formula = false;
            continue;
        //  THICKNESS keyword
		} else if( keyword == "THICKNESS" ) {
            float thickness_in = 0;
            if( records.size() > 1 ) {
                istringstream value( records[1] );
                value >> thickness_in;
                if ( ! value || thickness_in < 0 ) {
                    cout << "Invalid thickness on line " << line_number << ",  value " << ", " << records[1] << endl;
                    entry_error = true;
                    error = true;
                } else if( thickness_in > 0 ) {
                    thickness_save = thickness_in;
                }
            }
            continue;
        //  DENSITY keyword
		} else if( keyword == "DENSITY" ) {
            float density_in = 0;
            if( records.size() > 1 ) {
                istringstream value( records[1] );
                value >> density_in;
                if ( ! value || density_in < 0 ) {
                    cout << "Invalid thickness on line " << line_number << ",  value " << ", " << records[1] << endl;
                    entry_error = true;
                    error = true;
                } else if( density_in > 0 ) {
                    density_save = density_in;
                }
            }
            continue;
        //  COMMENT keyword
        } else if( keyword == "COMMENT" ) {
            //			save in comments list
            if( records.size() > 1 ) {
                standard_entry.comments.push_back( records[1] );
                //  If this is between the standard and spectrum keywords, capture it for the next standard entry
                if( between_standards ) stray_comments.push_back( records[1] );
                continue;
            }
        //  SPECTRUM keyword
        } else if( keyword == "SPECTRUM" ) {
            if( records.size() < 2 || records[1].size() < 1 || records[1] == BLANK_CHARACTER ) {
                cout << "*** No spectrum file name entry on line " << line_number << ". ***" << endl;
                entry_error = true;
                error = true;
                continue;
            }
            //  Get path from standards input file name if no path specified in spectrum file name
            int slash_pos = records[1].rfind( SLASH_CHARACTER );  //  Use \ for Windows, / for Linux and Mac
            if ( slash_pos <= 0 || slash_pos >= records[1].length()-1 ) slash_pos = records[1].rfind( BACKSLASH_CHARACTER );
            if ( slash_pos <= 0 || slash_pos >= records[1].length()-1 ) slash_pos = -1;
            string spectrumPathName = records[1];
            if( slash_pos < 0 ) spectrumPathName = standardsPathName + records[1];
            standard_entry.spectrumFileName = spectrumPathName;
            between_standards = false;

        //  If no keywords, assume an element information line with ordered entries

        } else if( records.size() > 0 && keyword != "STANDARD" ) {
			int entry_number = 0;
            //  Check for an element information line
            ElementListEntry el_entry;
            //  Carry over information about entry percents
            el_entry.stoichiometry.input_fractions_are_formula = standard_entry.input_fractions_are_formula;
            //  Element symbol (maybe with line and qualifier)
            bool element_error = parse_element_string( records[0], el_entry );
			if( element_error ) {
				cout << "Invalid element symbol or qualifier on line " << line_number << endl;
				entry_error = true;
                error = true;
				continue;
			};
			//  Emission line (KLMN)
			entry_number++;
            if( records.size() > entry_number && records[entry_number].length() > 0 ) {
                string em_line = upper_trim( records[entry_number].substr( 0, 1 ) );
                if( em_line == "K" ) el_entry.quant_level = K_LEVEL;
                else if( em_line == "L" ) el_entry.quant_level = L_LEVEL;
                else if( em_line == "M" ) el_entry.quant_level = M_LEVEL;
                else if( em_line == "N" ) el_entry.quant_level = N_LEVEL;
                else if( em_line != BLANK_CHARACTER ) {
                    cout << "Invalid emission line symbol on line " << line_number;
                    cout << ",  Element " << el_entry.element.symbol() << ", " << records[entry_number] << endl;
                    entry_error = true;
                }
            }
            //  Element qualifier (IXFM)
			entry_number++;
            if( records.size() > entry_number && records[entry_number].length() > 0 ) {
                string el_qual = upper_trim( records[entry_number].substr( 0, 1 ) );
                if( el_qual == "I" ) el_entry.qualifier = IGNORE;
                else if( el_qual == "F" ) el_entry.qualifier = FORCE;
                else if( el_qual == "X" ) el_entry.qualifier = EXCLUDE;
                else if( el_qual == "M" ) el_entry.qualifier = MATRIX;
                else if( el_qual != BLANK_CHARACTER ) {
                    cout << "Invalid element qualifier symbol on line " << line_number;
                    cout << ",  Element " << el_entry.element.symbol() << ", " << records[entry_number] << endl;
                    entry_error = true;
                }
            }
            //  Spectrum component symbol (optional, especially useful for scattered components)
			entry_number++;
			el_entry.type = ELEMENT;    //  Default for this input line
            if( records.size() > entry_number && records[entry_number].length() > 0 ) {
                string s_comp = upper_trim( records[entry_number].substr( 0, 3 ) );
                if( s_comp.substr( 0, 2 ) == "EL" ) el_entry.type = ELEMENT;
                else if( s_comp == "COM" ) el_entry.type = COMPTON;
                else if( s_comp == "INC" ) el_entry.type = COMPTON;
                else if( s_comp == "RAY" ) el_entry.type = RAYLEIGH;
                else if( s_comp == "COH" ) el_entry.type = RAYLEIGH;
                else if( s_comp != BLANK_CHARACTER ) {
                    cout << "Invalid spectrum component symbol on line " << line_number;
                    cout << ",  Element " << el_entry.element.symbol() << ", " << records[entry_number] << endl;
                    entry_error = true;
                }
            }
            //  Composition (percent or fraction or parts-per-million)
			entry_number++;
			enum ElementEntryNormalization { PPM = -1, FRACTION, PERCENT };
			ElementEntryNormalization composition_format = PERCENT;
            float percent_in = 0;
            if( records.size() > entry_number && records[entry_number].length() > 0 ) {
                istringstream value( strip_suffix(records[entry_number]) );
                value >> percent_in;
                if( 0 <= records[entry_number].rfind( "%" ) && records[entry_number].rfind( "%" ) < records[entry_number].length() ) {
                    composition_format = PERCENT;
                } else if( 0 <= records[entry_number].rfind( "f" ) && records[entry_number].rfind( "f" ) < records[entry_number].length() ) {
                    composition_format = FRACTION;
                } else if( 0 <= records[entry_number].rfind( "F" ) && records[entry_number].rfind( "F" ) < records[entry_number].length() ) {
                    composition_format = FRACTION;
                } else if( 0 <= records[entry_number].rfind( "p" ) && records[entry_number].rfind( "p" ) < records[entry_number].length() ) {
                    composition_format = PPM;
                } else if( 0 <= records[entry_number].rfind( "P" ) && records[entry_number].rfind( "P" ) < records[entry_number].length() ) {
                    composition_format = PPM;
                }
                if( composition_format == PPM ) percent_in *= PPM_PERCENT;
                else if( composition_format == FRACTION ) percent_in *= 100;
                if ( ! value || percent_in < 0 || percent_in > 100 ) {
                    cout << "Invalid composition on line " << line_number << ",  Element " << el_entry.element.symbol() << ", " << records[entry_number] << endl;
                    entry_error = true;
                } else {
                    el_entry.percent = percent_in;
                    el_entry.given = percent_in;
                }
            }
            //  Composition uncertainty, assumed relative (percent or fraction same as composition entry above) unless followed by a or A, then absolute uncertainty
            entry_number++;
            if( records.size() > entry_number && records[entry_number].length() > 0 ) {
                float u = 0;
                istringstream value( strip_suffix(records[entry_number]) );
                value >> u;
                if( composition_format == PPM ) u *= PPM_PERCENT;
                else if( composition_format == FRACTION ) u *= 100;
                if( ( 0 <= records[entry_number].rfind( "a" ) && records[entry_number].rfind( "a" ) < records[entry_number].length() )
                        || ( 0 <= records[entry_number].rfind( "A" ) && records[entry_number].rfind( "A" ) < records[entry_number].length() ) ) {
                    if( percent_in > 0 ) u = ( u / percent_in) * 100;    //  Convert from absolute uncertainty to relative percent
                }
                if ( ! value || u < 0 || u > 100 ) {
                    cout << "Invalid uncertainty on line " << line_number << ",  Element " << el_entry.element.symbol() << ", " << records[entry_number] << endl;
                    entry_error = true;
                } else {
                    el_entry.uncertainty = u;
                }
            }
            //  Oxide ratio (or carbonate ratio)
            entry_number++;
            if( records.size() > entry_number && records[entry_number].length() > 0 ) {
                float oxr;
                istringstream value( records[entry_number] );
                value >> oxr;
                if ( ! value ) {
                    cout << "Invalid oxide ratio on line " << line_number << ",  Element " << el_entry.element.symbol() << ", " << records[entry_number] << endl;
                    entry_error = true;
                } else {
                    //  Negative entry means use default value
                    if( oxr < 0 ) {
                        if( standard_entry.carbonates ) {
                            el_entry.stoichiometry.formula = CARBONATE;
                            el_entry.stoichiometry.formula_ratio = standard_entry.mat.default_carbonate_ratio( el_entry.element );
                            if( el_entry.stoichiometry.formula_ratio == 0 ) { //  This element won't form carbonates, so treat it as an oxide
                                el_entry.stoichiometry.formula_ratio = standard_entry.mat.default_oxide_ratio( el_entry.element );
                                if( el_entry.stoichiometry.formula_ratio ) el_entry.stoichiometry.formula = OXIDE;
                            }
                        } else {
                            el_entry.stoichiometry.formula_ratio = standard_entry.mat.default_oxide_ratio( el_entry.element );
                            if( el_entry.stoichiometry.formula_ratio ) el_entry.stoichiometry.formula = OXIDE;
                        }
                    } else {
                        el_entry.stoichiometry.formula_ratio = oxr;
                        if( oxr > 0 ) {
                            if( standard_entry.carbonates && standard_entry.mat.default_carbonate_ratio( el_entry.element ) > 0 ) el_entry.stoichiometry.formula = CARBONATE;
                            else el_entry.stoichiometry.formula = OXIDE;
                        }
                    }
                }
            }
            //  Weight of element correction factor from this standard for this element
            entry_number++;
            if( records.size() > entry_number && records[entry_number].length() > 0 ) {
                float w;
                istringstream value( records[entry_number] );
                value >> w;
                if ( ! value || w < 0 ) {
                    cout << "Invalid weight on line " << line_number << ",  Element " << el_entry.element.symbol() << ", " << records[entry_number] << endl;
                    entry_error = true;
                } else {
                    el_entry.weight = w;
                }
            }
            //  Value of element correction factor from this standard for this element
            entry_number++;
            if( records.size() > entry_number && records[entry_number].length() > 0 ) {
                float ecf;
                istringstream value( records[entry_number] );
                value >> ecf;
                if ( ! value || ecf < 0 ) {
                    cout << "Invalid element calibration factor on line " << line_number << ",  Element " << el_entry.element.symbol() << ", " << records[entry_number] << endl;
                    entry_error = true;
                } else {
                    el_entry.ecf = ecf;
                }
            }
            //  Error of element correction factor from this standard for this element (always relative percent)
            entry_number++;
            if( records.size() > entry_number && records[entry_number].length() > 0 ) {
                float ecf_s;
                istringstream value( records[entry_number] );
                value >> ecf_s;
                if ( ! value || ecf_s < 0 ) {
                    cout << "Invalid ecf sigma on line " << line_number << ",  Element " << el_entry.element.symbol() << ", " << records[entry_number] << endl;
                    entry_error = true;
                } else {
                    el_entry.ecf_sigma = ecf_s / 100;
                }
            }
            if ( ! entry_error ) {
                //  Put this new entry in the element list
                //      (function will automatically apply match criteria and transfer non-default info)
                add_element_list_entry( el_entry, standard_entry.element_list );
            } else {
                error = true;
            }
        }   //      if( keyword ==

        //  Process this standard if there has already been a standard being read in
        //      and there is an end of file, a SPECTRUM keyword, or a new STANDARD keyword
        if( ! between_standards && ( !StdCalEntryFile || keyword == "SPECTRUM" || keyword == "STANDARD" ) ) {
            //  Put in a separate entry for each spectrum name using the standard defined so far
            if ( ! entry_error ) {
                if( standard_entry.element_list.size() <= 0 ) {
                    cout << "*** No element list for spectrum file entry on line " << line_number << ". ***" << endl;
                    entry_error = true;
                    error = true;
                    continue;
                }
                //  Make a complete list of elements and given percents
                vector <Element> material_elements;
                vector <float> material_fractions;
                unsigned int iel;
                for( iel=0; iel<standard_entry.element_list.size(); iel++ ) {
                    //  Skip over any element list entries that should not be included in composition
                    if( standard_entry.element_list[iel].percent <= 0 ) continue;
                    if( standard_entry.element_list[iel].type != ELEMENT ) continue;
                    if( standard_entry.element_list[iel].qualifier != NO_QUALIFIER && standard_entry.element_list[iel].qualifier == FORCE
                            && standard_entry.element_list[iel].qualifier == MATRIX ) continue;
                    //  Replace the fraction if this element is already in list
                    bool element_found = false;
                    unsigned int ie;
                    for( ie=0; ie<material_elements.size(); ie++ ) {
                        if( ! ( material_elements[ie] == standard_entry.element_list[iel].element ) ) continue;
                        material_fractions[ie] = standard_entry.element_list[iel].percent / 100;
                        element_found = true;
                        break;
                    }
                    if( ! element_found ) {
                        material_elements.push_back( standard_entry.element_list[iel].element );
                        material_fractions.push_back( standard_entry.element_list[iel].percent / 100 );
                    }
                }
                //  Create an XrayMaterial object for this standard
                XrayMaterial temp_mat( material_elements, material_fractions );
                //  Now see if any oxide or carbonate ratios or uncertainties were entered and put them into the XrayMaterial object
                for( iel=0; iel<standard_entry.element_list.size(); iel++ ) {
                     if( standard_entry.element_list[iel].stoichiometry.formula != PURE_ELEMENT )
                        temp_mat.stoichiometry( standard_entry.element_list[iel].element, standard_entry.element_list[iel].stoichiometry );
                    if( standard_entry.element_list[iel].uncertainty >= 0 )
                        temp_mat.uncertainty( standard_entry.element_list[iel].element, standard_entry.element_list[iel].uncertainty / 100 );
                }
                temp_mat.thickness( thickness_save );
                temp_mat.density( density_save );
                standard_entry.mat = temp_mat;
                //  Add in any comments captured between SPECTRUM keywords
                if( stray_comments.size() > 0 ) {
                    unsigned int i_com;
                    for( i_com=0; i_com<stray_comments.size(); i_com++ ) {
                        standard_entry.preceding_comments.push_back( stray_comments[i_com] );
                    }
                }
                standards_out.push_back( standard_entry );
//                cout << "standard_entry.mat.toString   std: " << (standard_entry.names.size()>0?standard_entry.names[0]:"(no names)") << endl;
//                cout << standard_entry.mat.toString() << endl;
//                cout << endl;
            }
            //  Get rid of any preceding comments since they have now been included in the standards list
            standard_entry.preceding_comments.clear();
            stray_comments.clear();
            //  Set the flag to capture any more comments in case they are before the next standard keyword
            between_standards = true;
        }
        //  STANDARD keyword  (must be handled after processing any old standard above since data is reset)
        if( keyword == "STANDARD" ) {
            //  New standard, reset entry information
            StandardInformation new_entry;  //  Use defaults from StandardInformation struct definition
            //  Indicate that weights were read in ans specified by the user (not true for TXT standards input files)
            new_entry.user_weights = true;
            standard_entry = new_entry;
            //  Put all of the names for this standard into the entry (if any)
            unsigned int ir;
            for( ir=1; ir<records.size(); ir++ ) {
                //  Ignore blank names (sometimes Excel puts in extra commas at the end of a line)
                if( records[ir].length() == 0 ) continue;
                //  parse_records absorbs leading and trailing blanks
                if( records[ir].length() == 1 && records[ir] == " " ) continue;
                standard_entry.names.push_back( records[ir] );
            }
            //  Add in any captured comments
            unsigned int i_com;
            for( i_com=0; i_com<stray_comments.size(); i_com++ ) {
                standard_entry.preceding_comments.push_back( stray_comments[i_com] );
            }
            stray_comments.clear();
            between_standards = false;
        }

        if( entry_error ) error_number++;
        if( entry_error && error_number > MAX_ERROR_MESSAGES ) return -1;    //  Avoid many error messages if wrong format file is opened

	}   //		read and process lines

	StdCalEntryFile.close();

	if( error ) {
		return -2;
	} else {
		return 0;
	}

};
