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
#include <iomanip>
#include "read_EMSA_PIXL.h"
#include "write_conditions.h"

//  Written Sept. 29, 2017

	using namespace std;

//      Write configuration info used in FP calculations (from cond_in.conditionsVector with keywords and units)

int write_conditions ( const XRFconditionsInput &cond_in ) {

    cout << "Configuration parameters for fundamental parameters calculations:" << endl;
    const int cond_entries_per_line = 10;
    int ic_debug;
    for( ic_debug = 0; ic_debug < cond_in.conditionsVector.size(); ic_debug++ ) {
        cout << setw(15) << get_EMSA_keyword( ic_debug );
        if( ( ic_debug + 1 ) % cond_entries_per_line == 0 ) cout << endl;
    }
    if( ic_debug % cond_entries_per_line != 0 ) cout << endl;
    for( ic_debug = 0; ic_debug < cond_in.conditionsVector.size(); ic_debug++ ) {
        cout << setw(9) << cond_in.conditionsVector[ic_debug];
        cout << " " << setw(5) << get_EMSA_units( ic_debug, cond_in.conditionsVector[ic_debug] );
        if( ( ic_debug + 1 ) % cond_entries_per_line == 0 ) cout << endl;
    }
    if( ic_debug % cond_entries_per_line != 0 ) cout << endl;
    cout << "      optic file: " << cond_in.optic_file_name << endl;
    cout << "      external tube file: " << cond_in.tube_file_name << endl;
return 0;

};
