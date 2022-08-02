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

#include <vector>
#include <math.h>
#include "split_component.h"
#include <iostream>         //  **** Needed for debug only ****

// Splits a spectrum component (like background) into several components
//      that can be individually fit via least squares for a better shape  10/20/2020 W. T. Elam

using namespace std;

float split_weight( const float energy, const std::vector <float> &region_list, const int fcn_index ) {
//    cout << "split_weight 1   " << energy << "   " << region_list.size() << "   " << fcn_index << endl;
    //  This function returns weight values for a single component
    //      based on a split into regions delimited by the input array
    //  The argument fcn_index tells which split is being formed (which region)
    if( region_list.size() == 0 ) return 1;
    if( fcn_index < 0 || fcn_index >= region_list.size() ) return -1;
    //  Process left-hand part of region
    if( fcn_index == 0 ) {
        if( energy <= region_list[fcn_index] ) return 1; //  Below start of first region
    } else {
        if( energy < region_list[fcn_index-1] ) return 0;   //  Outside of a central region
        else if( region_list[fcn_index-1] <= energy && energy < region_list[fcn_index] ) {
            float denominator = region_list[fcn_index] - region_list[fcn_index-1];
            if( denominator == 0 ) return 1;    //  Zero-length region, energy must be equal to the boundary
            //  Positive linear slope within the left-hand part of the region
            float weight = ( energy - region_list[fcn_index-1] ) / denominator;
//    cout << "split_weight left   " << region_list[fcn_index-1] << "   " << region_list[fcn_index] << "   " << energy - region_list[fcn_index-1] << "   " << weight << endl;
            return weight;
        }
    }
    //  Process right-hand part of region
    if( fcn_index+1 == region_list.size() ) {
        if( energy > region_list[fcn_index] ) return 1; //  Beyond end of last region
    } else {
        if( region_list[fcn_index+1] < energy ) return 0;   //  Outside of a central region
        else if( region_list[fcn_index] <= energy && energy < region_list[fcn_index+1] ) {
            float denominator = region_list[fcn_index+1] - region_list[fcn_index];
            if( denominator == 0 ) return 1;    //  Zero-length region, energy must be equal to the boundary
            //  Negative linear slope within the right-hand part of the region
            float weight = 1 - ( energy - region_list[fcn_index] ) / denominator;
//    cout << "split_weight right   " << region_list[fcn_index] << "   " << region_list[fcn_index+1] << "   " << energy - region_list[fcn_index] << "   " << weight << endl;
            return weight;
        }
    }

    return 0;
}
