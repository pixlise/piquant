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

/****
 *
 * Routine to find energy calibration parameters
 *   for energy-dispersive X-ray spectra
 *
 * Inputs:  a vector with an X-ray spectrum from
 *    an energy-calibration specimen (It must have
 *    only two large, well-separated peaks.  The
 *    energy of the peaks and the location of a
 *    channel between the peaks are given in the
 *    constants below.)
 *
 * Outputs: offset for first channel and energy per channel
 *
 *
 * Header file to declare constants for energy_calibration.cpp
 *
 * Originally written on April 16, 2016 by
 * Aakash Sethi (asethi77@cs.washington.edu)
 * Allen Truong (truona2@uw.edu)
 ****/

// Modified Nov. 8, 2016
//  Change energies of calibration peaks to come from input element list
//  Use exclusion zone around first peak to search for second peak (not fixed split channel as before)
// Modified Nov. 8, 2016
//  Change minimum counts to accept a peak from 200 to 100 (for map spectra with short dwell times)
//  See also modifications in source code file this date

#ifndef __ENERGY_CALIBRATION_H
  const float MIN_COUNT_THRESHOLD = 100;  // Minimum count threshold before peaks are considered too small
  const int PEAK12_SPLIT = 1000;  // A channel number between the two peaks

  #include <vector>
  #include "parse_element_list.h"

  int energy_calibrate( const std::vector <float> &spectrum_xrf_anal,  const std::vector<ElementListEntry> &element_list,
                       float &energy_start, float &energy_per_channel );

#endif // __ENERGY_CALIBRATION_H

//  Returns:
//      -1  Not enough channels in spectrum
//      -2  Not enough counts in peak 1 (lower energy)
//      -3  Not enough counts in peak 2 (higher energy)
//      -4  First element has no K lines
//      -5  Second element has no K lines
