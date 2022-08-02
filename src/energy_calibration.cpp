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

/**
 * Algorithm to determine XRF energy per channel and energy start
 * based on XRF inputs. Developed at the request of Tim Elam
 * (wtelam@apl.uw.edu)
 *
 * April 16, 2016
 *
 * Aakash Sethi (asethi77@cs.washington.edu)
 * Allen Truong (truona2@uw.edu)
 */
 // Adapted for PIXL Ca/Zr calibration bead and changed asserts to error returns
 //     W. T. Elam   April 28, 2016
// Modified Nov. 8, 2016
//  Change energies of calibration peaks to come from input element list
//  Use exclusion zone around first peak to search for second peak (not fixed split channel as before)
// Modified Nov. 9, 2016 to fix bug in exclusion zones when using element list
// Modified Nov. 21, 2016 to remove writes to cout (accidentally left in during above change)
// Modified Dec. 8, 2016 to to increase size of exclusion around first peak (in case 2nd peak is small)
//  Modified Jan. 2, 2017 to calculate slope and use zero offset if only one element given in list
//  Modified Jan. 18, 2017 to put offset then slope into argument list (to be consistent throughout code)
//  Modified June 9, 2017 to check for not enough channels when element list has 1 or more elements
//  Modified Sept.13, 2017 fix normalization for alpha 1 & 2 peaks, fix indices for L lines, fix == in level 1 setup


//#include <iostream>
#include "energy_calibration.h"
#include "XrayEdge.h"
#include "XrayLines.h"

/**
 * Function to compute the first order least squares fit
 */
static void least_squares_fit(int num_elems, const float *const x_arr,
      const float *const y_arr, float &intercept, float &slope)
{
  //assert(num_elems != 0);
  //assert(sizeof(x_arr) == sizeof(y_arr));

  float sumx = 0;
  float sumy = 0;
  float sumx2 = 0;
  float sumxy = 0;

  int i;
  for (i = 0; i < num_elems; i++)
  {
    sumx += x_arr[i];
    sumx2 += x_arr[i] * x_arr[i];
    sumy += y_arr[i];
    sumxy += x_arr[i] * y_arr[i];
  }

  slope = (num_elems * sumxy - sumx * sumy) / (num_elems * sumx2 - sumx * sumx);
  intercept = (sumx2 * sumy - sumx * sumxy) / (num_elems * sumx2 - sumx * sumx);
}

/**
 * Finds the channel in data_array with the maximum number of counts and
 * returns the index of that channel and the counts.
 */
static void find_max_channel( const std::vector <float> &data_array, const int start_channel,
    const int end_channel, int &max_channel, int &max_counts) {
  //assert(start_channel <= end_channel);

    max_channel = start_channel;
    max_counts = data_array[start_channel];
    int i;
    for (i = start_channel; i <= end_channel; i++) {
        if (data_array[i] > max_counts) {
        max_counts = data_array[i];
        max_channel = i;
    }
    }
}

/**
 * Takes a data_array and a central channel index and returns the average
 * channel count value and norm of it and the surrounding four channels.
 */
static void five_point_average(const int max_chan_index, const std::vector <float> &data_array,
    float &channel_average, float &channel_norm) {
    channel_average = 0;
    channel_norm = 0;
    for (int k = max_chan_index - 2; k <= max_chan_index + 2; k++) {
        channel_average += data_array[k] * (k - 1);
        channel_norm += data_array[k];
    }
    if( channel_norm != 0 ) channel_average /= channel_norm;
}

/**
 * This routine assumes that the highest points in the spectrum are well-separated
 peaks from the first two elements in the list.  If there is only one element,
 then calculate an energy per channel with zero offset.
 If there are no elements, use the default energies given below.
 * For example, a fused bead of Ca and Zr and gives only two very large peaks.
 */
int energy_calibrate( const std::vector <float> &spectrum_xrf_anal,  const std::vector<ElementListEntry> &element_list,
                       float &energy_start, float &energy_per_channel ) {

    float num_channels = spectrum_xrf_anal.size();
    float pk1_chan, pk2_chan, chan_norm;
    int max_chan, max_counts, start_chan, end_chan;
    float Peak1_energy = 3691;  // Ca Ka, eV  (PIXL calibration bead)
    float Peak2_energy = 15776;  // Zr Ka, eV  (PIXL calibration bead)
    int lowest_channel = num_channels / 100;    // Avoid any big noise peaks
    int highest_channel = num_channels - 10;    // Avoid any pileup or registers stored in the end channels

    if( element_list.size() < 1 ) {

        if( num_channels < 3 * PEAK12_SPLIT / 2 ) return -1;

        // Find the search ranges for Zn from the number of channels
        start_chan = lowest_channel;  // Avoid any big noise peaks
        end_chan = PEAK12_SPLIT;

    } else {
        if( highest_channel < 0 ) return -1;
        //    Find the peak energies from the first two elements in the list
        EdgeIndex level1 = K1;
        if( element_list[0].quant_level == K_LEVEL ) level1 = K1;
        if( element_list[0].quant_level == L_LEVEL ) level1 = L3;
        if( element_list[0].quant_level == M_LEVEL ) level1 = M5;
        if( element_list[0].quant_level == N_LEVEL ) level1 = N5;
        XrayEdge edge1( element_list[0].element, level1 );
        XrayLines lines1( edge1 );
        //  Use K alpha 1 and K alpha 2 lines, line indices 0 and 1
        int line_index_1 = 0;
        int line_index_2 = 1;
        if( level1 == L3 ) {
            //  Line indices for L alpha 1 and L alpha 2
            line_index_1 = 1;
            line_index_2 = 2;
        }
        if( level1 == M5 || level1 == N5 ) {
            //  Only one M alpha line and only one N alpha line
            line_index_2 = -1;
        }
        if( lines1.numberOfLines() <= line_index_1 ) return -4;
        Peak1_energy = lines1.relative( line_index_1 ) * lines1.energy( line_index_1 );     //  alpha 1 line
        float norm = lines1.relative( line_index_1 );
        if( line_index_2 >= 0 && lines1.numberOfLines() > line_index_2 ) {
            Peak1_energy += lines1.relative( line_index_2 ) * lines1.energy( line_index_2 );   //  alpha 2 line
            norm += lines1.relative( line_index_2 );
        }
        Peak1_energy /= norm;
        Peak2_energy = -1;
//std::cout << "E1 levels " << element_list[0].element.symbol() << " " << element_list[0].quant_level << "  " << level1 << "  " << edge1.level() << "  " << lines1.energy( line_index_1 ) << std::endl;
//std::cout << "E1 energies " << lines1.relative( line_index_1 ) << "  " << lines1.energy( line_index_1 ) << "  " << lines1.relative( line_index_2 ) << "  " << lines1.energy( line_index_2 ) << "  " << norm << "  " << Peak1_energy << std::endl;
        if( element_list.size() > 1 ) {
            EdgeIndex level2 = K1;
            if( element_list[1].quant_level == K_LEVEL ) level2 = K1;
            if( element_list[1].quant_level == L_LEVEL ) level2 = L3;
            if( element_list[1].quant_level == M_LEVEL ) level2 = M5;
            if( element_list[1].quant_level == N_LEVEL ) level2 = N5;
            XrayEdge edge2( element_list[1].element, level2 );
            XrayLines lines2( edge2 );
            //  Use K alpha 1 and K alpha 2 lines, line indices 0 and 1
            int line_index_1 = 0;
            int line_index_2 = 1;
            if( level2 == L3 ) {
                //  Line indices for L alpha 1 and L alpha 2
                line_index_1 = 1;
                line_index_2 = 2;
            }
            if( level2 == M5 || level2 == N5 ) {
                //  Only one M alpha line and only one N alpha line
                line_index_2 = -1;
            }
            if( lines2.numberOfLines() <= line_index_1 ) return -5;
            Peak2_energy = lines2.relative( line_index_1 ) * lines2.energy( line_index_1 );     //  alpha 1 line
            norm = lines2.relative( line_index_1 );
            if( line_index_2 >= 0 && lines2.numberOfLines() > line_index_2 ) {
                Peak2_energy += lines2.relative( line_index_2 ) * lines2.energy( line_index_2 );   //  alpha 2 line
                norm += lines2.relative( line_index_2 );
            }
            Peak2_energy /= norm;
//std::cout << "E2 levels " << element_list[1].element.symbol() << " " << element_list[1].quant_level << "  " << level2 << "  " << edge2.level() << "  " << lines2.energy( line_index_1 ) << std::endl;
//std::cout << "E2 energies " << lines2.relative( line_index_1 ) << "  " << lines2.energy( line_index_1 ) << "  " << lines2.relative( line_index_2 ) << "  " << lines2.energy( line_index_2 ) << "  " << norm << "  " << Peak2_energy << std::endl;
        }
        //   Switch places if not in ascending order by energy
        if( Peak1_energy > Peak2_energy ) {
            float temp = Peak2_energy;
            Peak2_energy = Peak1_energy;
            Peak1_energy = temp;
        }

        // Use the full spectrum search range
        start_chan = lowest_channel;
        end_chan = highest_channel;

    }

    // Search the lower spectrum for the first peak
    find_max_channel(spectrum_xrf_anal, start_chan, end_chan, max_chan, max_counts);
//    std::cout << "1: " << start_chan << "  " << end_chan << "  " << max_chan << "  " << max_counts << std::endl;

    if( max_counts <= MIN_COUNT_THRESHOLD ) return -2;

    // Take a 5-point average to find the center channel of the peak
    pk1_chan = 0;
    chan_norm = 0;
    five_point_average(max_chan, spectrum_xrf_anal, pk1_chan, chan_norm);

    // Repeat for second peak
    if( element_list.size() < 1 ) {
        // Search the upper spectrum for the second peak
        start_chan = PEAK12_SPLIT;
        end_chan = highest_channel;  // Avoid any pileup in last channels
        // Search the spectrum for the second peak
        find_max_channel(spectrum_xrf_anal, start_chan, end_chan, max_chan, max_counts);
    } else if( Peak1_energy > 0 ) {
        //  Exclude a region around the peak that was found and search the rest of the spectrum
        //  Upper part of spectrum first
        start_chan = 1.15 * pk1_chan;    //  avoid other lines associated with first peak
        end_chan = highest_channel;
        int max_chan_upper, max_counts_upper;
        find_max_channel(spectrum_xrf_anal, start_chan, end_chan, max_chan_upper, max_counts_upper);
//std::cout << "L: " << start_chan << "  " << end_chan << "  " << max_chan_upper << "  " << max_counts_upper << std::endl;
       //  Now lower part of spectrum
        start_chan = lowest_channel;
        end_chan = 0.9 * pk1_chan;
        find_max_channel(spectrum_xrf_anal, start_chan, end_chan, max_chan, max_counts);
//std::cout << "U: " << start_chan << "  " << end_chan << "  " << max_chan << "  " << max_counts << std::endl;
        if( max_counts_upper > max_counts ) {
            max_counts = max_counts_upper;
            max_chan = max_chan_upper;
        }

    }

    //  If only one element given in list, calculate slope and use zero offset
    if( element_list.size() == 1 ) {
        energy_start = 0;
        energy_per_channel = Peak2_energy / pk1_chan;   //  Peak1_energy and Peak2_energy switched above
        return 1;   //  Flag to indicate that energy is from one peak only
    }

    //assert(max_counts >= MIN_COUNT_THRESHOLD);
    if( max_counts <= MIN_COUNT_THRESHOLD ) return -3;

    // Take a 5-point average to find the center channel of the peak
    pk2_chan = 0;
    chan_norm = 0;
    five_point_average(max_chan, spectrum_xrf_anal, pk2_chan, chan_norm);

    //   Switch places if not in ascending order by channel number
    if( pk1_chan > pk2_chan ) {
        float temp = pk2_chan;
        pk2_chan = pk1_chan;
        pk1_chan = temp;
    }
//std::cout << Peak1_energy << "  " << pk1_chan << "  " << Peak2_energy << "  " << pk2_chan << std::endl;
    float channel_avgs[2];
    float channel_energy_levels[2];

    channel_avgs[0] = pk1_chan;
    channel_avgs[1] = pk2_chan;

    channel_energy_levels[0] = Peak1_energy;
    channel_energy_levels[1] = Peak2_energy;

    least_squares_fit(sizeof(channel_avgs) / sizeof(float), (float*) channel_avgs, (float*) channel_energy_levels, energy_start, energy_per_channel);

//std::cout << Peak1_energy << "  " << pk1_chan << "  " << Peak2_energy << "  " << pk2_chan << std::endl;

    return 0;
}

