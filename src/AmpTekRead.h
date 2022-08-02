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

#ifndef AmpTekRead_h
#define AmpTekRead_h

#include <vector>

using namespace std;

//	reads asc files from AmpTek MCA software    Jan. 2006   W. T. Elam  APL/UW

// Modified Feb. 20, 2017
//  Put in defaults for most variables in AmpTekSpec
//  To return zero eV start and unity energy per channel if no calibration (originally returned nan)

struct AmpTekSpec {
	string fileID;
	string dataType;
	string description;
	float gain = 1;
	int threshold = 0;
	int live_mode = 0;
	float preset_time = 0;
	float live_time = 1;
	float real_time = 1;
	string start_time;
	int serial_number;
	string cal_label;
	vector <float> cal_channel;
	vector <float> cal_energy;
	float ev_ch = 0;
	float ev_start = 0;
	vector <float> spectrum;
};

int amptek_read ( string fileName, AmpTekSpec &sp );

//  integer returns
//		0	file read OK
//		-1	can't open file
//		-2	can't interpret file
//      -3  unexpected file read error or end of file

#endif
