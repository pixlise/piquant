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

#ifndef Sewell_tube_calc_h
#define Sewell_tube_calc_h

float Sewell_j( float z );
float Sewell_eta( float z, float energy );
float Sewell_r( float u0, float eta, float tilt = 0 );
float Sewell_S_lines( float u0, float j, float ee, float z_a );
float Sewell_S_continuum( float u0, float j, float ee, float z_a );
float Sewell_stopping( float energy, float j, float z_a );
float Sewell_Q( float energy, float ec );
float Sewell_pz( float j, float energy, float eta, float u0,
			   float z_a, float tilt = 0 );
float Sewell_h( float u0, float zBar, float eta, float tilt = 0 );
float Sewell_pz_m( float pz, float u0, float z, float tilt = 0 );
float Sewell_pz_r( float pz, float pz_m, float h );
float Sewell_phi_pz( float pz, float pz_m, float pz_r, float h );
float Sewell_f( float chi, float pz_m, float pz_r, float pz, float h );
float Reed_f( float cj, float tau_ij, float sigma_j, float ri, float omega_j, float ai,
				float aj, float ui, float uj, float chi_i, float sigmaLenard );

#endif
