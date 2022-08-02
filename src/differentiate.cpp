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

// function to differentiate array
//  Oct. 25, 2-13   W. T. ELam   APL/UW

#include "differentiate.h"


using namespace std;

void differentiate( vector <float> &d ) {

      const int n = d.size();
    float t1, t2, t3, t4;

// function to differentiate array
// Ref. Carl-Erik Froberg 'Intro. to Numerical Analysis'
// Addison-Wesley,1965

	if(n < 2) return;
	t1=d[0];
	d[0]=d[1]-t1;
	if(n > 2) {;
		t2=d[1];
		d[1]=0.5 * (d[2] -t1);
		if(n > 4) {;
			t3= d[2];
			d[2]=0.5* (d[3]-t2 -(1./6.)*(d[4]-2.0*t3 + t1) );
			if (n > 6) {;
				int j=n-3;
				int i;
				for( i=3; i<j; i++ ) {;
					t4=d[i];
					d[i]=0.5* (d[i+1]-t3-(1./6.)*(d[i+2]-2.0*t4+t2)
						+(1./30.)*(d[i+3]-3.0*d[i+1]+3.0*t3-t1));
					t1=t2;
					t2=t3;
					t3=t4;
				};
			};
			t4=d[n-3];
			if(n == 5) {;
				t3=t2;
				t2=t1;
			};
			d[n-3]=0.5*(d[n-2]-t3-(1./6.)*(d[n-1]-2.0*t4+t2));
			t2=t4;
		};
		if(n == 3) t2=t1;
		t1=d[n-2];
		d[n-2]=0.5*(d[n-1]-t2);
	};
	d[n-1]=d[n-1] -t1;
	return;
}
