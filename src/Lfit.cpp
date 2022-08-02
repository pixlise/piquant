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

#include "Lfit.h"
#include "XRFconstants.h"
#include <math.h>

//  Adapted from "Numerical Recipes in C"
//  Modified May 27, 2017
//      Added check for not-a-number and infinity to lowerUpperSubst
//      to prevent fit from returning unusable results
//      (should probably switch to using singular value decomposition)

using namespace std;

int lfit( const vector <float> &y, const vector <float> &sig, vector <float> &a,
	vector <float> &var, float &chisq, const vector <float> &funcs, const int np ) {
	int i,j,k;
	float ym,wt,sum,sig2i;
	int ma = a.size();
	vector <float> beta(ma);
	vector <float> afunc(ma);
	vector <float> covar(ma*ma);
	int ndat = y.size();
//		zero covariance and right-hand-side matrices (covar used to store alpha)
	for (j=0;j<ma;j++) {
		for (k=0;k<ma;k++) covar[j*ma+k]=0.0;
		beta[j]=0.0;
	}
//		load matrices with weighted sums
	for (i=0;i<ndat;i++) {
		int ia;
//			move basis function values at x[i] to short vector
		for ( ia=0; ia<ma; ia++ ) afunc[ia] = funcs[ia*np+i];
		ym=y[i];	//	measured value
		sig2i=1.0/(sig[i]*sig[i]);	//	standard deviation of measured value
		for (j=0;j<ma;j++) {
			wt=afunc[j]*sig2i;
			for (k=0;k<=j;k++)
				covar[j*ma+k] += wt*afunc[k];
			beta[j] += ym*wt;
		}
	}
//		fill in other half by symmetry
	for (j=1;j<ma;j++)
		for (k=0;k<j;k++)
			covar[k*ma+j]=covar[j*ma+k];
//		use Lower-Upper decomposition to find solution
	vector <int> index(ma);
	float d;
	int err = lowerUpperDecomp( covar, ma, index, d, ma );
	if ( err != 0 ) return err;
	lowerUpperSubst( covar, ma, index, beta, ma );
//		fill in output vector with fit coefficients
	for (j=0;j<ma;j++) a[j]=beta[j];
//		calculate chi squared
	chisq=0.0;
	for (i=0;i<ndat;i++) {
		int ia;
		for ( ia=0; ia<ma; ia++ ) afunc[ia] = funcs[ia*np+i];
		sum = 0;
		for ( j=0;j<ma;j++) sum += a[j]*afunc[j];
		float diff = (y[i]-sum)/sig[i];
		chisq += diff*diff;
	}
//		calculate variances of fit coefficients
	for (j=0;j<ma;j++) {
		for (k=0;k<ma;k++) beta[k]=0.0;
		beta[j] = 1;
//			use lowerUpperSubst to find one column of the inverse of original covar matrix
		lowerUpperSubst( covar, ma, index, beta, ma );
//			put diagonal element into coefficient variance
		var[j] = beta[j];
	};
	return 0;
};


int lowerUpperDecomp ( vector <float> &a, int n, vector <int> &index, float &d, const int np )
{
	int i,imax,j,k;
	float big,dummy,sum,temp;

	vector <float> pivot(n);
	d=1.0;
	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++) {
			if ((temp=fabs(a[i*np+j])) > big) big=temp;
		};
		if (big == 0.0) return -1;
		pivot[i]=1.0/big;
	};
	for (j=0;j<n;j++) {
		for (i=0;i<j;i++) {
			sum=a[i*np+j];
			for (k=0;k<i;k++) {	//	changed k=1 to k=0  Feb. 22, 2005
				sum -= a[i*np+k]*a[k*np+j];
			};
			a[i*np+j]=sum;
		};
		big=0.0;
		imax = j;
		for (i=j;i<n;i++) {
			sum=a[i*np+j];
			for (k=0;k<j;k++) {
				sum -= a[i*np+k]*a[k*np+j];
			};
			a[i*np+j]=sum;
			if ( (dummy=pivot[i]*fabs(sum)) >= big) {
				big=dummy;
				imax=i;
			};
		};
		if (j != imax) {
			for (k=0;k<n;k++) {
				dummy=a[imax*np+k];
				a[imax*np+k]=a[j*np+k];
				a[j*np+k]=dummy;
			};
			d = -d;
			pivot[imax]=pivot[j];
		};
		index[j]=imax;
		if (a[j*np+j] == 0.0) a[j*np+j]=MINIMUM;
		if (j != n-1) {	//	changed n to n-1  Feb. 22, 2005
			dummy=1.0/(a[j*np+j]);
			for (i=j+1;i<n;i++) {
				a[i*np+j] *= dummy;
			};
		};
	};
	return 0;
};


void lowerUpperSubst ( vector <float> &a, int n, vector <int> &index, vector <float> &rhs, const int np ) {
	int i,ii=-1,ip,j;
	float sum;

	for (i=0;i<n;i++) {
		ip=index[i];
		sum=rhs[ip];
		rhs[ip]=rhs[i];
		if (ii >= 0) {
			for (j=ii;j<=i-1;j++) {
				sum -= a[i*np+j]*rhs[j];
			};
		} else {
			if (sum) {
				ii=i;
			};
		};
		rhs[i]=sum;
		if( isnan( rhs[i] ) || isinf( rhs[i] ) ) rhs[i] = 0;
    };
	for (i=n-1;i>=0;i--) {
		sum=rhs[i];
		for (j=i+1;j<n;j++) {
			sum -= a[i*np+j]*rhs[j];
		};
		rhs[i]=sum/a[i*np+i];
        if( isnan( rhs[i] ) || isinf( rhs[i] ) ) rhs[i] = 0;
	};
};
