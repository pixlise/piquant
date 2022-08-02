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

#ifndef Lfit_h
#define Lfit_h

#include <vector>

//	Linear least squares fitting to functions stored in array funcs
//		from Press, Flannery, Teukosky, and Vetterling, Numerical Recipes
//			(Cambridge Univ. Press, Cambridge) 1986.  ISBN 0 521 30811 9.
//		(roughly translated from Fortran, pages 509 - 515)
//  Modified May 26, 2017 to remove using namespace std; from include file

int lfit( const std::vector <float> &y, const std::vector <float> &sig, std::vector <float> &a,
	std::vector <float> &var, float &chisq, const std::vector <float> &funcs, const int np );
//	y		data points to be fit (input)
//	sig		individual standard deviations of data points (input)
//	a		coefficients of fit (output)
//	var		variance matrix (output)
//	chisq	chi squared fit parameter (output)
//	funcs	functions to use in fit (input)
//			index of funcs is i * np + j, where i is index into a (and var) and j is index into y
//	np		physical size of 1st dimension of funcs vector (usually also the physical dimension of y and sig)

//	matrix solvers
int lowerUpperDecomp( std::vector <float> &a, int n, std::vector <int> &index, float &d, const int np );
void lowerUpperSubst( std::vector <float> &a, int n, std::vector <int> &index, std::vector <float> &rhs, const int np );





#endif
