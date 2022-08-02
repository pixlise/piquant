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


/************************************************************************
 * FUNCTION:    box_smth()
 *
 * PURPOSE:     calculates a smoothed spectrum using a boxcar filter
 *
 * CALLED BY:   snipbg, snipbg_lsq
 *
 * PARAMETERS:  y - original spectrum
 *              nchan - number of channels in spectrum
 *              ich1 - first channel number to be smoothed
 *              ich2 - last channel number to be smoothed
 *              iwid - width of the filter (2m+1), iwid<41
 *
 * RETURN(s):   s - smoothed spectrum, only defined between ich1 and ich2
 *
 * REVISION:    1.0
 *
 * NOTES:       written March 5, 2017   W. T. Elam   APL/UW
 ************************************************************************/

#include <vector>
#include <cstdlib>
#include <math.h>
#include "snip.h"

#include <algorithm> // Needed for some compilers (eg MSVC) if min/max is not available

#define _MIN min
#define _MAX max

using namespace std;

void box_smth( const std::vector <float> &y,  std::vector <float> &s,
            const int ich1, const int ich2, const int iwid )

{
 	const int nchan = y.size();
	int iw, m, jch1, jch2;
	float sum, m_norm;
	int i,j;

	iw = _MIN(iwid, nchan);
	m = iw/2;
	m_norm = float( 2 * m + 1 );

	//convolve spectrum with boxcar filter of width ~iw
	jch1 = _MAX(ich1, 0);
	jch2 = _MIN(ich2, nchan - 1);
	for(i = jch1; i <= jch2; i++)
	{
		sum = 0;
		for(j = -m; j <= m; j++)
		{
			int j1 = _MAX(i+j, jch1 );
			j1 = _MIN(j1, jch2 );
			sum += y[j1];
		}
		s[i] = sum / m_norm;
	}
	return;
}//end box_smth()


/************************************************************************
 * FUNCTION:    snipbg()
 *
 * PURPOSE:     calculates a continuum via peak stripping
 *
 * CALLED BY:
 *
 * PARAMETERS:  y - original spectrum
 *              nchan - number of channels in spectrum
 *              ich1 - first channel number of region to calculate the continuum
 *              ich2 - last channel number of region to calculate the continuum
 *              fwhm - width parameter for smoothing and stripping algorithm, set it
 *                     to average FWHM of peaks in the spectrum, typical value 8.0
 *              niter - number of iterations of SNIP algorithm, typical 24
 *
 * RETURN(s):   yback - calculated continuum in the region ich1-ich2
 *
 * REVISION:    2.0
 *
 * NOTES:       wte - 5 March 2017   change Savitsky-Golay smoothing to boxcar
 *                      to avoid occasional instabilities for large or narrow peaks
 *              wte - 18 Dec. 2005
 *					- created header, use vector size for nchan, and eliminate fortranRound
 *				ens - 3 August 2005: revised
 *					- snipbg() ich1 and ich2 values no longer passed by reference
 *				ens - 14 April 2005: this function makes use of the sgsmth() function
 *              ens - 12 April 2005: this algorithm, developed using FORTRAN,
 *              wte - 23 May 2007:  re-worked handling of ends of range
 *              originally appeared in "Handbook of X-Ray Spectrometry"
 *              by R.E.Van Grieken and A.A Markowicz
 ************************************************************************/

void snipbg( const  std::vector <float> &y,  std::vector <float> &yback, const int ich1,
            const int ich2, const int fwhm, const int niter )
{

	const int nchan = y.size();
	const float SQRT2 = 1.4142f;
	const int NREDUC = 8;
	int iw, i1, i2;
	int i,n;

	//smooth spectrum
	iw = fwhm;
	i1 = _MAX(ich1, 0);
	i2 = _MIN(ich2, nchan-1);
	box_smth(y, yback, i1, i2, iw); //data smoothing routine

	//square root transformation over required spectrum region
	for( i = i1; i <= i2; i++)
	{
		yback[i] = sqrt(_MAX(yback[i], (float)0.0));
	}//end for

	//peak stripping
	float redfac = 1;
	for( n = 1; n <= niter; n++)
	{
		//set width, reduce width for last nreduc iterations
		if(n > niter-NREDUC)
			redfac = redfac/SQRT2;
		iw = int( redfac * fwhm + 0.5f );
		for( i = ich1; i <= ich2; i++)
		{
			i1 = _MAX(i-iw, ich1);
			i2 = _MIN(i+iw, ich2);
			yback[i] = _MIN(yback[i], (float)(0.5*(yback[i1]+yback[i2])));
		}
	}

	//back transformation
	for(i = ich1; i <= ich2; i++)
		yback[i] = yback[i]*yback[i];

	return;
}//end snipbg


/************************************************************************
 * FUNCTION:    snipbg_lsq()
 *
 * PURPOSE:     calculates a continuum via peak stripping with improved
 *              estimation via least squares
 *
 * CALLED BY:
 *
 * PARAMETERS:  y - original spectrum
 *              nchan - number of channels in spectrum
 *              ich1 - first channel number of region to calculate the continuum
 *              ich2 - last channel number of region to calculate the continuum
 *              fwhm - width parameter for smoothing and stripping algorithm, set it
 *                     to average FWHM of peaks in the spectrum, typical value 8.0
 *              niter - number of iterations of SNIP algorithm, typical 24
 *
 * RETURN(s):   yback - calculated continuum in the region ich1-ich2
 *
 * REVISION:    3.0
 *
 * NOTES:       wte - 5 March 2017   change Savitsky-Golay smoothing to boxcar
 *                      to avoid occasional instabilities for large or narrow peaks
 *              wte - 6 June 2014   add least squares adjustment to SNIP result
 *              ens - 3 August 2005: revised
 *					- sgsmth() ich1 and ich2 values no longer passed by reference
 *				ens - 12 April 2005: this algorithm, developed using FORTRAN,
 *              wte - 23 May 2007:  re-worked handling of ends of range
 *              originally appeared in "Handbook of X-Ray Spectrometry"
 *              by R.E.Van Grieken and A.A Markowicz
 ************************************************************************/

void snipbg_lsq( const  std::vector <float> &y,  std::vector <float> &yback, const int ich1,
            const int ich2, const int fwhm, const int niter )
{

	const int nchan = y.size();
	const float SQRT2 = 1.4142f;
	const int NREDUC = 8;
	int iw, i1, i2;
	int i,n;

	//smooth spectrum
	iw = fwhm;
	i1 = _MAX(ich1, 0);
	i2 = _MIN(ich2, nchan-1);
	box_smth(y, yback, i1, i2, iw); //data smoothing routine

	//square root transformation over required spectrum region
	for( i = i1; i <= i2; i++)
	{
		yback[i] = sqrt(_MAX(yback[i], (float)0.0));
	}//end for

	//peak stripping
	float redfac = 1;
	for( n = 1; n <= niter; n++)
	{
		//set width, reduce width for last nreduc iterations
		if(n > niter-NREDUC)
			redfac = redfac/SQRT2;
		iw = int( redfac * fwhm + 0.5f );
		for( i = ich1; i <= ich2; i++)
		{
			i1 = _MAX(i-iw, ich1);
			i2 = _MIN(i+iw, ich2);
			yback[i] = _MIN(yback[i], (float)(0.5*(yback[i1]+yback[i2])));
		}
	}

	//back transformation
	for(i = ich1; i <= ich2; i++)
		yback[i] = yback[i]*yback[i];

    // the SNIP algorithm gives a background that is slightly below
    //  the average value for the spectrum, causing false positives
    //  in peak fits, so attempt to use least squares to adjust it
    // form least squares sums, ignoring peak regions (only include
    //  regions that are near the background found by SNIP above
    float y_sum = 0;
    float f_sum = 0;
	for(i = ich1; i <= ich2; i++) {
        // skip this channel if counts are more than 3 sigma different than background
        if( fabs( y[i] - yback[i] ) > 3 * sqrt( yback[i] ) ) continue;
        y_sum += y[i] * yback[i];
        f_sum += yback[i] * yback[i];
    }
    // now adjust yback for a least-squares fit to the measured spectrum
	for(i = ich1; i <= ich2; i++) yback[i] *= y_sum / f_sum;

	return;
}//end snipbg_lsq

/************************************************************************
 * FUNCTION:    snipbg()
 *
 * PURPOSE:     calculates a continuum via peak stripping
 *
 * CALLED BY:
 *
 * PARAMETERS:  y - original spectrum
 *              nchan - number of channels in spectrum
 *              ich1 - first channel number of region to calculate the continuum
 *              ich2 - last channel number of region to calculate the continuum
 *              fwhm - width parameter for smoothing and stripping algorithm, set it
 *                     to average FWHM of peaks in the spectrum, typical value 8.0
 *              niter - number of iterations of SNIP algorithm, typical 24
 *              ich1_2 - first channel number of region where fwhm2 is used
 *              ich2_2 - last channel number of region where fwhm2 is used
 *              fwhm2 - width parameter for smoothing and stripping algorithm n region 2
 *
 * RETURN(s):   yback - calculated continuum in the region ich1-ich2
 *
 * REVISION:    2.0
 *
 * NOTES:       wte - 27 July 2018   change to two-region stripping for better
 *                      background removal under peaks (inlcuding Compton peak)
 *                      and for better fit to broad hump in continuum near center
 *              wte - 5 March 2017   change Savitsky-Golay smoothing to boxcar
 *                      to avoid occasional instabilities for large or narrow peaks
 *              wte - 18 Dec. 2005
 *					- created header, use vector size for nchan, and eliminate fortranRound
 *				ens - 3 August 2005: revised
 *					- snipbg() ich1 and ich2 values no longer passed by reference
 *				ens - 14 April 2005: this function makes use of the sgsmth() function
 *              ens - 12 April 2005: this algorithm, developed using FORTRAN,
 *              wte - 23 May 2007:  re-worked handling of ends of range
 *              originally appeared in "Handbook of X-Ray Spectrometry"
 *              by R.E.Van Grieken and A.A Markowicz
 ************************************************************************/

void snipbg_2zone( const  std::vector <float> &y,  std::vector <float> &yback, const int ich1,
            const int ich2, const int fwhm, const int niter,
            const int ich1_2z, const int ich2_2z, const int fwhm2 )
{

	const int nchan = y.size();
	const float SQRT2 = 1.4142f;
	const int NREDUC = 8;
	int iw, i1, i2;
	int i,n;

	//smooth spectrum
	iw = fwhm;
	i1 = _MAX(ich1, 0);
	i2 = _MIN(ich2, nchan-1);
	box_smth(y, yback, i1, i2, iw); //data smoothing routine

	//square root transformation over required spectrum region
	for( i = i1; i <= i2; i++)
	{
		yback[i] = sqrt(_MAX(yback[i], (float)0.0));
	}//end for

	//peak stripping
	float redfac = 1;
	for( n = 1; n <= niter; n++)
	{
		//set width, reduce width for last nreduc iterations
		if(n > niter-NREDUC)
			redfac = redfac/SQRT2;
		iw = int( redfac * fwhm + 0.5f );
		int iw2 = int( redfac * fwhm2 + 0.5f ); //  2 zone
		for( i = ich1; i <= ich2; i++)
		{
			if( ich1_2z > 0 && ich2_2z > 0 && iw2 > 0 && i >= ich1_2z && i <= ich2_2z ) {
                i1 = i-iw2; //  2 zone
                i2 = i+iw2; //  2 zone
            } else {
                i1 = i-iw;
                i2 = i+iw;
            }
			i1 = _MAX(i1, ich1);
			i2 = _MIN(i2, ich2);
			yback[i] = _MIN(yback[i], (float)(0.5*(yback[i1]+yback[i2])));
		}
	}

	//back transformation
	for(i = ich1; i <= ich2; i++)
		yback[i] = yback[i]*yback[i];

	return;
}//end snipbg_2zone

