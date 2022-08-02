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

// Nick Yang
// 8/9/13
// Boxcar optic
// Modified May 1, 2017 to remove using from header file, to use enum for type, and to make file use more consistent
// Modified May 26, 2017 To add optic efficiency from new breadboard system

#ifndef XrayOptic_h
#define XrayOptic_h

#include <string>

enum XrayOpticType { NO_OPTIC = 0, BOXCAR, PIXL, TRANSMISSION_FILE, NEW_BB, ARRAY_INPUT, PIXL_FM_OPTIC_OLD, PIXL_FM_OPTIC };
//  PIXL_FM_OPTIC_OLD was incorrectly calculated with wrong Be window thickness for flight X-ray tube

struct OpticFileData {
	std::vector<float> vectorDataX;
	std::vector<float> vectorDataY;
	std::vector<float> vectorDataD; //  Derivatives for spline interpolation
};

class XrayOptic
{
private:
	float m_centerEnergy;
	float m_bandwidth;
	float m_maxTransmission;
	bool  m_defaultFlag;
	XrayOpticType	  m_type;
	std::string m_optic_file;
	OpticFileData m_vectorData;
	int GetLineCount( std::string fileName );
	int ParseFile( std::string fileName );

public:
	XrayOptic();
	XrayOptic( const float centerEnergy, const float bandwidth, const float maxTransmission = 1, const XrayOpticType type = NO_OPTIC );
	XrayOptic( const std::string optic_file_in );
	XrayOptic( const std::vector <float> energy_data_in, const std::vector <float> transmission_data_in );
	XrayOptic( const std::vector <float> energy_data_in, const std::vector <float> transmission_data_in, const std::vector <float> derivative_data_in );

	void SetType (XrayOpticType type_in);
	float CheckTransmission( const float energy ) const;

	const   float	GetCenterEnergy()		const	{ return m_centerEnergy; }
	const   float	GetBandwidth()			const	{ return m_bandwidth; }
	const   float	GetMaxTransmission()	const	{ return m_maxTransmission; }
	const   bool	DefaultCheck()			const	{ return m_defaultFlag; }
	const   XrayOpticType		GetType()	const	{ return m_type; }

	std::string toString() const;
};

#endif
