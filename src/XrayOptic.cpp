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
// 9-17-2013
//
// Adding file parsing support, prints data with a significant loss of precision
// Need to test struct data retention
// Can select the emumerator type via .xsp file
// Changed to zeros for center and bandwidth for file WTE Sep. 30, 2013
// Modified Nov. 21, 2016 to include generic PIXL optic transmission data from include file (no external file to find)
// Modified May 1, 2017 to use optic file in a consistent way and to use strongly-typed enum for optic type
// Modified May 26, 2017 to include efficiency function for new breadboard optic (from Chris Heirwegh)
//                          Also make sure transmission is not negative from interpolation
//  Modified June 26, 2017
//      Fix handling of transmission for optic function from file
//      Add throw if file not found to allow error return to GUI
//      Fix handling of last line in file during read
//  Modified Oct. 22, 2020
//      Add constructor and transmission calculation for data input in arrays, with spline interpolation
//  Modified Oct. 30, 2020
//      Changed to extrapolate below energy arrays to lowest transmission value (flat extrapolation)
//  Modified Nov. 2, 2020
//      Add PIXL FM optic as input from include file, with spline fit
//  Modified Dec. 28, 2020
//      Improvements to optic response function, new optic response clculations, results into FM_OpticResponse_Nov2020.h
//  Modified Apr. 28, 2021  Fix bug in checking size of derivative data in 3 vector constructor
//                          Change PIXL FM optic type (#7) to PIXL_FM_OPTIC_OLD, add new PIXL_FM_OPTIC (#8) calculated with correct Be window for X-ray tube


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include "stdlib.h"
#include "XrayOptic.h"
#include "spline.h"
#include "interp.h"

#include "Energy_vs_Efficiency_smoothed.h"
#include "Breadboard_optic_efficiency_May2017.h"
#include "FM_OpticResponse_Nov2020.h"

using namespace std;

void PrintVector( vector<float> vectorData );

// Default constructor, 100% transmission, zeros everywhere else, default flag set to true
XrayOptic::XrayOptic ( )
{
	m_defaultFlag = true;
	m_centerEnergy = 0;
	m_bandwidth = 0;
	m_maxTransmission = 1;
	m_type = NO_OPTIC;
}


XrayOptic::XrayOptic( const float centerEnergy, const float bandwidth, const float maxTransmission, const XrayOpticType type )
{


	//  error if this constructor is called with type TRANSMISSION_FILE;
    if (m_type == TRANSMISSION_FILE) return;


	m_defaultFlag = false;				// Toggling flag to false so we know it has been defined with real values
	m_centerEnergy = centerEnergy;
	m_bandwidth = bandwidth;
	m_maxTransmission = maxTransmission;
	m_type = type;

	// Load in the optic characteristics if old PIXL breadboard optic was selected
	if (m_type == PIXL)
	{
        m_centerEnergy = 0;
        m_bandwidth = 0;
		// Properly size the vectors to fit the data
		m_vectorData.vectorDataX.resize(N_Energy_vs_Efficiency_smoothed);
		m_vectorData.vectorDataY.resize(N_Energy_vs_Efficiency_smoothed);
        //  Move the data into vectors
		int io;
		for( io=0; io<N_Energy_vs_Efficiency_smoothed; io++ ) {
		    m_vectorData.vectorDataX.at(io) = X_Energy_vs_Efficiency_smoothed[io];
		    m_vectorData.vectorDataY.at(io) = Y_Energy_vs_Efficiency_smoothed[io];
		};
    }

	// Load in the optic characteristics if new PIXL breadboard optic was selected
	if (m_type == NEW_BB)
	{
        m_centerEnergy = 0;
        m_bandwidth = 0;
		// Properly size the vectors to fit the data
		m_vectorData.vectorDataX.resize(N_Energy_vs_Efficiency_smoothed_newBB);
		m_vectorData.vectorDataY.resize(N_Energy_vs_Efficiency_smoothed_newBB);
        //  Move the data into vectors
		int io;
		for( io=0; io<N_Energy_vs_Efficiency_smoothed_newBB; io++ ) {
		    m_vectorData.vectorDataX.at(io) = X_Energy_vs_Efficiency_smoothed_newBB[io];
		    m_vectorData.vectorDataY.at(io) = Y_Energy_vs_Efficiency_smoothed_newBB[io];
		};
	}
	//  Derivative data all zeros for linear interpolation
    vector <float> temp_vec( m_vectorData.vectorDataX.size(), 0 );
    m_vectorData.vectorDataD = temp_vec;

	// Load in the optic characteristics if PIXL FM optic was selected
    //  PIXL_FM_OPTIC_OLD was incorrectly calculated with wrong Be window thickness for flight X-ray tube
	if (m_type == PIXL_FM_OPTIC_OLD)
	{
        m_centerEnergy = 0;
        m_bandwidth = 0;
		// Properly size the vectors to fit the data
		m_vectorData.vectorDataX.resize(N_FM_OpticResp_7);
		m_vectorData.vectorDataY.resize(N_FM_OpticResp_7);
		m_vectorData.vectorDataD.resize( N_FM_OpticResp_7, 0 );
        //  Move the data into vectors
		int io;
		for( io=0; io<N_FM_OpticResp; io++ ) {
		    m_vectorData.vectorDataX.at(io) = X_FM_OpticResp_7[io];
		    m_vectorData.vectorDataY.at(io) = Y_FM_OpticResp_7[io];
		//  Use spline fit from include file (with discontinuity)
		    m_vectorData.vectorDataD.at(io) = D_FM_OpticResp_7[io];
		};
/*        //  Perform a spline fit to the calculated response - only after discontinuity at 3 keV
        float initial_slope = 0;
        if( m_vectorData.vectorDataX.size() > 3 && m_vectorData.vectorDataX[3] - m_vectorData.vectorDataX[2] != 0 )
            initial_slope = ( m_vectorData.vectorDataY[3] - m_vectorData.vectorDataY[2] ) / ( m_vectorData.vectorDataX[3] - m_vectorData.vectorDataX[2] ) / 4;
//        spline( m_vectorData.vectorDataX, m_vectorData.vectorDataY, initial_slope, 0, m_vectorData.vectorDataD );
        const int start_spline = 2;
        const int n_spline = m_vectorData.vectorDataX.size() - start_spline;
        vector <float> tempX( n_spline );
        vector <float> tempY( n_spline );
        int j;
        for( j=0; j<n_spline; j++ ) {
            tempX[j] = m_vectorData.vectorDataX[j+start_spline];
            tempY[j] = m_vectorData.vectorDataY[j+start_spline];
        }
        spline( tempX, tempY, initial_slope, 0, temp_vec );
        //  Put derivatives for spline interpolation into class data (zeros will cause linear interpolation)
        for( j=0; j<n_spline; j++ ) m_vectorData.vectorDataD[j+start_spline] = temp_vec[j];
        cout << "   Optic derivatives   ";
        cout.precision( 8 );
        for( j=0; j<m_vectorData.vectorDataD.size(); j++ ) cout << ",  " << m_vectorData.vectorDataD[j];
        cout << endl;
*/	}

	// Load in the optic characteristics if PIXL FM optic was selected
	if (m_type == PIXL_FM_OPTIC)
	{
        m_centerEnergy = 0;
        m_bandwidth = 0;
		// Properly size the vectors to fit the data
		m_vectorData.vectorDataX.resize(N_FM_OpticResp);
		m_vectorData.vectorDataY.resize(N_FM_OpticResp);
		m_vectorData.vectorDataD.resize( N_FM_OpticResp, 0 );
        //  Move the data into vectors
		int io;
		for( io=0; io<N_FM_OpticResp; io++ ) {
		    m_vectorData.vectorDataX.at(io) = X_FM_OpticResp[io];
		    m_vectorData.vectorDataY.at(io) = Y_FM_OpticResp[io];
		//  Use spline fit from include file (with discontinuity)
		    m_vectorData.vectorDataD.at(io) = D_FM_OpticResp[io];
		};
	}

}


XrayOptic::XrayOptic( const std::string optic_file_in )
{
    m_maxTransmission = 1;

//  Read optic transmission function from file if type = TRANSMISSION_FILE
	m_defaultFlag = false;				// Toggling flag to false so we know it has been defined with real values
	m_type = TRANSMISSION_FILE;
    m_centerEnergy = 0;
    m_bandwidth = 0;
    m_optic_file = optic_file_in;
    // Properly size the vectors to fit the text file
    int lineCount = GetLineCount( m_optic_file );
    if( lineCount <= 0 ) return;
    m_vectorData.vectorDataX.resize( lineCount );
    m_vectorData.vectorDataY.resize( lineCount );

    ParseFile( m_optic_file );

    //cout << "Optical parameters loaded via: " << fileName << endl;
    //cout << setprecision(5);
    //PrintVector(m_vectorData.vectorDataX);
    //PrintVector(m_vectorData.vectorDataY);
}

XrayOptic::XrayOptic( const std::vector <float> energy_data_in, const std::vector <float> transmission_data_in ) {
	m_defaultFlag = false;				// Toggling flag to false so we know it has been defined with real values
	m_centerEnergy = 0;
	m_bandwidth = 0;
	m_type = ARRAY_INPUT;
	// Load in the optic characteristics from the input arrays
	int max_size = energy_data_in.size();
	if( transmission_data_in.size() < max_size ) max_size = transmission_data_in.size();
	float max_trans = 0;
    // Properly size the vectors to fit the data
    m_vectorData.vectorDataX.resize(max_size);
    m_vectorData.vectorDataY.resize(max_size);
    //  Move the data into vectors
    int io;
    for( io=0; io<max_size; io++ ) {
        m_vectorData.vectorDataX.at(io) = energy_data_in[io];
        m_vectorData.vectorDataY.at(io) = transmission_data_in[io];
        if( max_trans < transmission_data_in[io] ) max_trans = transmission_data_in[io];
    };
	m_maxTransmission = max_trans;
	//  Derivative data all zeros for linear interpolation
    vector <float> temp_vec( m_vectorData.vectorDataX.size(), 0 );
    m_vectorData.vectorDataD = temp_vec;
};

XrayOptic::XrayOptic( const std::vector <float> energy_data_in, const std::vector <float> transmission_data_in, const std::vector <float> derivative_data_in ) {
	m_defaultFlag = false;				// Toggling flag to false so we know it has been defined with real values
	m_centerEnergy = 0;
	m_bandwidth = 0;
	m_type = ARRAY_INPUT;
	// Load in the optic characteristics from the input arrays
	int max_size = energy_data_in.size();
	if( transmission_data_in.size() < max_size ) max_size = transmission_data_in.size();
	float max_trans = 0;
    // Properly size the vectors to fit the data
    m_vectorData.vectorDataX.resize(max_size);
    m_vectorData.vectorDataY.resize(max_size);
    m_vectorData.vectorDataD.resize(max_size);
    //  Move the data into vectors
    int io;
    for( io=0; io<max_size; io++ ) {
        m_vectorData.vectorDataX.at(io) = energy_data_in[io];
        m_vectorData.vectorDataY.at(io) = transmission_data_in[io];
        if( max_trans < transmission_data_in[io] ) max_trans = transmission_data_in[io];
        if( io < derivative_data_in.size() ) m_vectorData.vectorDataD.at(io) = derivative_data_in[io];
        else m_vectorData.vectorDataD.at(io) = 0;
    };
	m_maxTransmission = max_trans;
};



void XrayOptic::SetType( XrayOpticType type_in )
{
	m_type = type_in;
}

// Returns the transmission depending on the type and energy
float XrayOptic::CheckTransmission( const float energy ) const
{
//    cout << "CheckTransmission " << energy << endl;

	switch( GetType() )
	{
		// Case 1: Vacuum, 100% transmission
		case NO_OPTIC:
			return 1;

		// Case 2: Boxcar, 0% transmission outside of boxcar boundries
		case BOXCAR:
			if (energy >= (m_centerEnergy - (m_bandwidth / 2.0)) && energy <= (m_centerEnergy + (m_bandwidth / 2.0)))
			{
//				cout << "Boxcar Transmission: " << m_maxTransmission << endl;
				return m_maxTransmission;
			}
			else
			{
//				cout << "Boxcar Transmission: 0" << endl;
				return 0;
			}
			break;
		// Case 3: Rover (based on old breadboard optic)
		case PIXL:
		//  Case 4: Transmission from file read in with ParseFile()
		case TRANSMISSION_FILE:
		//  Case 5: Use transmission function read in from file
		case NEW_BB:
		//  Case 6: Use transmission function found for new breadboard (May 2017)
		case ARRAY_INPUT:
		//  Case 7: Use transmission function found for new breadboard (May 2017)
		case PIXL_FM_OPTIC:
		//  For all these cases, return an interpolated value determined by the optic transmission function in the arrays
            if( m_vectorData.vectorDataX.size() <= 0 ) return 0;
			float energyUpperLimit = m_vectorData.vectorDataX.at(m_vectorData.vectorDataX.size() - 1);		// Use the last vector element to get the upper energy limit
			float energyLowerLimit = m_vectorData.vectorDataX.at(0);		// Use the first vector element to get the lower energy limit
			float energyKeV = energy / 1000;		// Convert to KeV from eV
			if( m_type == ARRAY_INPUT || m_type == PIXL_FM_OPTIC_OLD || m_type == PIXL_FM_OPTIC ) energyKeV = energy;		// Array input is in eV
			// Interp only if we are still in the range of the optic data, otherwise return 0
			if (energyKeV >= energyLowerLimit && energyKeV <= energyUpperLimit)
			{
//				cout << energyKeV << "\t";
                float interpResult = splint( m_vectorData.vectorDataX, m_vectorData.vectorDataY, m_vectorData.vectorDataD, energyKeV);
//				float interpResult = interp(energyKeV, m_vectorData.vectorDataX, m_vectorData.vectorDataY);
//				cout << setprecision(9);
//				cout << interpResult << endl;
				return interpResult > 0 ? interpResult : 0;
			}
			else if( energyKeV < energyLowerLimit )
			{
                return m_vectorData.vectorDataY.at(0);
            }
			else
			{
				return 0;
			}
	}

    return 0;
}

// Parses through the file to find how many lines it has, returns the integer count
int XrayOptic::GetLineCount(string fileName)
{
	ifstream file(fileName.c_str());

	// If we couldn't open the input file stream for reading
    if (!file)
    {
        // Print an error and exit
        cout << "**** Can't open optic transmission file: " << fileName << endl;
        throw exception();
    }

	int lineCount = 0;
	// Get the length of the file in total lines
    while( file )
	{
        string line;
        getline(file, line);
        if( line.length() <= 0 ) break;
        // Make a stream for the line
        istringstream in(line);
        // Load line contents into two floats to check validity
        float x, y;
        in >> x >> y;
        if( ! in ) {
            cout << "**** Error reading optic file, line " << lineCount + 1 << endl;
            throw exception();
        }
        lineCount++;

        if (file.eof())	break;
	}
	file.close();

	return lineCount;
}

/*
void PrintVector(std::vector<float> vectorData)
{
	using namespace std;
	// Load values into vector for storage. Even indicies store X column, Odd indicies store Y column, print them to validate the data
	int linePosition = 0;
	int lineCount = GetLineCount(fileName);

	for (int i = 0; i < lineCount; i++)
	{
		try
		{
			cout << (i + 1) << "\t";					// Print line number as seen in the file
			cout << vectorData.at(i) << endl;
		}

		catch (exception& e)
		{
			cout << "Element: " << i << ", index beyond vector dimensions" << endl;
		}
	}
};
*/
// Reads the file and loads contents into floats to be printed
// Need to store the contents in something, and write something that interpolates
// Want to have a resolution selector, a parameter or argument that chooses how many points get interpolated
int XrayOptic::ParseFile(string fileName)
{
	// Open stream for file to be read into
	ifstream file(fileName.c_str());

    // If we couldn't open the input file stream for reading
    if (!file)
    {
        // Print an error and exit
        cout << "*** Optic file not found or can't open" << endl;
        return 0;
    }

    // Parse through the file line by line and store the contents in float vectors
    int lineNumber = 0;
	while (file)
    {
		// Read file line by line
        string line;
		getline(file, line);
        if( line.length() <= 0 ) break;

		// Make a stream for the line
		istringstream in(line);

		// Load line contents into two floats and print
		float x, y;
		in >> x >> y;

		// Load values into vectors for storage
		m_vectorData.vectorDataX.at(lineNumber) = x;		// Store X
		m_vectorData.vectorDataY.at(lineNumber) = y;		// Store Y
//	cout << " x " << m_vectorData.vectorDataX.at(lineNumber) << "   y " << m_vectorData.vectorDataY.at(lineNumber) << endl;

		lineNumber++;

		if (file.eof())	break;							// Catch the end of file before it goes through the while loop one extra time
    }

	//  Derivative data all zeros for linear interpolation
    vector <float> temp_vec( m_vectorData.vectorDataX.size(), 0 );
    m_vectorData.vectorDataD = temp_vec;

    return 0;
}

/*
void InterpPrintTest()
{
	using namespace std;

	const int arrX[] = { 0, 1, 2, 3, 4 };
	vector<float> vectorX (arrX, arrX + sizeof(arrX) / sizeof(arrX[0]) );
	const int arrY[] = { 0, 10, 20, 30, 400 };
	vector<float> vectorY (arrY, arrY + sizeof(arrY) / sizeof(arrY[0]) );

	cout << "Interp Test" << endl;
	cout << "Vector X:\t";
	for (int i = 0; i < 5; i++)
	{
		cout << vectorX.at(i) << "\t";
	}
	cout << "\nVector Y:\t";
	for (int i = 0; i < 5; i++)
	{
		cout << vectorY.at(i) << "\t";
	}
	cout << "\n\nInterpolated:\t";

	float floatArray[] = { .5, 1.5, 2.5, 3.5, 4.5 };
	for (int i = 0; i < 5; i++)
	{
		cout << floatArray[i] << "\t";
	}
	cout << "\n\t";
	for (int i = 0; i < 5; i++)
	{
		cout << "\t" << interp(floatArray[i], vectorX, vectorY) << " ";
	}
	cout << "\n" << endl;
}
*/

std::string XrayOptic::toString() const
{
    std::ostringstream os;
    os << "XrayOptic:" << endl;
    os << "  m_optic_file: " << m_optic_file << endl;
    os << "  m_type: ";
    switch(m_type)
    {
    case NO_OPTIC:
        os << "NO_OPTIC";
        break;
    case BOXCAR:
        os << "BOXCAR";
        break;
    case PIXL:
        os << "PIXL";
        break;
    case TRANSMISSION_FILE:
        os << "TRANSMISSION_FILE";
        break;
    case NEW_BB:
        os << "NEW_BB";
        break;
    case ARRAY_INPUT:
        os << "ARRAY_INPUT";
        break;
    case PIXL_FM_OPTIC:
        os << "PIXL_FM";
        break;
    }
    os << endl;

    os << "  m_centerEnergy: " << m_centerEnergy << endl;
    os << "  m_bandwidth: " << m_bandwidth << endl;
    os << "  m_maxTransmission: " << m_maxTransmission << endl;
    os << "  m_defaultFlag: " << m_defaultFlag << endl;
    os << "  m_vectorData X size: " << m_vectorData.vectorDataX.size() << endl;
    os << "  m_vectorData Y size: " << m_vectorData.vectorDataY.size() << endl;
    os << "  m_vectorData D size: " << m_vectorData.vectorDataD.size() << endl;

    return os.str();
}
