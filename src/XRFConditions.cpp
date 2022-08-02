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

#include "XRFconditions.h"
#include <sstream>
#include <string>

using namespace std;


string toString(const XRFconditions &cond)
{
    ostringstream os;
    os << "XRFconditions:" << endl;
    os << "  source:" << endl << cond.source.toString() << endl;
    os << "  filter:" << endl << cond.filter.toString() << endl;
    os << "  optic:" << endl << cond.optic.toString() << endl;
    os << "  dust_on_optic:" << endl << cond.dust_on_optic.toString() << endl;
    os << "  incidentPath:" << endl << cond.incidentPath.toString() << endl;

    os << "  solidAngleSource:" << cond.solidAngleSource << endl;
    os << "  excitAngle:" << cond.excitAngle << endl;
    os << "  excitCosecant:" << cond.excitCosecant << endl;
    os << "  geometryFactor:" << cond.geometryFactor << endl;

    os << "  dust_on_specimen:" << endl << cond.dust_on_specimen.toString() << endl;
    os << "  window:" << endl << cond.window.toString() << endl;

    os << "  emergAngle:" << cond.emergAngle << endl;
    os << "  emergCosecant:" << cond.emergCosecant << endl;

    os << "  emergentPath:" << endl << cond.emergentPath.toString() << endl;
    os << "  dust_on_detector:" << endl << cond.dust_on_detector.toString() << endl;

    os << "  solidAngleDetector:" << cond.solidAngleDetector << endl;

    os << "  detector:" << endl << cond.detector.toString() << endl;

    os << "  eMin:" << cond.eMin << endl;

    return os.str();
}
