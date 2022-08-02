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


#include <chrono>

class time_code
{
public:
    time_code(const char *name, bool log = true)
    {
        _name = name;
        _log = log;
        _start = _last = std::chrono::high_resolution_clock::now();
    }

    ~time_code()
    {
        split(0);
    }

    void split(const char *label)
    {
        auto end = std::chrono::high_resolution_clock::now();

        if(_log)
        {
            if(label)
            {
                auto elapsed = elapsedSince(true);
                cout << "  >>" << _name << "(" << label << ") took: " << elapsed << " sec" << endl;
            }
            else
            {
                auto elapsed = elapsedSince(false);
                cout << "====================================" << endl;
                cout << _name << "(RUN) took: " << elapsed << " sec" << endl;
            }
        }

        _last = end;
    }

    double elapsedSince(bool sinceSplit)
    {
        auto from = _start;
        if(sinceSplit)
        {
            from = _last;
        }

        auto now = std::chrono::high_resolution_clock::now();
        return std::chrono::duration_cast<std::chrono::duration<double>>(now-from).count();
    }


protected:
    std::chrono::high_resolution_clock::time_point _start, _last;
    const char *_name;
    bool _log;
};