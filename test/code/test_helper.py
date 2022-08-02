# Copyright (c) 2018-2022 California Institute of Technology (“Caltech”) and
# University of Washington. U.S. Government sponsorship acknowledged.
# All rights reserved.
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# * Neither the name of Caltech nor its operating division, the Jet Propulsion
#   Laboratory, nor the names of its contributors may be used to endorse or
#   promote products derived from this software without specific prior written
#   permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

import unittest

from helper import *


class TesterTester(unittest.TestCase):
    def test_getCSVDifferences(self):
        outfile = [ "some PIQUANT 1.2.1 file", "A, B, C", "100, 200, 300" ]
        expfile = [ "some PIQUANT 1.2.4 file", "A, B, E", "105, 200, 280" ]
        self.assertEqual(getCSVDifferences(outfile, expfile, 2, HEADER_COMPARE, 0, []), ([], [0, 1, 2], 0.047619047619047616))
        self.assertEqual(getCSVDifferences(outfile, expfile, 2, HEADER_IGNORE_ALL, 0, []), ([0], [1, 2], 0.047619047619047616))
        self.assertEqual(getCSVDifferences(outfile, expfile, 2, HEADER_IGNORE_PIQUANT_VERSION_DIFF, 0.048, []), ([0], [1, 2], 0.07142857142857142))
        self.assertEqual(getCSVDifferences(outfile, expfile, 2, HEADER_IGNORE_PIQUANT_VERSION_DIFF, 0.072, []), ([0,2], [1], 0))
        self.assertEqual(getCSVDifferences(outfile, expfile, 2, HEADER_IGNORE_PIQUANT_VERSION_DIFF, 0.072, []), ([0,2], [1], 0))
        self.assertEqual(getCSVDifferences(outfile, expfile, 2, HEADER_IGNORE_PIQUANT_VERSION_DIFF, 0, [2]), ([0,2], [1], 0))

    def test_firstDiffIdx(self):
        self.assertEqual(firstDiffIdx("aa", "ba"), 0)
        self.assertEqual(firstDiffIdx("aa", "ab"), 1)
        self.assertEqual(firstDiffIdx("aaa", "aab"), 2)
        self.assertEqual(firstDiffIdx("aaa", "aaab"), 3)
        self.assertEqual(firstDiffIdx("aaaba", "aaab"), 4)
        self.assertEqual(firstDiffIdx("this is a string", "this isn't the same"), 7)
        self.assertEqual(firstDiffIdx("same string", "same string"), -1)

    def test_isVersionDiff(self):
        self.assertEqual(isVersionDiff("Something   PIQUANT 2.6.310-compiling  path/", "Something   PIQUANT 2.6.311-compiling  path/"), True)
        self.assertEqual(isVersionDiff("Something   PIQUANT 2.6.310-compiling  path/", "Something   PIQUANT 2.6.0  path/"), True)

        self.assertEqual(isVersionDiff("Something   PIQUANT 2.6.310-compiling  path/", "Something else   PIQUANT 2.6.0  path/"), False)
        self.assertEqual(isVersionDiff("Something   PIQUANT 2.6.310-compiling  path/", "Totally, different, line"), False)

    def test_getNumValue(self):
        self.assertEqual(getNumValue("hello"), None)
        self.assertEqual(getNumValue("1.2"), 1.2)
        self.assertEqual(getNumValue("12"), 12)
        self.assertEqual(getNumValue("12.0"), 12.0)
        self.assertEqual(getNumValue("1.2e03"), 1200)
        self.assertEqual(getNumValue("1.2e003"), 1200)
        self.assertEqual(getNumValue("1.2345e006"), 1234500.0)
        self.assertEqual(getNumValue("1.2345e002"), 123.45)
        self.assertEqual(getNumValue("1.2345e2"), 123.45)

    def test_isDataRowDiffWithinLimit(self):
        self.assertEqual(isDataRowDiffWithinLimit("100, 1.5, 1.2e02", "100, 1.5, 1.2e02", 0.00), (True, None))
        self.assertEqual(isDataRowDiffWithinLimit("100, 1.5, 1.2e02", "100, 1.5, 1.2e02", 0.01), (True, None))

        self.assertEqual(isDataRowDiffWithinLimit("100.1, 1.5, 1.2e02", "100, 1.5, 1.2e02", 0.01), (True, None))
        self.assertEqual(isDataRowDiffWithinLimit("101.0, 1.5, 1.2e02", "100, 1.5, 1.2e02", 0.01), (True, None))
        self.assertEqual(isDataRowDiffWithinLimit("101.1, 1.5, 1.2e02", "100, 1.5, 1.2e02", 0.01), (False, 0.010999999999999944))
        self.assertEqual(isDataRowDiffWithinLimit("101.1, 1.5, 1.2e02", "100, 1.5, 1.2e02", 0.02), (True, None))
        self.assertEqual(isDataRowDiffWithinLimit("102.0, 1.5, 1.2e02", "100, 1.5, 1.2e02", 0.02), (True, None))
        self.assertEqual(isDataRowDiffWithinLimit("102.001, 1.5, 1.2e02", "100, 1.5, 1.2e02", 0.02), (False, 0.02001000000000005))

        self.assertEqual(isDataRowDiffWithinLimit("100.1, 1.6, 1.2e02", "100, 1.5, 1.2e02", 0.06), (False, 0.06666666666666672))
        self.assertEqual(isDataRowDiffWithinLimit("100.1, 1.6, 1.2e02", "100, 1.5, 1.2e02", 0.07), (True, None))
        self.assertEqual(isDataRowDiffWithinLimit("100.1, 1.7, 1.2e02", "100, 1.5, 1.2e02", 0.07), (False, 0.1333333333333333))
        self.assertEqual(isDataRowDiffWithinLimit("100.1, 1.7, 1.2e02", "100, 1.5, 1.2e02", 0.13), (False, 0.1333333333333333))
        self.assertEqual(isDataRowDiffWithinLimit("100.1, 1.7, 1.2e02", "100, 1.5, 1.2e02", 0.14), (True, None))

        self.assertEqual(isDataRowDiffWithinLimit("100, 1.5, 1.234e02", "100, 1.5, 1.234e02", 0.00), (True, None))
        self.assertEqual(isDataRowDiffWithinLimit("100, 1.5, 1.236e02", "100, 1.5, 1.234e02", 0.0015), (False, 0.0016207455429496647))
        self.assertEqual(isDataRowDiffWithinLimit("100, 1.5, 1.236e02", "100, 1.5, 1.234e02", 0.0017), (True, None))

        self.assertEqual(isDataRowDiffWithinLimit("100, 1.5, 1.236e02", "100, Hello, 1.234e02", 0.0017), (False, None))
        self.assertEqual(isDataRowDiffWithinLimit("100, 1.5, 1.236e02", "100, NaN, 1.234e02", 0.0017), (False, None))
        self.assertEqual(isDataRowDiffWithinLimit("100, 1.5, 1.236e02", "100, inf, 1.234e02", 0.0017), (False, None))

        self.assertEqual(isDataRowDiffWithinLimit("100.1, 1.5, 1.2e02", "100, 1.5, 1.2e02, 88", 0.01), (False, None))

        self.assertEqual(isDataRowDiffWithinLimit("100.1, -1.5, 1.2e02", "100, 1.5, 1.2e02", 0.01), (False, 2.0))
        self.assertEqual(isDataRowDiffWithinLimit("100.1, 1.5, 1.2e02", "100, -1.5, 1.2e02", 0.01), (False, 2.0))

        self.assertEqual(isDataRowDiffWithinLimit("-100, 1.5, 1.2e02", "-100, 1.5, 1.2e02", 0.00), (True, None))
        self.assertEqual(isDataRowDiffWithinLimit("-100, 1.5, 1.2e02", "-100, 1.5, 1.2e02", 0.01), (True, None))

        self.assertEqual(isDataRowDiffWithinLimit("-100.1, 1.5, 1.2e02", "-100, 1.5, 1.2e02", 0.01), (True, None))
        self.assertEqual(isDataRowDiffWithinLimit("-101.0, 1.5, 1.2e02", "-100, 1.5, 1.2e02", 0.01), (True, None))
        self.assertEqual(isDataRowDiffWithinLimit("-101.1, 1.5, 1.2e02", "-100, 1.5, 1.2e02", 0.01), (False, 0.010999999999999944))
        self.assertEqual(isDataRowDiffWithinLimit("101.1, 1.5, 1.2e02", "-100, 1.5, 1.2e02", 0.01), (False, 2.011))
        self.assertEqual(isDataRowDiffWithinLimit("-101.1, 1.5, 1.2e02", "100, 1.5, 1.2e02", 0.01), (False, 2.011))

        self.assertEqual(isDataRowDiffWithinLimit("101, 1.5, 1.2e02", "100, 0 1.2e02", 0.01), (False, None))

        # Catches first one (100-105)/105
        self.assertEqual(isDataRowDiffWithinLimit("100, 200, 300", "105, 200, 30", 0.0), (False, 0.047619047619047616))
        # Skips first one, finds (300-30)/30
        self.assertEqual(isDataRowDiffWithinLimit("100, 200, 300", "105, 200, 30", 0.05), (False, 9))

        # Some real examples
        self.assertEqual(isDataRowDiffWithinLimit(
            "0.22, 3.43318e-43, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.43318e-43, 0, 0, 0, 0",
            "0.22, 3.43318e-043, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.43318e-043, 0, 0, 0, 0",
            0.0), (True, None))

        self.assertEqual(isDataRowDiffWithinLimit(
            "13.8145, 92, 93.4549, 7, 92.9913, 9.69536, -1.45494, 93.4549, 92.9913, 92.9913, 92.9913, 92.9913, 92.9913, 92.9913, 92.9913, 92.9913, 92.9913, 92.9913",
            "13.8145, 92, 93.4549, 7, 92.9913, 9.69536, -1.45493, 93.4549, 92.9913, 92.9913, 92.9913, 92.9913, 92.9913, 92.9913, 92.9913, 92.9913, 92.9913, 92.9913",
            0.000007), (True, None))

        self.assertEqual(isDataRowDiffWithinLimit(
            "-3.63827e-05, 0, 0, 0, 0, 1.41421, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0",
            "-3.63884e-005, 0, 0, 0, 0, 1.41421, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0",
            0.004), (True, None))

        self.assertEqual(isDataRowDiffWithinLimit(
            "V, , , el, 0.0416%, 3.4%, 0.0,   1.00, 1.0488, 0.6%, 232588.7, 23",
            "V, , , el, 0.0416%, 3.4%, 0.0,   1.00, 1.0488, 0.6%, 232588.1, 23",
            0.00000257), (False, 2.579667661440206e-06))

        self.assertEqual(isDataRowDiffWithinLimit(
            "V, , , el, 0.0416%, 3.4%, 0.0,   1.00, 1.0488, 0.6%, 232588.7, 23",
            "V, , , el, 0.0416%, 3.4%, 0.0,   1.00, 1.0488, 0.6%, 232588.1, 23",
            0.000003), (True, None))


        self.assertEqual(isDataRowDiffWithinLimit(
            "252, 0.0000, 0.0000, 0.0000, 0.0000, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8561, 4.41, -nan, 0.0, 7.9978, 155, 5, Normal_A_0634651638_000001C5_000252.msa, 8561, 9635, 0, 0",
            "252, 0.0000, 0.0000, 0.0000, 0.0000, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8561, 4.41, nan, 0.0, 7.9978, 155, 5, Normal_A_0634651638_000001C5_000252.msa, 8561, 9635, 0, 0",
            0.000003), (True, None))

        self.assertEqual(isDataRowDiffWithinLimit(
            "352, 0.0622, 0.0082, 0.0164, 0.0893, 34.1, 3.9, 7.4, 29.5, 44.9, 440.1, 244.2, 67.2, 8193, 4.48, 0.37, -96.4, 8.1376, 155, 40, Normal_B_0634653145_000001C5_000352.msa, 8192, 9093, 0, 0",
            "352, 0.0622, 0.0082, 0.0165, 0.0893, 34.1, 3.9, 7.5, 29.5, 44.9, 440.5, 243.6, 67.2, 8193, 4.48, 0.37, -96.2, 8.1373, 155, 13, Normal_B_0634653145_000001C5_000352.msa, 8192, 9093, 0, 0",
            2.08), (True, None))
