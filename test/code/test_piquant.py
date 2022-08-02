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


class PiquantTester(unittest.TestCase):
    config_file = './test-data/config/PIXL/Config_PIXL_FM_SurfaceOps_Rev1_Jul2021.msa'
    calibration_file = './test-data/config/PIXL/Calibration_PIXL_FM_SurfaceOps_5minECFs_Rev1_Jul2021.csv'
    piquant = './Piquant'

    def test_multiple_maps(self):
        cmd = make_cmd(self, 'map', './test-data/msa/files.txt', 'Fe,Ca,Ti,K', 'multi_map.csv', '-t,6')
        log = run_piquant(self, cmd)
        compare_outputs(self, 'multi_map.csv', 'multi_map.csv', log)

    def test_3PMC_map(self):
        cmd = make_cmd(self, 'map', './test-data/msa/6files.txt', 'Fe,Ca,Ti,K', '6map.csv', '')
        log = run_piquant(self, cmd)
        compare_outputs(self, '6map.csv', '6map.csv', log)

    # This test is does the same as test_3PMC_map but with the input source being a PIXLISE binary file. We expect the same output
    def test_using_pmcs_AB(self):
        cmd = make_cmd(self, 'map', './test-data/pixlise-datasets/list.pmcs', 'Fe,Ca,Ti,K', 'multi_pmc_map.csv', '-t,1')
        log = run_piquant(self, cmd)
        compare_outputs(self, 'multi_pmc_map.csv', '6map_pmcsfile_AB.csv', log)

    def test_using_pmcs_Combined(self):
        cmd = make_cmd(self, 'map', './test-data/pixlise-datasets/list_combined.pmcs', 'Fe,Ca,Ti,K', 'multi_pmc_map_combined.csv', '-t,1')
        log = run_piquant(self, cmd)
        compare_outputs(self, 'multi_pmc_map_combined.csv', '6map_pmcsfile_combined.csv', log)

    def test_single_map(self):
        cmd = make_cmd(self, 'map', './test-data/msa/1file.txt', 'Fe,Ca,Ti,K', '1map.csv', '')
        log = run_piquant(self, cmd)
        compare_outputs(self, '1map.csv', '1map.csv', log)

    def test_bulksum(self):
        cmd = make_cmd(self, 'sum', './test-data/msa/sum.txt', 'Fe,Ca,Ti,K', 'bulksummax.msa', '')
        log = run_piquant(self, cmd)
        compare_outputs(self, 'bulksummax.msa', 'bulksummax.msa', log)
