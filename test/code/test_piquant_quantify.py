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


class PiquantQuantifyTester(unittest.TestCase):
    config_file = './test-data/config/PIXL/Config_PIXL_FM_SurfaceOps_Rev1_Jul2021.msa'
    piquant = './Piquant'

    def test_quantify(self):
        cmd = [
            self.piquant,
            'quant',
            self.config_file,
            './test-data/PIQUANT_test_data_May2020/Input_files_PIQUANT_test_data_May2020/Calibrate_Master_ECF_new_BB_01_08_2020.csv',
            './test-data/PIQUANT_test_data_May2020/Input_files_PIQUANT_test_data_May2020/Calibration_box_BHVO-2G_28kV_230uA_03_28_2019_bulk_sum.msa',
            'Si_K K_K P_K Ca_K Ti_K Cr_K Mn_K Fe_K Sr_K Ar_I',
            make_output_path('output_quantify.csv')
            ]
        log = run_piquant(self, cmd)

        compare_output_csvs(self,
            make_output_path('output_quantify.csv'),
            './test-data/PIQUANT_test_data_May2020/ExpectedOutputs/exp_quantify.csv',
            2,
            HEADER_IGNORE_ALL,
            0.0,
            [],
            log)
