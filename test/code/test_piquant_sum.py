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


class PiquantSumTester(unittest.TestCase):
    config_file = './test-data/config/PIXL/Config_PIXL_FM_SurfaceOps_Rev1_Jul2021.msa'
    calibration_file = './test-data/config/PIXL/Calibration_PIXL_FM_SurfaceOps_5minECFs_Rev1_Jul2021.csv'
    piquant = './Piquant'

    def test_sum_puck(self):
        cmd = make_cmd(self, 'sum', './test-data/PIQUANT_test_data_May2020/Data_files_Pucks_20200210_200210210814/msa_file_list_DetA_BHVO.txt', 'Si_K K_K P_K Ca_K Ti_K Cr_K Mn_K Fe_K Sr_K Ar_I', 'output_sum_puck.csv', '')
        log = run_piquant(self, cmd)
        compare_output_csvs(self,
            make_output_path('output_sum_puck.csv'),
            './test-data/PIQUANT_test_data_May2020/ExpectedOutputs/exp_sum_puck.csv',
            2,
            HEADER_IGNORE_ALL,
            0.0,
            [],
            log)

# Not a real test, just used this to generate a bulk sum for crosshair dataset in test-data repo
#    def test_sum_crosshair(self):
#        cmd = make_cmd(self, 'sum', './test-data/PIQUANT_test_data_May2020/Data_files_Crosshair_20200211_200211010146/msa_file_list.txt', 'Ti_K Cr_K Fe_K Ni_K Ar_I', 'output_sum_crosshair.msa', '')
#        log = run_piquant(self, cmd)
#        compare_output_csvs(self,
#            make_output_path('output_sum_crosshair.msa'),
#            './test-data/PIQUANT_test_data_May2020/Output_files_PIQUANT_test_data_V260_May2020/Bulk_sum_and_max_value_Plot_PIQUANT_test_data_May2020.csv',
#            2,
#            HEADER_IGNORE_ALL,
#            0.004,
#            [],
#            log)
