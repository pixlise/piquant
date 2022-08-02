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
import os
import subprocess


class PiquantVersionTester(unittest.TestCase):
    piquant = './Piquant'

    # Makes sure we're testing the version that was built from the code base we thought...
    # Code version is set in CMakeLists.txt, we check that the exe reports the same version as is set
    # in the file
    def test_aversion(self):
        ran = subprocess.run([self.piquant, 'version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Get the version number, ensure it's what we're expecting
        version = ran.stdout.decode("utf-8").rstrip()

        # We used to get the version number from CMakeLists.txt but it now looks for a CI build env var.
        # Locally now, we default to 0.0.0, so here if we can't get the build env var, we compare to that
        expectedVersion = '0.0.0'
        if 'BUILD_VERSION' in os.environ and len(os.environ['BUILD_VERSION']) > 0:
            expectedVersion = os.environ['BUILD_VERSION']

        # Get the branch we're on, this is passed into our docker container as an env variable (see root test.sh)
        expectedBranch = os.environ['GIT_BRANCH_ENV']
        if len(expectedBranch) <= 0:
            expectedBranch = os.environ['GIT_BRANCH_GIT']

        # Expecting to get back version-branch
        expectedStr = '{}-{}'.format(expectedVersion, expectedBranch)
        self.assertEqual(version, expectedStr)

