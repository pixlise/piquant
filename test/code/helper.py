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

import os
import subprocess
import math

output_root_path = '/build/test/output'

def make_output_path(to_file):
    return os.path.join(output_root_path, to_file)

def make_cmd(testSelf, cmd, files, elements, output, params):
    cmd = [ testSelf.piquant, cmd, testSelf.config_file, testSelf.calibration_file, files, elements, make_output_path(output) ]
    if len(params) > 0:
        cmd.extend(params.split(' '))
    return cmd

def run_piquant(testSelf, cmdlist):
    #print("=========================\nRunning Piquant\n\n")
    ran = subprocess.run(cmdlist, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    stdout = ran.stdout.decode('utf-8')
    stderr = ran.stderr.decode('utf-8')

    # If error code returned, show stdout at least
    if ran.returncode != 0:
        print("CMD \"{}\" returned code: {}".format(' '.join(cmdlist), ran.returncode))
        print("STDOUT: \"{}\"".format(stdout))
        print("STDERR: \"{}\"\n".format(stderr))

    testSelf.assertEqual(ran.returncode, 0)

    #if ran.returncode != 0:
    #    raise ValueError('Running: "{}" result was: {}', cmd, ran.returncode)
    #print("=========================\n")

    # Print the output file
    logfile = ''

    logpath = cmdlist[1]+'_log.txt'
    if os.path.exists(logpath):
        with open(logpath) as f:
            logfile = f.read()

    return [logfile, stdout, stderr]


def compare_outputs(testSelf, out_file, expected_file, outputLogs):
    exppath = './test-data/expected-output/{}'.format(expected_file)

    compare_outputs_exppath(testSelf, out_file, exppath, outputLogs)


def compare_outputs_exppath(testSelf, out_file, exppath, outputLogs):
    outpath = make_output_path(out_file)

    # Pull in files, sort rows alphabetically, then diff the 2 lists
    with(open(outpath)) as f:
        outfile_lines = sorted(f.readlines())
    with(open(exppath)) as f:
        expfile_lines = sorted(f.readlines())

    # Now loop through them
    equal = True
    for c in range(len(outfile_lines)):
        if outfile_lines[c] != expfile_lines[c]:
            equal = False
            print('\n\nFOUND NON-MATCHING LINES:')
            print('Output:   '+outfile_lines[c])
            print('Expected: '+expfile_lines[c]+'\n')
            break

    if not equal:
        print_logs(outputLogs)

    testSelf.assertEqual(len(outfile_lines), len(expfile_lines))
    testSelf.assertEqual(equal, True)

def print_logs(outputLogs):
    print('\nDUMPING PIQUANT LOG:')
    print(outputLogs[0])
    print('\nDUMPING PIQUANT STDOUT:')
    print(outputLogs[1])
    print('\nDUMPING PIQUANT STDERR:')
    print(outputLogs[2])


def firstDiffIdx(a, b):
    shortestLen = min(len(a), len(b))
    for c in range(shortestLen):
        if a[c] != b[c]:
            return c
    # Otherwise diff must be due to their length
    if len(a) == len(b):
        return -1
    return shortestLen


def isVersionDiff(out_line, exp_line,):
    diffIdx = firstDiffIdx(out_line, exp_line)

    PIQUANT_VER_START = 'PIQUANT '

    # Find where PIQUANT version is
    outVerIdx = out_line.find(PIQUANT_VER_START)
    expVerIdx = exp_line.find(PIQUANT_VER_START)

    # If they're in the same place, and before the difference point, it's probably just the version diff
    if outVerIdx == expVerIdx and outVerIdx > -1 and diffIdx > (outVerIdx+len(PIQUANT_VER_START)):
        return True
    return False


def getNumValue(str):
    try:
        return int(str)
    except ValueError:
        try:
            return float(str)
        except ValueError:
            return None


def isDataRowDiffWithinLimit(out_line, exp_line, maxVariancePercentDiff):
    dbg = False
    diffIdx = firstDiffIdx(out_line, exp_line)

    out_items = out_line.split(',')
    exp_items = exp_line.split(',')

    if len(out_items) != len(exp_items):
        # Different item count, can't argue there!
        if dbg:
            print("item count diff: '{}' vs '{}'".format(len(exp_items), len(out_items)))
        return False, None

    # Check each field
    for c in range(len(out_items)):
        out_item = out_items[c].strip()
        exp_item = exp_items[c].strip()

        # Are they both numbers
        out_num = getNumValue(out_item)
        exp_num = getNumValue(exp_item)

        # If both are NOT numbers, just string compare
        if out_num is None and exp_num is None:
            if out_item != exp_item:
                if dbg:
                    print("str diff: '{}' vs '{}'".format(out_item, exp_item))
                return False, None
            # else it's a match so we continue on...
        else:
            # If either of them is NOT a number, it's a string difference, so
            # cant check numerical difference is within limit
            out_weird = out_num is None or math.isnan(out_num) or math.isinf(out_num)
            exp_weird = exp_num is None or math.isnan(exp_num) or math.isinf(exp_num)

            if not out_weird and not exp_weird:
                # Both are numbers that can be compared
                if exp_num != out_num:
                    if exp_num == 0: # Can't calculate a diff variance
                        if dbg:
                            print("cant calc variance: {} vs {}".format(exp_num, out_num))
                        return False, None

                    # Check their % difference
                    variance = abs((out_num-exp_num)/exp_num)
                    if variance > maxVariancePercentDiff:
                        if dbg:
                            print("differs with variance: {} vs {} = {}".format(exp_num, out_num, variance))
                        return False, variance
            elif out_weird != exp_weird:
                if dbg:
                    print("weirdness: out={}({}), exp={}({})".format(out_item, out_num, exp_item, exp_num))
                return False, None

    return True, None

HEADER_COMPARE = 0
HEADER_IGNORE_PIQUANT_VERSION_DIFF = 1
HEADER_IGNORE_ALL = 2

def getCSVDifferences(outfile_lines, expfile_lines, firstDataRowIdx, headerLineCompare, maxVariancePercentDiff, ignoreRows):
    numRows = len(outfile_lines)
    if len(expfile_lines) < numRows:
        numRows = len(expfile_lines)

    skippedDiffLineIdxs = []
    foundDiffLineIdxs = []
    maxVariance = 0

    for c in range(numRows):
        # If there is a difference
        if outfile_lines[c] != expfile_lines[c]:
            out_line = outfile_lines[c]
            exp_line = expfile_lines[c]

            # If it's on the first line, it might be the piquant version number
            if c == 0:
                if headerLineCompare == HEADER_COMPARE:
                    # We were told to compare the header line, and it differs!
                    foundDiffLineIdxs.append(c)
                elif headerLineCompare == HEADER_IGNORE_ALL:
                    # Header line ignore, we add it to the skipped list
                    skippedDiffLineIdxs.append(c)
                elif headerLineCompare == HEADER_IGNORE_PIQUANT_VERSION_DIFF and isVersionDiff(out_line, exp_line):
                    # Told to check version, and that differed, so add to skip list
                    skippedDiffLineIdxs.append(c)
                else: # Otherwise it's a difference!
                    foundDiffLineIdxs.append(c)

            # If it's a data row, check each field, might be < rounding error
            elif c >= firstDataRowIdx:
                # If it's an ignore row...
                if c in ignoreRows:
                    skippedDiffLineIdxs.append(c)
                else:
                    ok, variance = isDataRowDiffWithinLimit(out_line, exp_line, maxVariancePercentDiff)
                    if variance is not None and variance > maxVariance:
                        maxVariance = variance
                    #print("variance: {} for {}/{}".format(maxVariance, out_line, exp_line))
                    if ok:
                        skippedDiffLineIdxs.append(c)
                    else:
                        #print("Found difference: '{}' vs '{}'. Variance: {}".format(out_line, exp_line, variance))
                        foundDiffLineIdxs.append(c)
            else:
                foundDiffLineIdxs.append(c)

    return skippedDiffLineIdxs, foundDiffLineIdxs, maxVariance

def compare_output_csvs(testSelf, outputPath, expPath, firstDataRowIdx, headerLineCompare, maxVariancePercentDiff, ignoreRows, outputLogs):
    with(open(outputPath)) as f:
        outfile_lines = f.read().splitlines()
    with(open(expPath)) as f:
        expfile_lines = f.read().splitlines()

    #print('"{}"'.format(outfile_lines[0]))
    #print('"{}"'.format(expfile_lines[0]))

    #if len(outfile_lines) != len(expfile_lines):
    #    print('output files is '+str(len(outfile_lines))+' lines long, expected '+str(len(expfile_lines))+' lines')
    testSelf.assertEqual(len(outfile_lines), len(expfile_lines))

    skipped, diffs, maxVariance = getCSVDifferences(outfile_lines, expfile_lines, firstDataRowIdx, headerLineCompare, maxVariancePercentDiff, ignoreRows)

    lines = 10
    if len(diffs) > 0:
        print('\n\nFOUND NON-MATCHING LINES, SHOWING FIRST {}:'.format(lines))
        for c in diffs:
            print('[{}] Output:   {}'.format(c, outfile_lines[c]))
            print(' {}  Expected: {}'.format(' '*len('{}'.format(c)), expfile_lines[c]))
            if c >= 10:
                break

        print_logs(outputLogs)

    if maxVariance > 0:
        print("Max variance seen between {} and {}: {:.10f}".format(outputPath, expPath, maxVariance))
    testSelf.assertEqual(len(diffs), 0)



