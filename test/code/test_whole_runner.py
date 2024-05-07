import unittest
import os
import subprocess
import boto3
import sys
import json

if not sys.warnoptions:
    import warnings

from helper import *

def prep_cloud(bucket, jobID, envName, nodeName):
    warnings.filterwarnings("ignore", category=ResourceWarning, message="unclosed.*<ssl.SSLSocket.*>")

    # Before we do anything, reset the bucket
    session = boto3.Session()

    # NOTE: if this is failing, set AWS_DEFAULT_REGION env variable!
    s3 = session.resource('s3')

    try:
        logGroup = '/run-piquant/'+envName
        print('\n\nDeleting log group before test run: '+logGroup)
        # Delete log group/stream
        cw = boto3.client('logs')
        cw.delete_log_group(logGroupName=logGroup)
        #cw.delete_log_stream(logGroupName='/run-piquant/'+envName+'/'+jobID, logStreamName='job-'+jobID+'-node-'+nodeName)
    except Exception as e:
        print('\nLog group probably exist yet when running this test, this is OK. Here is the detailed error message:')
        print(e)

    s3Bucket = s3.Bucket(bucket)
    s3Bucket.objects.all().delete()

    bucket_start_contents = [
        'Jobs/'+jobID+'/node00012.pmcs',
        'Jobs/'+jobID+'/node00012-roi.pmcs',
        'Jobs/'+jobID+'/node00000.pmcs',
        'PiquantConfig/PIXL/v6/Config_PIXL_FM_SurfaceOps_Rev1_Jul2021.msa',
        'PiquantConfig/PIXL/v6/Calibration_PIXL_FM_SurfaceOps_5minECFs_Rev1_Jul2021.csv',
        'PiquantConfig/PIXL/v6/config.json',
        'Datasets/test-rtt-123/5x11dataset.bin'
    ]

    # Upload our test files to the bucket
    local = './test-data/s3-test-bucket-contents/'
    for f in bucket_start_contents:
        path = local+f
        s3Bucket.upload_file(path, f)


class WholeRunnerTester(unittest.TestCase):
    def makeEnvVars(self, paramJson):
        envVars = {
            "QUANT_PARAMS": paramJson,
            "AWS_ACCESS_KEY_ID": os.getenv('AWS_ACCESS_KEY_ID'),
            "AWS_SECRET_ACCESS_KEY": os.getenv('AWS_SECRET_ACCESS_KEY'),
            "AWS_DEFAULT_REGION": os.getenv('AWS_DEFAULT_REGION')
        }
        return envVars

    def test_run_piquant_s3(self):
# Only here because for some reason they're not bothering fixing this warning...
# https://github.com/boto/boto3/issues/454
        warnings.filterwarnings("ignore", category=ResourceWarning, message="unclosed.*<ssl.SSLSocket.*>")

        bucket = 'test-piquant'
        jobID = 'test-job-123'
        envName = 'piquant-unit-test'
        nodeName = 'node00012'

        prep_cloud(bucket, jobID, envName, nodeName)

        paramObj = {
            'runtimeEnv': envName,
            'datasetsBucket': bucket,
            'jobBucket': bucket,
            'configBucket': bucket,
            'jobId': jobID,
            'jobsPath': "Jobs",
            'datasetPath': "Datasets/test-rtt-123",
            'detectorConfig': "PiquantConfig/PIXL/v6",
            'elements': ["Fe","Ca","Ti","K"],
            'parameters': "",
            'pmcListName': nodeName+'.pmcs',
            'command': "map"
        }

        paramJson = json.dumps(paramObj, indent=4, sort_keys=True)
        envVars = self.makeEnvVars(paramJson)

        ran = subprocess.run(['./PiquantRunner'],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=envVars)

        if ran.returncode != 0:
            print(ran.stdout.decode('utf-8'))
            print(ran.stderr.decode('utf-8'))

        self.assertEqual(ran.returncode, 0)

        # Check that the output made it to S3, where we expect, and that it's valid
        s3 = boto3.resource('s3')

        s3.meta.client.download_file(bucket, 'Jobs/'+jobID+'/output/node00012.pmcs_result.csv', make_output_path('downloaded_output.csv'))

        compare_outputs(self, 'downloaded_output.csv', '6map_pmcsfile_combined.csv', ['', '', ''])

# Same as the above test but uses the modified pmcs files which have ROI tags in them. We expect the output file to contain these
    def test_run_piquant_s3_roi(self):
# Only here because for some reason they're not bothering fixing this warning...
# https://github.com/boto/boto3/issues/454
        warnings.filterwarnings("ignore", category=ResourceWarning, message="unclosed.*<ssl.SSLSocket.*>")

        bucket = 'test-piquant'
        jobID = 'test-job-124'
        envName = 'piquant-unit-test'
        nodeName = 'node00012-roi'

        prep_cloud(bucket, jobID, envName, nodeName)

        paramObj = {
            'runtimeEnv': envName,
            'datasetsBucket': bucket,
            'jobBucket': bucket,
            'configBucket': bucket,
            'jobId': jobID,
            'jobsPath': "Jobs",
            'datasetPath': "Datasets/test-rtt-123",
            'detectorConfig': "PiquantConfig/PIXL/v6",
            'elements': ["Fe","Ca","Ti","K"],
            'parameters': "",
            'pmcListName': nodeName+'.pmcs',
            'command': "map"
        }

        paramJson = json.dumps(paramObj, indent=4, sort_keys=True)
        envVars = self.makeEnvVars(paramJson)

        ran = subprocess.run(['./PiquantRunner'],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=envVars)

        if ran.returncode != 0:
            print(ran.stdout.decode('utf-8'))
            print(ran.stderr.decode('utf-8'))

        self.assertEqual(ran.returncode, 0)

        # Check that the output made it to S3, where we expect, and that it's valid
        s3 = boto3.resource('s3')

        s3.meta.client.download_file(bucket, 'Jobs/'+jobID+'/output/node00012-roi.pmcs_result.csv', make_output_path('downloaded_output_roi.csv'))

        compare_outputs(self, 'downloaded_output_roi.csv', 'whole_runner_roi_map.csv', ['', '', ''])

    def test_run_piquant_s3_fit(self):
# Only here because for some reason they're not bothering fixing this warning...
# https://github.com/boto/boto3/issues/454
        warnings.filterwarnings("ignore", category=ResourceWarning, message="unclosed.*<ssl.SSLSocket.*>")

        bucket = 'test-piquant'
        jobID = 'test-job-125'
        envName = 'piquant-unit-test'
        nodeName = 'node00012-roi'

        prep_cloud(bucket, jobID, envName, nodeName)

        paramObj = {
            'runtimeEnv': envName,
            'datasetsBucket': bucket,
            'jobBucket': bucket,
            'configBucket': bucket,
            'jobId': jobID,
            'jobsPath': "Jobs",
            'datasetPath': "Datasets/test-rtt-123",
            'detectorConfig': "PiquantConfig/PIXL/v6",
            'elements': ["Fe","Ca","Ti","K"],
            'parameters': "",
            'pmcListName': nodeName+'.pmcs',
            'command': "quant"
        }

        paramJson = json.dumps(paramObj, indent=4, sort_keys=True)
        envVars = self.makeEnvVars(paramJson)

        ran = subprocess.run(['./PiquantRunner'],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=envVars)

        if ran.returncode != 0:
            print(ran.stdout.decode('utf-8'))
            print(ran.stderr.decode('utf-8'))

        self.assertEqual(ran.returncode, 0)

        # Check that the output made it to S3, where we expect, and that it's valid
        s3 = boto3.resource('s3')

        s3.meta.client.download_file(bucket, 'Jobs/'+jobID+'/output/node00012-roi.pmcs_result.csv', make_output_path('downloaded_fit_output.csv'))

        #with open(make_output_path('downloaded_fit_output.csv'), "r") as fit_file:
        #    print(fit_file.readline())
        #    print(fit_file.readline())

        #with open('./test-data/expected-output/quant_pmcsfile.csv', "r") as fit_file:
        #    print(fit_file.readline())
        #    print(fit_file.readline())

        compare_outputs(self, 'downloaded_fit_output.csv', 'quant_pmcsfile.csv', ['', '', ''], 1)

# In case we just wanna prep the S3 data & see how it looks, we can run this script on its own...
if __name__ == '__main__':
    WholeRunnerTester().prep_cloud()
