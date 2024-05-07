package main

import (
	"bytes"
	"encoding/json"
	"fmt"
	"io"
	"os"
	"path"
	"path/filepath"
	"strings"
	"testing"

	"github.com/aws/aws-sdk-go/aws"
	"github.com/aws/aws-sdk-go/service/s3"
	"github.com/pixlise/core/v4/api/quantification/quantRunner"
	"github.com/pixlise/core/v4/core/awsutil"
	"github.com/pixlise/core/v4/core/fileaccess"
	"github.com/pixlise/core/v4/core/logger"
)

func Example_getDatasetFileFromPMCList() {
	fmt.Println(getDatasetFileFromPMCList("../test/data/pixlise-datasets/list.pmcs"))

	// Output:
	// 5x11dataset.bin <nil>
}

func Example_prepConfigsForPiquant() {
	params := quantRunner.PiquantParams{
		RunTimeEnv:     "test",
		JobID:          "job-123",
		JobsPath:       "Jobs",
		DatasetPath:    "Downloads/SOL-00001/Experiment-00002",
		DetectorConfig: "PiquantConfig/PIXL",
		Elements: []string{
			"Al",
			"Ti",
			"Ca",
			"Fe",
		},
		Parameters:        "-q,pPIETXCFsr -b,0,12,60,910,280,16 -t,6",
		DatasetsBucket:    "dataset-bucket",
		PiquantJobsBucket: "job-bucket",
		ConfigBucket:      "config-bucket",
		QuantName:         "config test",
		PMCListName:       "file1.pmcs",
	}

	var mockS3 awsutil.MockS3Client
	l := &logger.StdOutLogger{}

	// Listing returns 1 item, get status returns error, check that it still requests 2nd item, 2nd item will fail to parse
	// but the func should still upload a blank jobs.json
	mockS3.ExpListObjectsV2Input = []s3.ListObjectsV2Input{
		{
			Bucket: aws.String(params.ConfigBucket), Prefix: aws.String(params.DetectorConfig),
		},
	}
	mockS3.QueuedListObjectsV2Output = []*s3.ListObjectsV2Output{
		{
			Contents: []*s3.Object{
				{Key: aws.String("PiquantConfig/PIXL/AP04_LVCMFM01_Teflon_1800s_011919_0602_28kV_20uA_10C_Efficiency.txt")},
				{Key: aws.String("PiquantConfig/PIXL/Config_PIXL_FM_ElemCal_CMH_May2019.msa")},
				{Key: aws.String("PiquantConfig/PIXL/FM_3_glasses_ECF_06_19_2019.txt")},
				{Key: aws.String("PiquantConfig/PIXL/FM_Efficiency_profile_Teflon_05_26_2019.txt")},
				{Key: aws.String("PiquantConfig/PIXL/config.json")},
			},
		},
	}

	mockS3.ExpGetObjectInput = []s3.GetObjectInput{
		{
			Bucket: aws.String(params.ConfigBucket), Key: aws.String("PiquantConfig/PIXL/AP04_LVCMFM01_Teflon_1800s_011919_0602_28kV_20uA_10C_Efficiency.txt"),
		},
		{
			Bucket: aws.String(params.ConfigBucket), Key: aws.String("PiquantConfig/PIXL/Config_PIXL_FM_ElemCal_CMH_May2019.msa"),
		},
		{
			Bucket: aws.String(params.ConfigBucket), Key: aws.String("PiquantConfig/PIXL/FM_3_glasses_ECF_06_19_2019.txt"),
		},
		{
			Bucket: aws.String(params.ConfigBucket), Key: aws.String("PiquantConfig/PIXL/FM_Efficiency_profile_Teflon_05_26_2019.txt"),
		},
		{
			Bucket: aws.String(params.ConfigBucket), Key: aws.String("PiquantConfig/PIXL/config.json"),
		},
		{
			Bucket: aws.String(params.PiquantJobsBucket), Key: aws.String("Jobs/job-123/file1.pmcs"),
		},
		{
			Bucket: aws.String(params.DatasetsBucket), Key: aws.String("Downloads/SOL-00001/Experiment-00002/TheDataset.bin"),
		},
	}
	mockS3.QueuedGetObjectOutput = []*s3.GetObjectOutput{
		{
			Body: io.NopCloser(bytes.NewReader([]byte(`Contents of AP04_LVCMFM01_Teflon_1800s_011919_0602_28kV_20uA_10C_Efficiency.txt`))),
		},
		{
			Body: io.NopCloser(bytes.NewReader([]byte(`Contents of Config_PIXL_FM_ElemCal_CMH_May2019.msa`))),
		},
		{
			Body: io.NopCloser(bytes.NewReader([]byte(`Contents of FM_3_glasses_ECF_06_19_2019.txt`))),
		},
		{
			Body: io.NopCloser(bytes.NewReader([]byte(`Contents of FM_Efficiency_profile_Teflon_05_26_2019.txt`))),
		},
		{
			Body: io.NopCloser(bytes.NewReader([]byte(`{
	"description": "PIXL configuration",
	"config-file": "Config_PIXL_FM_ElemCal_CMH_May2019.msa",
	"optic-efficiency": "AP04_LVCMFM01_Teflon_1800s_011919_0602_28kV_20uA_10C_Efficiency.txt",
	"calibration-file": "FM_3_glasses_ECF_06_19_2019.txt",
	"standards-file": ""
}`))),
		},
		{
			Body: io.NopCloser(bytes.NewReader([]byte(`TheDataset.bin
1
2
3
4`))),
		},
		{
			Body: io.NopCloser(bytes.NewReader([]byte(`Contents of TheDataset.bin`))),
		},
	}

	fs := fileaccess.MakeS3Access(&mockS3)
	paths, err := prepConfigsForPiquant(params, l, &fs)
	fmt.Println(err)

	fmt.Println(filepath.Base(paths.config))
	fmt.Println(filepath.Base(paths.calibration))
	fmt.Println(filepath.Base(paths.dataset))
	fmt.Println(filepath.Base(paths.pmcList))

	// Check the optic file setup worked
	data, err := os.ReadFile(paths.config)
	fmt.Println(err)
	datastr := string(data)
	fmt.Printf("%v\n", strings.Contains(datastr, fmt.Sprintf("\n##OPTICFILE : %v", filepath.Dir(paths.config))))

	fmt.Println(err)
	fmt.Println(mockS3.FinishTest())

	// Output:
	// <nil>
	// Config_PIXL_FM_ElemCal_CMH_May2019.msa
	// FM_3_glasses_ECF_06_19_2019.txt
	// TheDataset.bin
	// file1.pmcs
	// <nil>
	// true
	// <nil>
	// <nil>
}

func Example_saveOutputs() {
	const jobBucket = "job-bucket"
	const jobPath = "jobs/some-job-123"

	var mockS3 awsutil.MockS3Client
	l := &logger.StdOutLogger{}

	// Set up expected S3 calls
	mockS3.ExpPutObjectInput = []s3.PutObjectInput{
		{
			Bucket: aws.String(jobBucket), Key: aws.String(path.Join(jobPath, "output", "files_001.pmcs_result.csv")), Body: bytes.NewReader([]byte(`csv file contents`)),
		},
		{
			Bucket: aws.String(jobBucket), Key: aws.String(path.Join(jobPath, "piquant-logs", "files_001.pmcs_piquant.log")), Body: bytes.NewReader([]byte(`log file contents`)),
		},
		{
			Bucket: aws.String(jobBucket), Key: aws.String(path.Join(jobPath, "piquant-logs", "files_001.pmcs_stdout.log")), Body: bytes.NewReader([]byte(`piquant stdout`)),
		},
	}

	mockS3.QueuedPutObjectOutput = []*s3.PutObjectOutput{
		{},
		{},
		{},
	}

	// Create some "output" files
	tmpPath := os.TempDir()
	csvFile, err := os.CreateTemp(tmpPath, "output.csv")
	fmt.Printf("%v\n", err)
	_, err = csvFile.WriteString("csv file contents")
	fmt.Printf("%v\n", err)
	csvFile.Close()

	logFile, err := os.CreateTemp(tmpPath, "threads.log")
	fmt.Printf("%v\n", err)
	_, err = logFile.WriteString("log file contents")
	fmt.Printf("%v\n", err)
	logFile.Close()

	fs := fileaccess.MakeS3Access(&mockS3)
	saveOutputs(jobBucket, jobPath, "files_001.pmcs", csvFile.Name(), logFile.Name(), "piquant stdout", l, fs)

	fmt.Println(mockS3.FinishTest())

	// Cleanup
	os.RemoveAll(tmpPath)

	// Output:
	// <nil>
	// <nil>
	// <nil>
	// <nil>
	// <nil>
}

func Test_loadParams(t *testing.T) {
	// Set default QuantParams
	want := "new-pmc"
	var defaultParams quantRunner.PiquantParams
	defaultParams.PMCListName = want
	paramsJson, _ := json.Marshal(defaultParams)
	os.Setenv("QUANT_PARAMS", string(paramsJson))

	// Test that default params are loaded
	t.Run("Default", func(t *testing.T) {
		if got := loadParams(); got.PMCListName != want {
			t.Errorf("loadParams.PMCListName = %v, want %v", got, want)
		}
	})

	// Override PMCListName via JOB_COMPLETION_INDEX env var
	os.Setenv("JOB_COMPLETION_INDEX", "41")
	want = "node00042.pmcs"
	t.Run("NodeIndexOverride", func(t *testing.T) {
		if got := loadParams(); got.PMCListName != want {
			t.Errorf("loadParams.PMCListName = %v, want %v", got, want)
		}
	})

	// Cleanup
	os.Unsetenv("JOB_COMPLETION_INDEX")
	os.Unsetenv("QUANT_PARAMS")
}
