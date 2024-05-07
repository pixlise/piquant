// Licensed to NASA JPL under one or more contributor
// license agreements. See the NOTICE file distributed with
// this work for additional information regarding copyright
// ownership. NASA JPL licenses this file to you under
// the Apache License, Version 2.0 (the "License"); you may
// not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing,
// software distributed under the License is distributed on an
// "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
// KIND, either express or implied.  See the License for the
// specific language governing permissions and limitations
// under the License.

package main

import (
	"fmt"
	"log"
	"sync"
	"time"

	"github.com/aws/aws-sdk-go/aws"
	"github.com/aws/aws-sdk-go/aws/awserr"
	"github.com/aws/aws-sdk-go/aws/session"
	"github.com/aws/aws-sdk-go/service/cloudwatchlogs"
	"github.com/pixlise/core/v4/core/logger"
)

// This is heavily based on: github.com/jcxplorer/cwlogger and github.com/mathisve/golang-cloudwatch-logs-example

// CloudwatchLogger - Structure holding API logger internals
type CloudwatchLogger struct {
	cwClient       *cloudwatchlogs.CloudWatchLogs
	logGroupName   string
	logStreamName  string
	logLevel       logger.LogLevel
	sequenceToken  string
	queue          []string
	lock           sync.Mutex
	retentionDays  int
	logIntervalSec time.Duration
	running        bool
}

// InitCloudWatchLogger - initialises the logger, given settings and AWS session
func initCloudWatchLogger(sess *session.Session, logGroupName string, logStreamName string, logLevel logger.LogLevel, retentionDays int, logIntervalSec time.Duration) (logger.ILogger, error) {
	result := CloudwatchLogger{
		cwClient:       cloudwatchlogs.New(sess),
		logGroupName:   logGroupName,
		logStreamName:  logStreamName,
		logLevel:       logLevel,
		sequenceToken:  "",
		queue:          []string{},
		lock:           sync.Mutex{},
		retentionDays:  retentionDays,
		logIntervalSec: logIntervalSec,
		running:        true,
	}
	/*
		tok := ""
		totalDel := 0
		for i := 0; i < 100; i++ {
			var err error
			var delcount int
			tok, delcount, err = result.deleteOldGroups("/aws/lambda/", tok)
			if err != nil {
				return &result, err
			}
			totalDel += delcount
		}
	*/
	// Make sure log group exists
	err := result.ensureLogGroupExists(logGroupName, int64(retentionDays))
	if err != nil {
		return &result, err
	}

	// Now open the log stream
	go result.processQueue(logIntervalSec)

	return &result, nil
}

func (l *CloudwatchLogger) Close() {
	// Sleep this thread for long enough that the other thread pumps any messages left in its queue to cloudwatch
	time.Sleep(time.Second * l.logIntervalSec * 2)
}

/*
func (l *CloudwatchLogger) deleteOldGroups(prefix string, contToken string) (string, int, error) {

		in := cloudwatchlogs.DescribeLogGroupsInput{
			Limit: aws.Int64(50),
		}
		if len(contToken) > 0 {
			in.NextToken = aws.String(contToken)
		}

		resp, err := l.cwClient.DescribeLogGroups(&in)
		if err != nil {
			return "", 0, err
		}

		delcount := 0
		for _, logGroup := range resp.LogGroups {
			if strings.HasPrefix(*logGroup.LogGroupName, prefix) {
				_, err := l.cwClient.DeleteLogGroup(&cloudwatchlogs.DeleteLogGroupInput{
					LogGroupName: logGroup.LogGroupName,
				})
				if err != nil {
					fmt.Printf("DeleteLogGroup %v failed: %v\n", *logGroup.LogGroupName, err)
				} else {
					delcount++
				}
			}
		}

		return *resp.NextToken, delcount, nil
	}
*/
func (l *CloudwatchLogger) ensureLogGroupExists(name string, retentionDays int64) error {
	_, err := l.cwClient.CreateLogGroup(&cloudwatchlogs.CreateLogGroupInput{
		LogGroupName: aws.String(name),
	})

	if err != nil {
		// If it already exists, don't fail!
		if aerr, ok := err.(awserr.Error); ok {
			if aerr.Message() == "The specified log group already exists" {
				return nil
			}
		}
		return err
	}

	_, err = l.cwClient.PutRetentionPolicy(&cloudwatchlogs.PutRetentionPolicyInput{
		RetentionInDays: aws.Int64(retentionDays),
		LogGroupName:    aws.String(name),
	})

	return err
}

func (l *CloudwatchLogger) createLogStream(name string) error {
	// Ensure it exists here because it may be deleted at runtime from cloudwatch, if we're creating or re-creating
	// our log stream, it's good to know that the group is there
	err := l.ensureLogGroupExists(l.logGroupName, int64(l.retentionDays))
	if err != nil {
		return err
	}

	_, err = l.cwClient.CreateLogStream(&cloudwatchlogs.CreateLogStreamInput{
		LogGroupName:  aws.String(l.logGroupName),
		LogStreamName: aws.String(name),
	})

	return err
}

// processQueue will process the log queue
func (l *CloudwatchLogger) processQueue(logIntervalSec time.Duration) error {
	var logQueue []*cloudwatchlogs.InputLogEvent

	for l.running {
		l.lock.Lock()
		if len(l.queue) > 0 {
			for _, item := range l.queue {
				logQueue = append(logQueue, &cloudwatchlogs.InputLogEvent{
					Message:   aws.String(item),
					Timestamp: aws.Int64(time.Now().UnixNano() / int64(time.Millisecond)),
				})
			}

			l.queue = []string{}
		}

		l.lock.Unlock()

		if len(logQueue) > 0 {
			input := cloudwatchlogs.PutLogEventsInput{
				LogEvents:    logQueue,
				LogGroupName: aws.String(l.logGroupName),
			}

			if l.sequenceToken == "" {
				err := l.createLogStream(l.logStreamName)
				if err != nil {
					// Write to stderr
					log.Printf("createLogStream failed: %v", err)
				}
			} else {
				input = *input.SetSequenceToken(l.sequenceToken)
			}

			input = *input.SetLogStreamName(l.logStreamName)

			resp, err := l.cwClient.PutLogEvents(&input)
			if err != nil {
				// Write to stderr
				log.Printf("PutLogEvents failed: %v", err)
			}

			if resp != nil {
				if resp.NextSequenceToken != nil {
					l.sequenceToken = *resp.NextSequenceToken
				} else {
					l.sequenceToken = ""
				}
			}

			logQueue = []*cloudwatchlogs.InputLogEvent{}
		}

		time.Sleep(time.Second * logIntervalSec)
	}

	// Might be useful seeing this on shutdown...
	fmt.Println("logger processQueue complete")
	return nil
}

// Duplicated from logger because it wasn't public:
var logLevelPrefix = map[logger.LogLevel]string{
	logger.LogDebug: "DEBUG",
	logger.LogInfo:  "INFO",
	logger.LogError: "ERROR",
}

// Log enqueues a log message to be written to a log stream.
//
// The log message must be less than 1,048,550 bytes, and the time must not be
// more than 2 hours in the future, 14 days in the past, or older than the
// retention period of the log group.
//
// This method is safe for concurrent access by multiple goroutines.
func (l *CloudwatchLogger) Printf(level logger.LogLevel, format string, a ...interface{}) {
	// If we're not on this log level, skip
	if l.logLevel > level {
		return
	}

	txt := logLevelPrefix[level] + ": " + fmt.Sprintf(format, a...)

	defer l.lock.Unlock()
	l.lock.Lock()

	// Add to the log queue
	l.queue = append(l.queue, txt)

	// Also write to local stdout
	log.Println(txt)
}

// Debugf - Print debug to log, with format string
func (l *CloudwatchLogger) Debugf(format string, a ...interface{}) {
	l.Printf(logger.LogDebug, format, a...)
}

// Infof - Print info to log, with format string
func (l *CloudwatchLogger) Infof(format string, a ...interface{}) {
	l.Printf(logger.LogInfo, format, a...)
}

// Errorf - Print error to log, with format string
func (l *CloudwatchLogger) Errorf(format string, a ...interface{}) {
	l.Printf(logger.LogError, format, a...)
}

func (l *CloudwatchLogger) SetLogLevel(level logger.LogLevel) {
	l.logLevel = level
}
func (l *CloudwatchLogger) GetLogLevel() logger.LogLevel {
	return l.logLevel
}
