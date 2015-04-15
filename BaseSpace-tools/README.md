BaseSpace-tools
===============

Install "basespace-python-sdk" first from this link:  
https://developer.basespace.illumina.com/docs/content/documentation/sdk-samples/python-sdk-overview

Get your accessToken by folowing the steps 1-5 from this link:
https://support.basespace.illumina.com/knowledgebase/articles/403618-python-run-downloader

#### DownloadData.py
download sequencing data from BaseSpace account, need your aceesToken and project name (the client\_key and client\_secret are not necessary)

```
usage: DownloadData.py [-h] [-k KEY] [-s SECRET] [-t TOKEN] -p PROJECT
                       [-d DIRECTORY]

optional arguments:
  -h, --help            show this help message and exit
  -k KEY, --key KEY     the client_key, default:
                        4f7366779990451799fc491d4f6f51b5
  -s SECRET, --secret SECRET
                        the client_serect, default:
                        c88e8c58e5814d04b3c662dc199693ab
  -t TOKEN, --token TOKEN
                        the accessToken for the app [NECESSARY], default:
                        639f9af2f031415cb89bfeef18716b72
  -p PROJECT, --project PROJECT
                        query project full name.
  -d DIRECTORY, --directory DIRECTORY
                        local download file folder, default: current folder
``` 
