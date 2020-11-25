#!/bin/bash
echo "Starting job"
echo "copying proxy file to /tmp area"
cp x509up_u123627 /tmp/x509up_u123627
echo "copy done..."
echo "start running"
python codeFile/GetDatavariables.py $*
echo "running done"
