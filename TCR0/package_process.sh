#!/bin/bash

if [[ $0 == bash ]]; then
	wd=$(which $1 | sed "s/$(basename $1)//")
else
	wd=$(which $0 | sed "s/$(basename $0)//")
fi

cd $wd

tar -czf TCR_tmp.tar.gz TCR_utilities.py TCR_diff_analysis.r TCR_stats_tmp.py vdjtools-1.2.1.jar vdjtools-patch.sh 
cp TCR_analysis_pipeline_model TCR_analysis_pipeline
cat TCR_tmp.tar.gz >> TCR_analysis_pipeline

# Now get the final package 'TCR_analysis_pipeline' !

# test
./TCR_analysis_pipeline -h
