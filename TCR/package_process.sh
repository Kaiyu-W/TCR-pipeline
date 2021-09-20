#!/bin/bash
wd=$(which $0 | sed "s/$(basename $0)//")
cd $wd

tar -czf TCR_tmp.tar.gz TCR_diff_analysis.r TCR_stats_tmp.py vdjtools-1.2.1.jar vdjtools-patch.sh TCR_utilities.py
cp TCR_analysis_pipeline_model TCR_analysis_pipeline
cat TCR_tmp.tar.gz >> TCR_analysis_pipeline

# Now get the final package 'TCR_analysis_pipeline' !

# test
./TCR_analysis_pipeline -h