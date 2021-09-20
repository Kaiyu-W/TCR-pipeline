# `TCR-pipeline`
This pipeline is designed to process `TCR-seq` data, with input data generating from `MiXCR` or `CellRanger` (the latter of which file format needs to transform by hand).

It needs `Java`, `Python`, `R`, each with extra packages `VDJTools`, `Numpy`/`Pandas`/`Matplotlib`/`Seaborn`, `DESeq2`/`ggplot2`.   
Please make sure Java/Python/Rscript works, and install the packages in advance, otherwise the pipeline will failed to work. Recommended to create `Conda` environment to manager and call the packages.

The way this script works is to unzip the temporary script into __`/tmp/`__, and delete that at the end of the program.  
The `TCR_analysis_pipeline (version 0)` generated from package_process.sh in `TCR0/` includes `vdjtools`, which means it'll have no need for configuration of vdjtools.  
However it needs to install `Rpackages/` into `/tmp/` __by hand__, or the plot function will not work successfully.  

To solve this problem, it recommends that the users assign the path of vdjtools.jar by `-v option`, with `Rpackages/` that vdjtools depends on __in the same directory__.   
This version is named as `TCR_analysis_pipeline` in `main`, with processed codes in `TCR/`.   

Usage:
 TCR_analysis_pipeline [options: -iotmc/h]  

 TCR analysis, with MiXCR/Cellranger result inputting, VDJTools and Custom R/Python analysis.  
 Make sure python packages 'Matplotlib/Seaborn/Numpy/Pandas' and R packages 'DESeq2/ggplot2' installed in your environment.  
 Here needs grouped Clonotypes files, as well as corresponding Metadata file for sample information.  

Options:  
 -i <data_dir>                 Directory of input files, such like MixCR / Cellranger(processed) results  
 -o <res_dir>                  Directory of output of TCR pipeline  
 -t <VDJ_type>                 VDJ type for analysis, choose from 'ALL IGH IGK IGL TRA TRB TRD TRG'  
 -m <metadata_file>            Absolute Path of metadata file, including group infomation of samples  
 -v <VDJTools_jar>             Path of VDJTools code file (.jar), and make sure Rpackages/ existed in the same directory  
 -c                            use -c when format transformation from input to what vdjtools needs  
 -h                            Help  

Example:  
 cd /mnt/e/Cryo-TCR/data/TCR_data/TCR_Raw_mixcr  
 TCR_analysis_pipeline -i ./ -o ./temp_res/ -c -t ALL -m ../metadata.txt -v /mnt/c/vdjtools-1.2.1/vdjtools-1.2.1.jar  

More details please see the script. Please contact me if you need help.   
I'll complete the workflow in the future if I have more time as well as skills, since it lacks the code that helps to install the dependent packages.

Pipeline written by `Kaiyu Wang`, in SIBCB of CAS, Shanghai, in Sep 2021.
