# `TCR-pipeline`
This pipeline is designed to process `TCR-seq` data, with input data generating from `MiXCR` or `CellRanger` (the latter of which file format needs to transform by hand).
  
It needs `Java`, `Python`, `R`, each with extra packages `VDJTools`, `Numpy`/`Pandas`/`Matplotlib`/`Seaborn`, `DESeq2`/`ggplot2`.   
Please make sure Java/Python/Rscript works, and install the packages in advance, otherwise the pipeline will failed to work. Recommended to create `Conda` environment to manager and call the packages.
  
The way this script works is to unzip the temporary script into __`/tmp/`__, and delete that at the end of the program.  
  
The `TCR_analysis_pipeline (version 0)` generated from package_process.sh in `TCR0/` includes `vdjtools`, which means it'll have no need for configuration of vdjtools (-v option do not work).  
However it needs to install `Rpackages/` into `/tmp/` __by hand__, or the plot function will not work successfully.  
  
To solve this problem, it recommends that the users assign the path of vdjtools.jar by `-v option`, with `Rpackages/` that vdjtools depends on __in the same directory__.   
This version is named as `TCR_analysis_pipeline` in `main`, with preprocessed codes in `TCR/`.   
  
  
___Usage___:  
```Shell
 TCR_analysis_pipeline [options: -iotmvc/h]  
```
 
 TCR analysis, with MiXCR/Cellranger result inputting, VDJTools and Custom R/Python analysis.  
 Make sure python packages 'Matplotlib/Seaborn/Numpy/Pandas' and R packages 'DESeq2/ggplot2' installed in your environment.  
 Here needs grouped Clonotypes files, as well as corresponding Metadata file for sample information.  
  
___Options___:  
| option | format | detail |  
| :----: | :----: | ------ |  
| -i | <data_dir> | Directory of input files, such like MixCR / Cellranger (processed) results |  
| -o | <res_dir> | Directory of output of TCR pipeline |  
| -t | <VDJ_type> | VDJ type for analysis, choose from 'ALL IGH IGK IGL TRA TRB TRD TRG' |  
| -m | <metadata_file> | Absolute Path of metadata file, including group infomation of samples |  
| -v | <VDJTools_jar> | Path of VDJTools code file (.jar), and make sure directory 'Rpackages/' existed in the same directory |  
| -c |  | use -c when format transformation from input to what vdjtools needs |  
| -h |  | Help |  
  
___Example___:  
```Shell
# to prepare the input files from cellranger and mixcr output
<!-- 
# cellranger vdj process fastq into filter.fastq
cd xxx
for i in *_S1_L001_I1_001.fastq.gz; do
	name=$( echo $i | sed "s/\(.*\)\(_S1_L001_I1_001.fastq.gz\)/\1/" )
	cellranger vdj --id=$name \
		--reference=refdata-cellranger-vdj-GRCm38-alts-ensembl-5.0.0 \
		--fastqs=./ \
		--sample=$name &
done
# mixcr process all filter.fastq
for i in 1 2 3 4; do
  for j in CryoTCR NonCryoTCR; do
        analysis_name=${j}_${i}
        input_file=${j}_${i}.filter.fastq
        mixcr analyze amplicon \
            -s mmu \
            --starting-material rna \
            --5-end no-v-primers --3-end c-primers \
            --adapters adapters-present \
            $input_file ./res/$analysis_name
        echo $analysis_name OVER!
  done
done
 -->
 
$cd /mnt/e/Cryo-TCR/data/TCR_data/TCR_Raw_mixcr  
$ls -l *txt  
   
#-rwxrwxrwx 1 wky wky 5460960 Sep  8 20:54 CryoTCR_1.clonotypes.ALL.txt  
#-rwxrwxrwx 1 wky wky 5491474 Sep  8 20:54 CryoTCR_2.clonotypes.ALL.txt  
#-rwxrwxrwx 1 wky wky 5321901 Sep  8 20:54 CryoTCR_3.clonotypes.ALL.txt  
#-rwxrwxrwx 1 wky wky 5462734 Sep  8 20:54 CryoTCR_4.clonotypes.ALL.txt  
#-rwxrwxrwx 1 wky wky 4713741 Sep  8 20:54 NonCryoTCR_1.clonotypes.ALL.txt  
#-rwxrwxrwx 1 wky wky 4776696 Sep  8 20:54 NonCryoTCR_2.clonotypes.ALL.txt  
#-rwxrwxrwx 1 wky wky 4946934 Sep  8 20:54 NonCryoTCR_3.clonotypes.ALL.txt  
#-rwxrwxrwx 1 wky wky 4939794 Sep  8 20:54 NonCryoTCR_4.clonotypes.ALL.txt  
  
$cat ../metadata.txt # tab delimiter  
#Sample  Group  
#CryoTCR_1       Cryo  
#CryoTCR_2       Cryo  
#CryoTCR_3       Cryo  
#CryoTCR_4       Cryo  
#NonCryoTCR_1    NonCryo  
#NonCryoTCR_2    NonCryo  
#NonCryoTCR_3    NonCryo  
#NonCryoTCR_4    NonCryo  
  
$TCR_analysis_pipeline -i ./                                       \
                       -o ./temp_res/                              \
                       -m ../metadata.txt                          \
                       -v /mnt/c/vdjtools-1.2.1/vdjtools-1.2.1.jar \
                       -t ALL                                      \
                       -c  
```
  
More details please see the script. Please contact me if you need help.   
I'll complete the workflow in the future if I have more time as well as skills, since it lacks the code that helps to install the dependent packages.
  
Pipeline written by __`Kaiyu Wang`__, in SIBCB of CAS, Shanghai, in Sep 2021.  
