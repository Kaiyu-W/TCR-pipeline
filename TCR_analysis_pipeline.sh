#!/bin/bash

# configure tools
vdjtools=/mnt/c/vdjtools-1.2.1/vdjtools-1.2.1.jar
vdjtools_patch=/mnt/c/vdjtools-1.2.1/vdjtools-patch.sh
stats_py=/mnt/c/my_script/TCR_stats.py

# configure parameters received from shell
Help() {
	echo -e "\nUsage:\n TCR_analysis_pipeline [options: -iotmc/h]\n"
	echo " TCR analysis, with MiXCR/Cellranger result inputting, VDJTools and Custom R/Python analysis."
	echo -e " Here needs grouped Clonotypes files, as well as corresponding Metadata file for sample information.\n"
	echo "Options:"
	printf " %-30sDirectory of input files, such like MixCR / Cellranger(processed) results\n" "-i <data_dir>"
	printf " %-30sDirectory of output of TCR pipeline\n" "-o <res_dir>"
	printf " %-30sVDJ type for analysis, choose from 'ALL IGH IGK IGL TRA TRB TRD TRG'\n" "-t <VDJ_type>"
	printf " %-30sAbsolute Path of metadata file, including group infomation of samples\n" "-m <metadata_file>"
	printf " %-30suse -c when format transformation from input to what vdjtools needs\n" "-c"
	printf " %-30sHelp\n" "-h"
	echo -e "\nPipeline Written by Kaiyu Wang, in SIBCB of CAS, Shanghai, in Sep 2021.\n"
	exit 0
}
n_para=0
Convert=FALSE
while getopts 'hci:o:t:m:' OPT; do
	case ${OPT} in
       	c)
			Convert=TRUE
			echo -e "\nFirst transfer the input file to vdjtools format"
			;;
      	i)
        	data_dir=${OPTARG}
        	n_para=$[$n_para + 1]
        	;;
      	o)
	        res_dir=${OPTARG}
	        n_para=$[$n_para + 1]
	        ;;
      	t)
	        VDJ_type=${OPTARG}
	        n_para=$[$n_para + 1]
	        ;;
      	m)
			metadata=${OPTARG}
			n_para=$[$n_para + 1]
		;;
      	h)
       		Help
       		exit 0
        	;;
      	*)
			echo -e "\nERROR: Unrecognized parameter. Please check it carefully!"
			Help
	        exit 1
	        ;;
    esac
done
if [ $n_para -ne 4 ]; then
	echo -e "\nERROR: Received not enough parameters!"
	Help
	exit 1
fi
if [ ! -d $data_dir ]; then
	echo -e "\nERROR: Data_dir not exist!"
	Help
	exit 1
fi
if [[ "ALL IGH IGK IGL TRA TRB TRD TRG" =~ "$VDJ_type" ]]; then
	Nothing=0
else
	echo -e "\nERROR: VDJ_type error! Only choose from 'ALL IGH IGK IGL TRA TRB TRD TRG'!"
	Help
	exit 1
fi
if [ ! -f $metadata ]; then
	echo -e "\nERROR: Metadata not exist!"
	Help
	exit 1
fi
if [ $(echo $metadata | grep -c "^/") -ne 1 ]; then
	metadata=$(pwd)/$metadata
fi
if [ ! -d $res_dir ]; then
	mkdir -p $res_dir
fi

#####################
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/home/wky/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/wky/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/home/wky/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/wky/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<
#####################

conda activate vdjtools

# process
# 1. MiXCR process
# not shown

# 2. VDJTools process and plot
# format convert
if [[ $Convert == TRUE ]]; then
	java -jar $vdjtools Convert -S mixcr $data_dir/*.clonotypes.${VDJ_type}.txt $res_dir
else
	cp $data_dir/*.clonotypes.${VDJ_type}.txt $res_dir
fi

if [ ! -f $res_dir/*.clonotypes.${VDJ_type}.txt ]; then
	echo -e "Data didn't exist! \n\tFiles should be as format of {SAMPLE}TCR.clonotypes.{VDJ_type}.txt\n"
	exit 1
fi

# process and plot
cd $res_dir

java -jar $vdjtools CalcBasicStats *.clonotypes.${VDJ_type}.txt .
java -jar $vdjtools CalcPairwiseDistances -p --low-mem *.clonotypes.${VDJ_type}.txt .
java -jar $vdjtools JoinSamples -p *.clonotypes.${VDJ_type}.txt . # only the first 5 sample at most
java -jar $vdjtools RarefactionPlot *.clonotypes.${VDJ_type}.txt .

for Sample in $(gawk -F"\t" '{if(NR!=1)print $1}' $metadata); do
	java -jar $vdjtools PlotQuantileStats -t 5 ${Sample}.clonotypes.${VDJ_type}.txt ./$Sample
	java -jar $vdjtools PlotFancySpectratype -t 20 ${Sample}.clonotypes.${VDJ_type}.txt ./$Sample
	java -jar $vdjtools PlotSpectratypeV -u ${Sample}.clonotypes.${VDJ_type}.txt ./$Sample
	java -jar $vdjtools PlotSpectratypeV ${Sample}.clonotypes.${VDJ_type}.txt ./$Sample
	$vdjtools_patch java -jar $vdjtools PlotFancyVJUsage ${Sample}.clonotypes.${VDJ_type}.txt ./$Sample
done

java -jar $vdjtools ClusterSamples -p -l . .
java -jar $vdjtools CalcSegmentUsage -p *.clonotypes.${VDJ_type}.txt .

# 3. Python and R-DESeq2 analysis / plot

python $stats_py -i ./ -o ./res_python/ -m $metadata -t $VDJ_type

echo ${VDJ_type} OVER!
