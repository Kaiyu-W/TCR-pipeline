import os
import os.path
import sys 
sys.path.append('/tmp/')
import getopt

try:
    opts, args = getopt.getopt(sys.argv[1:], "hi:o:m:t:", ["help","InPath=","OutPath=", "MetaFile=", "VDJType="])
except getopt.GetoptError:
    print('\nOptions:\n\tpython TCR_stats.py -i <inputpath/> \n\t\t\t    -o <outputpath/> \n\t\t\t    -m <metafile> \n\t\t\t    -t <VDJType> \n\t\t\t    -h --help\nDefault sample names are "Cryo" and "NonCryo", please make sure metafile corresponding to that or change the script to receive more parameters.\n')
    sys.exit(2)

n_arg=0
for opt, arg in opts:
    if opt in ("-h", "--help"):
        print('\nOptions:\n\tpython TCR_stats.py -i <inputpath/> \n\t\t\t    -o <outputpath/> \n\t\t\t    -m <metafile> \n\t\t\t    -t <VDJType> \n\t\t\t    -h --help\n')
        sys.exit()
    elif opt in ("-i", "--InPath"):
        input_path = arg
        n_arg = n_arg + 1
    elif opt in ("-o", "--OutPath"):
        output_path = arg
        n_arg = n_arg + 1
    elif opt in ("-m", "--MetaFile"):
        sample_info = arg
        n_arg = n_arg + 1
    elif opt in ("-t", "--VDJType"):
        VDJ_type = arg
        n_arg = n_arg + 1
    else:
        assert Fasle, "Unhandled Option"

if n_arg != 4:
    print('\nAll options needed to input!\nOptions:\n\tpython TCR_stats.py -i <inputpath/> \n\t\t\t    -o <outputpath/> \n\t\t\t    -m <metafile> \n\t\t\t    -t <VDJType> \n\t\t\t    -h --help\nDefault sample names are "Cryo" and "NonCryo", please make sure metafile corresponding to that or change the script to receive more parameters.\n')
    sys.exit(2)

if not os.path.exists(input_path):
    print('\nUnknown input-path!\n')
    sys.exit(2)
if not os.path.exists(sample_info):
    print('\nUnknown metafile(sample-info)!\n')
    sys.exit(2)
if VDJ_type not in ['ALL','IGH','IGK','IGL','TRA','TRB','TRD','TRG']:
    print('\nUnknown VDJ-Type! Only choose from "ALL IGH IGK IGL TRA TRB TRD TRG"!\n')
    sys.exit(2)
if not os.path.exists(output_path):
    os.makedirs(output_path)

import matplotlib as mpl
mpl.use('Agg')

import TCR_utilities as myTCR

import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = (10.0, 10.0)

# from IPython import get_ipython
# use 'os.system' to run cmd
diff_analysis_r='/tmp/TCR_diff_analysis.r'

myTCR.Richness_BarPlot(
    basicstats_file = input_path + 'basicstats.txt',
    meta_file = sample_info,
    output_path = output_path,
    show = False,
    save = True,
    meta_Sample = 'Sample', 
    meta_Group = 'Group', 
    Group1 = 'Cryo', 
    Group2 = 'NonCryo'
    )

myTCR.Length_BarPlot(
    input_path = input_path,
    meta_file = sample_info,
    output_path = output_path,
    VDJ_type = VDJ_type,
    show = False,
    save = True,
    meta_Sample = 'Sample', 
    meta_Group = 'Group', 
    Group1 = 'Cryo', 
    Group2 = 'NonCryo'
    )

myTCR.CountProportion_LappedBarPlot(
    input_path = input_path,
    meta_file = sample_info,
    output_path = output_path,
    VDJ_type = VDJ_type, 
    show = False, 
    save = True, 
    meta_Sample = 'Sample', 
    meta_Group = 'Group', 
    Group1 = 'Cryo', 
    Group2 = 'NonCryo'
    )

myTCR.Abundance_Proportion_relatplot(
    input_path = input_path,
    meta_file = sample_info,
    output_path = output_path,
    VDJ_type = VDJ_type, 
    show = False, 
    save = True, 
    meta_Sample = 'Sample', 
    meta_Group = 'Group', 
    Group1 = 'Cryo', 
    Group2 = 'NonCryo'
    )

myTCR.overlapped_CDR3_process(
    intersect_batch_aa = input_path + 'intersect.batch.aa.txt' , 
    output_overlapped_CDR3 = output_path + 'overlapped_CDR3.txt'
    )

myTCR.overlapped_CDR3_clustermap(
    input_overlapped_CDR3 = output_path + 'overlapped_CDR3.txt', 
    output_path = output_path, 
    show = False, 
    save = True
    )

myTCR.diversity_evenness_clonality(
    input_path = input_path,
    meta_file = sample_info,
    output_path = output_path,
    VDJ_type = VDJ_type, 
    show = False, 
    save = True, 
    meta_Sample = 'Sample', 
    meta_Group = 'Group', 
    Group1 = 'Cryo', 
    Group2 = 'NonCryo'
    )

if os.path.exists(diff_analysis_r):
    for VDJ_gene in ['V', 'D', 'J']:
        myTCR.Build_VDJ_ExpMtx(
            input_path = input_path,
            VDJ_gene = VDJ_gene, 
            meta_file = sample_info,
            output_path = output_path,
            VDJ_type = VDJ_type, 
            meta_Sample = 'Sample'
            )
        VDJ_gene_x = " " + VDJ_gene + " "
        os.system('Rscript ' + diff_analysis_r + VDJ_gene_x + output_path + ' ' + output_path + ' ' + sample_info)
        myTCR.Diffexp_VolcanoPlot(
            diffexp_result_VDJ = output_path + 'diffexp_result_' + VDJ_gene + '.txt', 
            VDJ_gene = VDJ_gene, 
            output_path = output_path, 
            show = False, 
            save = True
            )
        myTCR.Diffexp_Heatmap(
            diffexp_result_VDJ = output_path + 'diffexp_result_' + VDJ_gene + '.txt', 
            VDJ_gene = VDJ_gene, 
            meta_file = sample_info,
            output_path = output_path,
            padj_thre = 0.05, 
            logFC_thre = 0, 
            show = False, 
            save = True, 
            meta_Sample = 'Sample',
            meta_Group = 'Group', 
            Group1 = 'Cryo', 
            Group2 = 'NonCryo'
            )
