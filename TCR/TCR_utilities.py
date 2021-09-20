import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# plt.rcParams['figure.figsize'] = (10.0, 10.0)

def Richness_BarPlot(basicstats_file, meta_file, output_path, show = True, save = True, meta_Sample = 'Sample', meta_Group = 'Group', Group1 = 'Cryo', Group2 = 'NonCryo'):
    list_pos = []
    for line in open(meta_file, 'r'):
        info = line[:-1].split('\t')
        if info[0] != meta_Sample:
            if str(info[1]) == Group1:
                list_pos.append(info[0])

    list_class = []
    list_sample = []
    list_richness = []
    for line in open(basicstats_file, 'r'):
        info=line[:-1].split('\t')
        if info[0] != 'sample_id':
            if info[0].split('.')[0] in list_pos:
                list_class.append(Group1)
            else:
                list_class.append(Group2)
            list_sample.append(info[0].split('.')[0])
            list_richness.append(float(info[3]))

    ax1 = sns.barplot(x = list_class, y = list_richness, capsize = .2, palette = (sns.xkcd_rgb["dark red"], sns.xkcd_rgb["marine blue"]))
    ax1_1 = sns.stripplot(x = list_class, y = list_richness, size = 5, palette = (sns.xkcd_rgb["dark red"], sns.xkcd_rgb["marine blue"]))
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    if save:
        plt.savefig(output_path + 'Richness_barplot_' + meta_Group + '.pdf')
    if show:
        plt.show()
    plt.close()

    ax2 = sns.barplot(x = list_sample, y = list_richness)
    ax2.set_xticklabels(ax2.get_xticklabels(), rotation = 45)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    if save:
        plt.savefig(output_path + 'Richness_barplot_' + meta_Sample + '.pdf')
    if show:
        plt.show()
    plt.close()


def Length_BarPlot(input_path, meta_file, output_path, VDJ_type = 'ALL', show = True, save = True, meta_Sample = 'Sample', meta_Group = 'Group', Group1 = 'Cryo', Group2 = 'NonCryo'):
    dict_sample_LengthCount = {}
    list_pos = []
    for line in open(meta_file, 'r'):
        info = line[:-1].split('\t')
        id = info[0]
        if id != meta_Sample:
            if str(info[1]) == Group1:
                list_pos.append(id)
            read_file = open(input_path + id + '.clonotypes.' + VDJ_type + '.txt', 'r')
            dict_length_count = {}
            for line in read_file:
                if line[0] != 'c':
                    length = len(line[:-1].split('\t')[3])
                    if length < 30:
                        if length not in dict_length_count.keys():
                            dict_length_count[length] = 1
                        else:
                            dict_length_count[length] += 1
            dict_sample_LengthCount[id] = dict_length_count
            read_file.close()

    list_length=[]
    list_class=[]
    list_sample=[]
    list_count=[]
    for sample in dict_sample_LengthCount.keys():
        dict_length_count=dict_sample_LengthCount[sample]
        for length in dict_length_count.keys():
            list_length.append(length)
            list_sample.append(sample)
            if sample in list_pos:
                list_class.append(Group1)
            else:
                list_class.append(Group2)
            list_count.append(dict_length_count[length])

    ax1 = sns.barplot(x = list_length, y = list_count, hue = list_class, errwidth = 0.5, capsize = 0.5, palette = (sns.xkcd_rgb["dark red"], sns.xkcd_rgb["marine blue"]))
    ax1_1 = sns.stripplot(x = list_length, y = list_count, hue = list_class, size = 2, dodge = True, palette = (sns.xkcd_rgb["dark red"], sns.xkcd_rgb["marine blue"]))
    ax1.set_xticklabels(ax1.get_xticklabels(), rotation = 0, fontsize = 8)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    if save:
        plt.savefig(output_path + 'Length_barplot_' + meta_Group + '.pdf')
    if show:
        plt.show()
    plt.close()

    ax2 = sns.barplot(x = list_length, y = list_count, hue = list_sample, errwidth = 0.5, capsize = 0.5)
    ax2.set_xticklabels(ax2.get_xticklabels(), rotation = 0, fontsize = 8)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    if save:
        plt.savefig(output_path+'Length_barplot_' + meta_Sample + '.pdf')
    if show:
        plt.show()
    plt.close()

def CountProportion_LappedBarPlot(input_path, meta_file, output_path, VDJ_type = 'ALL', show = True, save = True, meta_Sample = 'Sample', meta_Group = 'Group', Group1 = 'Cryo', Group2 = 'NonCryo'):
    dict_sample_ReadsCounts = {}
    list_pos = []
    for line in open(meta_file, 'r'):
        info = line[:-1].split('\t')
        id = info[0]
        if id != meta_Sample:
            if str(info[1]) == Group1:
                list_pos.append(id)
            read_file = open(input_path + id + '.clonotypes.' + VDJ_type + '.txt', 'r')
            dict_read_count = {'1':0,'2-3':0,'4-10':0,'11-30':0,'31-100':0,'101-MAX':0}
            number_clonotypes = 0
            for line in read_file:
                if line[0] != 'c':
                    number_clonotypes += 1
                    info = line[:-1].split('\t')
                    count = info[0]
                    if int(count) == 1:
                        dict_read_count['1'] += 1
                    elif int(count) in [2,3]:
                        dict_read_count['2-3'] += 1
                    elif 4 <= int(count) <= 10:
                        dict_read_count['4-10'] += 1
                    elif 11 <= int(count) <= 30:
                        dict_read_count['11-30'] += 1
                    elif 31 <= int(count) <= 100:
                        dict_read_count['31-100'] += 1
                    elif int(count) >= 101:
                        dict_read_count['101-MAX'] += 1
            for key in dict_read_count:
                dict_read_count[key] = dict_read_count[key] / number_clonotypes * 100
            dict_sample_ReadsCounts[id] = dict_read_count
            read_file.close()

    list_sample = []
    list_class = []
    list_read = []
    list_count = []
    dict_read_SampleCount = {'1': {},'2-3': {},'4-10': {},'11-30': {},'31-100': {},'101-MAX': {}}
    print(dict_sample_ReadsCounts)

    for sample in dict_sample_ReadsCounts.keys():
        dict_read_count = dict_sample_ReadsCounts[sample]
        for read in dict_read_count.keys():
            list_sample.append(sample)
            if sample in list_pos:
                list_class.append(Group1)
            else:
                list_class.append(Group2)
            list_read.append(read)
            count = float(dict_read_count[read])
            list_count.append(count)
            dict_sample_count = dict_read_SampleCount[read]
            dict_sample_count[sample] = count
            dict_read_SampleCount[read] = dict_sample_count

    ax1 = sns.barplot(x=list_read,y=list_count,hue=list_class,errwidth=1,capsize=0.2,palette=(sns.xkcd_rgb["dark red"],sns.xkcd_rgb["marine blue"]))
    ax1_1 = sns.stripplot(x=list_read,y=list_count,hue=list_class,dodge=True,palette=(sns.xkcd_rgb["dark red"],sns.xkcd_rgb["marine blue"]))
    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=0,fontsize=10)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    if save:
        plt.savefig(output_path+'CountProportion_lappedplot_' + meta_Group + '.pdf')
    if show:
        plt.show()
    plt.close()

    fig,ax = plt.subplots()
    list_1 = []
    list_2_3 = []
    list_4_10 = []
    list_11_30 = []
    list_31_100 = []
    list_101_MAX = []
    for read in dict_read_SampleCount.keys():
        dict_sample_count = dict_read_SampleCount[read]
        list_sample = []
        for sample in dict_sample_count.keys():
            list_sample.append(sample)
            if read == '1':
                list_1.append(dict_sample_count[sample])
            elif read == '2-3':
                list_2_3.append(dict_sample_count[sample])
            elif read == '4-10':
                list_4_10.append(dict_sample_count[sample])
            elif read == '11-30':
                list_11_30.append(dict_sample_count[sample])
            elif read == '31-100':
                list_31_100.append(dict_sample_count[sample])
            elif read == '101-MAX':
                list_101_MAX.append(dict_sample_count[sample])
    list_1 = np.array(list_1)
    list_2_3 = np.array(list_2_3)
    list_4_10 = np.array(list_4_10)
    list_11_30 = np.array(list_11_30)
    list_31_100 = np.array(list_31_100)
    list_101_MAX = np.array(list_101_MAX)
    width = 0.35
    print(list_1)
    ax.bar(list_sample, list_1,width, label = '1')
    ax.bar(list_sample, list_2_3,width, label = '2_3', bottom = list_1)
    ax.bar(list_sample, list_4_10,width, label = '4-10', bottom = list_1 + list_2_3)
    ax.bar(list_sample, list_11_30,width, label = '11-30', bottom = list_1 + list_2_3 + list_4_10)
    ax.bar(list_sample, list_31_100, width, label = '31-100', bottom = list_1 + list_2_3 + list_4_10 + list_11_30)
    ax.bar(list_sample, list_101_MAX, width, label = '101-MAX', bottom = list_1 + list_2_3 + list_4_10 + list_11_30 + list_31_100)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend()
    ax.set_xticklabels(list_sample, rotation = 45, fontsize = 10)
    if save:
        plt.savefig(output_path+'CountProportion_barplot_' + meta_Sample + '.pdf')
    if show:
        plt.show()
    plt.close()

def Abundance_Proportion_relatplot(input_path, meta_file, output_path, VDJ_type = 'ALL', show = True, save = True, meta_Sample = 'Sample', meta_Group = 'Group', Group1 = 'Cryo', Group2 = 'NonCryo'):
    dict_sample_AbundanceRichness = {}
    list_pos = []
    for line in open(meta_file, 'r'):
        info = line[:-1].split('\t')
        id = info[0]
        if id != meta_Sample:
            if str(info[1]) == Group1:
                list_pos.append(id)
            dict_abundance_richness = {}
            read_file = open(input_path + id + '.clonotypes.' + VDJ_type + '.txt', 'r')
            for line in read_file:
                if line[0] != 'c':
                    info = line[:-1].split('\t')
                    abundance = info[0]
                    if abundance not in dict_abundance_richness.keys():
                        dict_abundance_richness[abundance] = 1
                    else:
                        dict_abundance_richness[abundance] += 1
            dict_sample_AbundanceRichness[id] = dict_abundance_richness
            read_file.close()

    list_abundance = []
    list_richniess = []
    list_sample = []
    list_class = []
    for sample in dict_sample_AbundanceRichness.keys():
        dict_abundance_richness = dict_sample_AbundanceRichness[sample]
        for abundance in dict_abundance_richness.keys():
            list_sample.append(sample)
            if sample in list_pos:
                list_class.append(Group1)
            else:
                list_class.append(Group2)
            list_abundance.append(int(abundance))
            list_richniess.append(int(dict_abundance_richness[abundance]))

    ax = sns.lineplot(x = list_abundance, y = list_richniess, hue = list_sample, style = list_class)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.xscale('log')
    plt.yscale('log')
    if save:
        plt.savefig(output_path+'Abundance_Proportion_relatplot.pdf')
    if show:
        plt.show()
    plt.close()


def overlapped_CDR3_process(intersect_batch_aa, output_overlapped_CDR3):
    # first run vdjtools CalcPairwiseDistances code
    input_file_temp = open(intersect_batch_aa,'r')
    save_file_temp = open(output_overlapped_CDR3,'w')
    list_samples = []
    diction = {}
    row = 1
    for line in input_file_temp:
        if row != 1:
            info = line[:-1].split('\t')
            for sample in [info[0],info[1]]:
                if sample not in list_samples:
                    list_samples.append(sample)
            key = info[0] + 'and' + info[1]
            value = info[4]
            diction[key] = value
        row += 1
    input_file_temp.close()
    word = ''
    for sample in list_samples:
        sample = sample.split('.')[0]
        word += sample + '\t'
    save_file_temp.write(word[:-1] + '\n')
    for sample_row in list_samples:
        word = ''
        for sample_column in list_samples:
            if sample_row + 'and' + sample_column in diction.keys():
                word += str(diction[sample_row + 'and' + sample_column]) + '\t'
            else:
                if sample_column + 'and' + sample_row in diction.keys():
                    word += str(diction[sample_column + 'and' + sample_row]) + '\t'
                else:
                    word += '0' + '\t'
        word = word[:-1] + '\n'
        save_file_temp.write(word)
    save_file_temp.close()

def overlapped_CDR3_clustermap(input_overlapped_CDR3, output_path, show = True, save = True):
    list_samples = []
    list_value_list = []
    number_row = 1
    for line in open(input_overlapped_CDR3):
        if number_row != 1:
            list_value = []
            index = 0
            for value in line[:-1].split('\t'):
                list_value.append(float(value))
                index += 1
            list_value_list.append(list_value)
        else:
            for sample in line[:-1].split('\t'):
                list_samples.append(sample)
        number_row += 1
    data = pd.DataFrame(data = list_value_list, index = list_samples, columns = list_samples)
    ax = sns.clustermap(data, standard_scale = 1, figsize = (10, 10))
    plt.setp(ax.ax_heatmap.xaxis.get_majorticklabels(), rotation = 45)
    if save:
        plt.savefig(output_path + 'overlapped_CDR3.pdf')
    if show:
        plt.show()
    plt.close()

def Shannon_entropy(list_frequency):
    sum = 0
    for frequency in list_frequency:
        sum += frequency * math.log(frequency)
    H = -sum
    return H

def Pielou_evenness(H, richness):
    E = H / math.log(richness)
    return E

def Clonality(E):
    C = 1 - E
    return C

def diversity_evenness_clonality(input_path, meta_file, output_path, VDJ_type = 'ALL', show = True, save = True, meta_Sample = 'Sample', meta_Group = 'Group', Group1 = 'Cryo', Group2 = 'NonCryo'):
    list_pos = []
    list_characteristic = []
    list_value = []
    list_sample = []
    list_class = []
    for line in open(meta_file, 'r'):
        info = line[:-1].split('\t')
        id = info[0]
        if id != meta_Sample:
            if info[1] == Group1:
                list_pos.append(id)
            read_file = open(input_path + id + '.clonotypes.' + VDJ_type + '.txt', 'r')
            richness = 0
            TotalCount = 0
            for line in read_file:
                if line[0] != 'c':
                    count = int(line[:-1].split('\t')[0])
                    richness += 1
                    TotalCount += count
            read_file.close()
            read_file = open(input_path + id + '.clonotypes.' + VDJ_type + '.txt', 'r')
            list_frequency = []
            for line in read_file:
                if line[0] != 'c':
                    count = int(line[:-1].split('\t')[0])
                    list_frequency.append(count / TotalCount)
            diversity = Shannon_entropy(list_frequency)
            evenness = Pielou_evenness(diversity, richness)
            clonality = Clonality(evenness)
            list_characteristic += ['diversity', 'evenness', 'clonality']
            list_value += [diversity, evenness, clonality]
            for time in range(0, 3):
                list_sample.append(id)
                if id in list_pos:
                    list_class.append(Group1)
                else:
                    list_class.append(Group2)
            read_file.close()

    factor_max = max(list_value[0::3])
    for i in range(len(list_value)):
        if i%3 == 0:
            list_value[i] = list_value[i] / factor_max

    print(list_class)
    ax1 = sns.barplot(x = list_characteristic, y = list_value, hue = list_class, errwidth = 0.5, capsize = 0.3, palette = (sns.xkcd_rgb["dark red"], sns.xkcd_rgb["marine blue"]))
    ax1_1 = sns.stripplot(x = list_characteristic, y = list_value, size = 2, hue = list_class, dodge = True, palette = (sns.xkcd_rgb["dark red"], sns.xkcd_rgb["marine blue"]))
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1_2 = ax1.twinx()
    ax1.yaxis.set_ticks_position('right')
    ax1_2.yaxis.set_ticks_position('left')
    ax1.set_ylim(0, 1.05)
    ax1_2.set_ylim(0, factor_max * 1.05)
    if save:
        plt.savefig(output_path + 'diversity_evenness_clonality_barplot_' + meta_Group + '.pdf')
    if show:
        plt.show()
    plt.close()

    ax2 = sns.violinplot(x = list_characteristic, y = list_value, hue = list_class, palette = ('white', 'white'))
    ax2_1 = sns.stripplot(x = list_characteristic, y = list_value, size = 4, hue = list_class, dodge = True, palette = (sns.xkcd_rgb["dark red"], sns.xkcd_rgb["marine blue"]))
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2_2 = ax2.twinx()
    ax2.yaxis.set_ticks_position('right')
    ax2_2.yaxis.set_ticks_position('left')
    ax2.set_ylim(0, 1.05)
    ax2_2.set_ylim(0, factor_max * 1.05)
    if save:
        plt.savefig(output_path + 'diversity_evenness_clonality_violinplot_' + meta_Group + '.pdf')
    if show:
        plt.show()
    plt.close()

    ax3 = sns.barplot(x = list_characteristic, y = list_value, hue = list_sample, errwidth = 0.5, capsize = 0.5)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3_1 = ax3.twinx()
    ax3.yaxis.set_ticks_position('right')
    ax3_1.yaxis.set_ticks_position('left')
    ax3.set_ylim(0, 1.05)
    ax3_1.set_ylim(0, factor_max * 1.05)
    if save:
        plt.savefig(output_path + 'diversity_evenness_clonality_' + meta_Sample + '.pdf')
    if show:
        plt.show()
    plt.close()

def Build_VDJ_ExpMtx(input_path, VDJ_gene, meta_file, output_path, VDJ_type = 'ALL', meta_Sample = 'Sample'):
    matrix_VDJ = output_path + '/matrix_' + VDJ_gene + '.txt'
    GeneIndex_dict = {'V':4, 'D':5, 'J':6}
    GeneIndex = GeneIndex_dict.get(VDJ_gene)

    list_files = []
    for line in open(meta_file, 'r'):
        if line[0] != meta_Sample[0]:
            file=line.split('\t')[0]
            list_files.append(file + '.clonotypes.' + VDJ_type + '.txt')

    list_sample = []
    list_VDJ = []
    list_dict_VDJ = []
    for file in list_files:
        read_file = open(input_path + file,'r')
        sample = file.split('.')[0]
        list_sample.append(sample)
        dict_CDR = {}
        dict_VDJ = {}
        for line in read_file:
            if line[0] != 'c':
                info = line[:-1].split('\t')
                name_VDJ = info[GeneIndex]
                if name_VDJ not in list_VDJ:
                    list_VDJ.append(name_VDJ)
                if name_VDJ not in dict_VDJ.keys():
                    dict_VDJ[name_VDJ] = 1
                else:
                    dict_VDJ[name_VDJ] += 1
        list_dict_VDJ.append(dict_VDJ)
        read_file.close()

    word = 'item'
    for sample in list_sample:
        word += '\t' + sample
    word += '\n'
    VDJ_file = open(matrix_VDJ, 'w')
    VDJ_file.write(word)
    for VDJ in list_VDJ:
        check = []
        info = str(VDJ)
        for dict_VDJ in list_dict_VDJ:
            if VDJ in dict_VDJ.keys():
                info += '\t'+str(dict_VDJ[VDJ])
                check.append(float(dict_VDJ[VDJ]))
            else:
                info += '\t'+'0'
                check.append(0)
        info += '\n'
        if max(check) != 0:
            VDJ_file.write(info)
    VDJ_file.close()
    print(matrix_VDJ + ' expression matrix have been build!')

def Diffexp_VolcanoPlot(diffexp_result_VDJ, VDJ_gene, output_path, show = True, save = True):
    # first run R code with Deseq2 analysis
    output_file = output_path + 'diffexp_VolcanoPlot_' + VDJ_gene + '.pdf'

    list_log2FC = []
    list_padj = []
    index_log2FC = 2
    index_padj = 6
    number_row = 1
    for line in open(diffexp_result_VDJ):
        if number_row == 1:
            index_log2FC = line[:-1].split('\t').index('log2FoldChange')
            index_padj = line[:-1].split('\t').index('padj')
        else:
            if line.split('\t')[index_log2FC] != 'NA' and line.split('\t')[index_padj] != 'NA':
                list_log2FC.append(float(line.split('\t')[index_log2FC]))
                list_padj.append(float(line.split('\t')[index_padj]))
        number_row += 1
    foldchange = np.asarray(list_log2FC)
    pvalue = np.asarray(list_padj)

    result = pd.DataFrame({'pvalue': pvalue, 'FoldChange': foldchange})
    result['log(pvalue)'] = -np.log10(result['pvalue'])
    result['sig'] = 'normal'
    result['size'] = np.abs(result['FoldChange']) / 10

    result.loc[(result.FoldChange > 0) & (result.pvalue < 0.05), 'sig'] = 'up'
    result.loc[(result.FoldChange < 0) & (result.pvalue < 0.05), 'sig'] = 'down'

    ax = sns.scatterplot(x = "FoldChange", y = "log(pvalue)", hue = 'sig', hue_order = ('up', 'down', 'normal'), palette = ("#E41A1C", "#377EB8", "grey"), data = result)
    ax.set_ylabel('-log10Padj', fontweight = 'bold')
    ax.set_xlabel('log2FoldChange', fontweight = 'bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if save:
        plt.savefig(output_file)
    if show:
        plt.show()
    plt.close()

def Diffexp_Heatmap(diffexp_result_VDJ, VDJ_gene, meta_file, output_path, padj_thre = 0.05, logFC_thre = 0, show = True, save = True, meta_Sample = 'Sample', meta_Group = 'Group', Group1 = 'Cryo', Group2 = 'NonCryo'):
    output_file = output_path + 'diffexp_heatmap_' + VDJ_gene + '.pdf'

    list_pos = []
    for line in open(meta_file,'r'):
        info = line[:-1].split('\t')
        if info[0] != meta_Sample:
            if str(info[1]) == Group1:
                list_pos.append(info[0])
                
    list_ids = []
    list_samples = []
    list_value_list = []
    list_sample_colors = []
    list_gene_colors = []
    begin_samples = 7
    index_log2FC = 2
    index_p = 6
    number_row = 1
    for line in open(diffexp_result_VDJ):
        if number_row != 1:
            list_ids.append(line[:-1].split('\t')[0])
            list_value = []
            for value in line[:-1].split('\t')[begin_samples:]:
                list_value.append(float(value))
            list_value_list.append(list_value)
            log2FC = float(line[:-1].split('\t')[index_log2FC])
            p_temp = line[:-1].split('\t')[index_p]
            if p_temp == 'NA':
                p = 1
            else:
                p = float(p_temp)
            if p < padj_thre and log2FC > logFC_thre:
                list_gene_colors.append(sns.xkcd_rgb["dark red"])
                print(line[:-1].split('\t')[begin_samples:])
            elif p < padj_thre and log2FC <= -logFC_thre:
                list_gene_colors.append(sns.xkcd_rgb["marine blue"])
            else:
                list_gene_colors.append(sns.xkcd_rgb["grey"])
        else:
            begin_samples = line[:-1].split('\t').index('padj') + 1
            index_log2FC = line[:-1].split('\t').index('log2FoldChange')
            index_padj = line[:-1].split('\t').index('padj')
            for sample in line[:-1].split('\t')[begin_samples:]:
                list_samples.append(sample)
                if sample in list_pos:
                    list_sample_colors.append(sns.xkcd_rgb["dark red"])
                else:
                    list_sample_colors.append(sns.xkcd_rgb["marine blue"])
        number_row += 1
    print(len(list_sample_colors))
    print(len(list_gene_colors))
    data = pd.DataFrame(data = list_value_list, index = list_ids, columns = list_samples)
    print(data)
    ax = sns.clustermap(data, standard_scale = 1, row_colors = list_gene_colors, col_colors = list_sample_colors)
    plt.setp(ax.ax_heatmap.xaxis.get_majorticklabels(), rotation = 45)
    if save:
        plt.savefig(output_file)
    if show:
        plt.show()
    plt.close()

