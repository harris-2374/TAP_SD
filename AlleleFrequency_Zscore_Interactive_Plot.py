import pandas as pd
import matplotlib.pyplot as plt
# import numpy as np
import argparse



"""

NOTES:
    1. Due to this program taking in large data sets, it may take a couple of minutes to produce and open the graphs.


"""

parser = argparse.ArgumentParser(description="Provides interactive Chromosome graphs with Z-score and Delta allele frequencies.")
parser.add_argument('-I', '--input', required=True, help='Please indicate the pathway to input .bed file created with "1_Fst_to_Z-score_transformer".')
parser.add_argument('-f', '--frequency', required=True, help="Indicate pathway to _pwc delta allele frequency file provided from Popoolation2.0.")
parser.add_argument('-r', '--reference', required=True, help='Provide pathway to reference .gff file. NOTE: This program currenlty only supports .gff, not .gtf.')
parser.add_argument('-t', '--threshold', required=True, help='Indicate the minimum transformed Z-score used in Gene identification.')
parser.add_argument('-c', '--chromosome', required=True, help="'all' = every chromosome,i.e 'NC_009144.3, NC_009145.4' = graphs for indicated chromosomes, or give name for singular chromosome. ")
args = parser.parse_args()


#-----------------------------------------------------------------------------------------------------------------------
#DATA PREP
#
def main_call(ref1, freq, input1, thresh, chrom_choice):
    ref_read = pd.read_csv(ref1, sep='\t', engine='python', index_col=False, names=['Chromosome', 'Feature', 'Type', 'Start', 'End', 'Score', 'Strand', 'Frame', 'Attribute'])
    ref_df = pd.DataFrame(data=ref_read)
    ref_df_cut = ref_df[8:]
    ref_df_cut_drop = ref_df_cut.drop(['Feature', 'Score', 'Strand', 'Frame', 'Attribute', 'Type'], axis=1)

    ref_start = ref_df_cut_drop['Start']

    ref_end = ref_df_cut_drop['End']

    Chromlist = []

    Chromlist_lowstart_highend = []

    for chrm in ref_df_cut['Chromosome']:
        if chrm in Chromlist:
            continue
        elif '#' in chrm:
            continue
        elif 'NW_' in chrm:
            continue
        else:
            Chromlist.append(chrm)

    for chrm1 in Chromlist:
        holder = 0
        for place in ref_df_cut_drop['Chromosome']:
            temp = []
            if chrm1 is place:
                temp.append(chrm1)
                strt = list(ref_start)
                temp.append(strt[holder])
                eend = list(ref_end)
                temp.append(eend[holder])
                Chromlist_lowstart_highend.append(temp)
                break
            else:
                holder += 1
                continue

    data_read = pd.read_csv(input1, sep='\t', index_col=False)
    data_df = pd.DataFrame(data=data_read)

    data_start = list(data_df['Start'])
    data_zscore = list(data_df['Z-score'])

    freq_read = pd.read_csv(freq, sep='\t', index_col=False)
    freq_df = pd.DataFrame(data=freq_read)

    freq_pos = list(freq_df['pos'])
    freq_diff = list(freq_df['diff:1-2'])




    #-----------------------------------------------------------------------------------------------------------------------

    longest_chrom = 0

    for chrom1 in Chromlist_lowstart_highend:
        if float(chrom1[2]) > float(longest_chrom):
            longest_chrom = chrom1[2]
            continue

    highzscr = 0.0

    for zscr in data_zscore:
        if zscr == 'na':
            continue
        else:
            if float(zscr) > float(highzscr):
                highzscr = float(zscr)
                continue
            else:
                continue

    highzscr2 = highzscr + 2.0

    def scatter(x, y, titl, x2, y2):
        f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
        ax2.scatter(x2, y2, s=5)
        ax1.scatter(x, y, s=5)
        ax1.axhline(float(thresh), 0, highzscr2, color='r')
        # ax2.axhline(0.75, 0, highzscr2, color='r')
        ax1.set_title(titl)
        ax1.set_ylim(0, float(highzscr2))
        ax2.set_ylim(0, 1.1)
        ax1.set_xlim(0, (longest_chrom + 1.0))
        ax2.set_xlim(0, (longest_chrom + 1.0))
        ax1.set_ylabel('Z-score')
        ax2.set_ylabel('Allele Freq. Difference')
        plt.xlabel('Chromosomal Position')
        plt.show()

    def datasetup(chrom):
        pos = []
        pos2 = []
        zscore = []
        freq1 = []
        holder2 = 0
        holder3 = 0
        for j in data_df['Chromosome']:
            if str(chrom) in str(j):
                if data_zscore[holder2] == 'na':
                    zscore.append(0)
                    pos.append(data_start[holder2])
                    holder2 += 1
                    continue
                elif float(data_zscore[holder2]) < 0:
                    zscore.append(0)
                    pos.append(data_start[holder2])
                    holder2 += 1
                    continue
                else:
                    zscore.append(float(data_zscore[holder2]))
                    pos.append(data_start[holder2])
                    holder2 += 1
                    continue
            else:
                holder2 += 1
                continue
        for j in freq_df['##chr']:
            if str(chrom) in str(j):
                if freq_diff[holder3] == 'na':
                    freq1.append(0)
                    pos2.append(freq_pos[holder3])
                    holder3 += 1
                    continue
                elif float(freq_diff[holder3]) < 0:
                    freq1.append(0)
                    pos2.append(freq_pos[holder3])
                    holder3 += 1
                    continue
                else:
                    freq1.append(float(freq_diff[holder3]))
                    pos2.append(freq_pos[holder3])
                    holder3 += 1
                    continue
            else:
                holder3 += 1
                continue
        scatter(pos, zscore, chrom, pos2, freq1)

    which_chrom = chrom_choice

    if which_chrom == 'all':
        for p in Chromlist:
            datasetup(p)
    elif ',' in which_chrom:
        for jj in list(which_chrom):
            print(jj)
            datasetup(jj)
    else:
        datasetup(which_chrom)


if __name__ == '__main__':
    main_call(args.reference, args.frequency, args.input, args.threshold, args.chromosome)
