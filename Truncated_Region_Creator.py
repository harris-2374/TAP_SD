import pandas as pd
import csv
import itertools
import argparse

parser = argparse.ArgumentParser(description="Takes Gene list and truncates list into regions of  continuoud windows "
                                             "with Z-scores above threshold.")
parser.add_argument('-I', '--input', type=str, required=True, help='Please indicate the pathway to input .txt gene file created '
                                                         'with "Up_Down_Stream_Gene_Identification.py"')
parser.add_argument('-o', '--output', type=str, required=True, help="Indicate pathway to output location, including the output "
                                                          "file name w/ .txt as file type.")
parser.add_argument('-s', '--step', type=int, required=True, help='Indicate the step size of the .fst file as a numeric value '
                                                        '(i.e 1000)')

args = parser.parse_args()
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

def main_call(step, input1, output):
    data_read = pd.read_csv(input1, sep='\t')
    data_df = pd.DataFrame(data=data_read)
    region = ['Region', ]
    max_fst = ['Max_Fst', ]
    max_zscore = ['Max_Zscore', ]
    max_pos = ["Max_position"]
    chrom = ['Chromosome', ]
    start = ['Start', ]
    end = ['End', ]
    genes = ['Genes', ]
    win_len = []

    def getstartandstopandchrom():
        holder1 = 0
        curr = 0
        place = 0
        regcall = 1
        genelist = []
        win_len1 = 0
        for strt in data_df['SNP_Pos.']:
            try:
                if holder1 == 0:  # sets start position
                    start.append(strt)
                    holder1 += 1
                    curr = strt
                    chrompos = data_df['Chromosome']
                    chrom.append(chrompos[place])
                    region.append(regcall)
                    gene1 = data_df['Within_gene']
                    if gene1[place] not in genelist:
                        genelist.append(gene1[place])
                    else:
                        continue
                    win_len1 += 1
                    place += 1
                    continue
                elif (strt - curr) != step:       # sets end position
                    chrompos2 = data_df['Chromosome']
                    throws = ['nan', 'NA', '']
                    if genelist in throws:
                        pass
                    else:
                        genes.append(genelist)
                    genelist = []
                    end.append(curr)
                    start.append(strt)
                    chrom.append(chrompos2[place])
                    place += 1
                    gene3 = data_df['Within_gene']
                    try:
                        if gene3[place] not in genelist:
                            genelist.append(gene3[place])
                        else:
                            continue
                    except KeyError:
                        continue
                    win_len.append(win_len1)
                    win_len1 = 1
                    regcall += 1
                    region.append(regcall)
                    curr = strt
                    continue
                else:
                    win_len1 += 1
                    curr = strt
                    gene2 = data_df['Within_gene']
                    if gene2[place] not in genelist:
                        genelist.append(gene2[place])
                        place += 1
                    else:
                        place += 1
                        continue

                continue
            except KeyError:
                continue
        win_len.append(win_len1)
        throws = ['nan', 'NA', '']
        if genelist in throws:
            pass
        else:
            genes.append(genelist)
        end.append(curr)

    def getfst_zscore(df):
        placehold = 0
        z_scr = df['Z-score']
        f_val = df['Fst_score']
        pos1 = df['SNP_Pos.']
        for ii in win_len:
            temp2 = list(z_scr[placehold:placehold + ii])
            f_temp = list(f_val[placehold:placehold + ii])
            p_temp = list(pos1[placehold:placehold + ii])
            top1 = 0.0
            position = 0
            fst = 0.0
            for jj in enumerate(temp2):
                if jj[1] > top1:
                    top1 = jj[1]
                    fst = f_temp[jj[0]]
                    position = p_temp[jj[0]]
                else:
                    continue
            max_zscore.append(top1)
            max_fst.append(fst)
            max_pos.append(position)
            placehold += ii

    ### CALLS ###
    getstartandstopandchrom()
    getfst_zscore(data_df)
    itertools.zip_longest(region, chrom, start, end, genes, max_fst, max_zscore, max_pos)

    with open(output, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(itertools.zip_longest(region, chrom, start, end, genes, max_fst, max_zscore, max_pos))


if __name__ =='__main__':
    main_call(args.step, args.input, args.output)
