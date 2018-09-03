import pybedtools
import pandas as pd
import csv
import itertools
import argparse
import time

# ----------------------------------------------------------------------------------------------------------------------
# Argparse #

parser = argparse.ArgumentParser(description="Takes .bed output file, and determines the windows position within the genome (i.e intergenic or specific gene).")
parser.add_argument('-I', '--input', required=True, help='Indicate pathway to input .bed file created with "1_Fst_to_Z-score_transformer".')
parser.add_argument('-o', '--output', required=True, help="Indicate pathway to output location, including the output filename.")
parser.add_argument('-r', '--reference', required=True, help='Provide pathway to reference .gff file. NOTE: This program currenlty only supports .gff, not .gtf.')
parser.add_argument('-t', '--threshold', required=True, help='Indicate the minimum transformed Z-score to be included in analysis. (i.e 7.5)')
args = parser.parse_args()

# ----------------------------------------------------------------------------------------------------------------------


def main_call(input1, output1, ref1, thresh1):
    total_windows = 0
    gff_reader = pd.read_csv(ref1, sep='\t', engine='python', names=['Chromosome', 'Feature', 'Type', 'Start', 'End', 'Score', 'Strand', 'Frame', 'Attribute'])
    print(list(gff_reader))
    gff_dataframe = pd.DataFrame(data=gff_reader)
    print(gff_dataframe.head(10))
    gff_dataframe_cut = gff_dataframe

    def genelist(line):
        '''
        :param line: list of attributes for each entry
        :return: list containing the gene for each entry
        '''
        tag = 'gene='
        gene_list = []
        for i2 in line:
            split_line = str(i2)
            need = split_line.split(';')
            if tag in str(i2):
                for j in need:
                    if tag in j:
                        gene_list.append(j)
                    else:
                        continue
            else:
                gene_list.append('N/A')
        return gene_list

    def mastergenelist():
        """
        
        :return: List of lists containing Chromosome, Start, Stop, and Gene positional data.
        """
        gene_list = list(genelist(gff_dataframe_cut['Attribute']))
        start_list = list(gff_dataframe_cut['Start'])
        stop_list = list(gff_dataframe_cut['End'])
        chromo_list = list(gff_dataframe_cut['Chromosome'])
        mastergenefile2 = []
        for val3 in range(0, len(gene_list)):
            temp = []
            tag = '##'
            if tag not in chromo_list[val3]:
                temp.append(chromo_list[val3])
                temp.append(start_list[val3])
                temp.append(stop_list[val3])
                temp.append(gene_list[val3])
                mastergenefile2.append(temp)
                continue
            else:
                continue
        return mastergenefile2

    mastergenefile = mastergenelist() #variable for recursive calling
    x = pybedtools.example_bedtool(input1) #Open data bed file
    df = pd.read_table(x.fn, low_memory=False)

    def gatherinterest(pos, cutoff):
        '''
        This filters through bed file to get all windows above indicated threshold. 
        :param pos: position of snp
        :param cutoff: z-score cut off
        :return: list of SNP's above threshold
        '''
        sampofint = []
        for interes in enumerate(pos):
            if interes[1] == 'na':
                continue
            else:
                if float(interes[1]) >= cutoff:
                    sampofint.append(interes)
                else:
                    continue
        return sampofint

    def fststartstop(interest):
        '''

        :param interest: list of SNP's above threshold
        :return: list of chromosome, position, and SNP z-score
        '''
        fststartstop_list = []
        for val3 in interest:
            temp = []
            if val3[0] in df['Z-score']:
                temp.clear()
                temp.append(df.iloc[val3[0]]['Chromosome']) 
                temp.append(int(df.iloc[val3[0]]['Position']))
                temp.append(float(val3[1]))
                temp.append(df.iloc[val3[0]]['Fst_Score'])
                fststartstop_list.append(temp)
                continue

        return fststartstop_list

    def downstreamV2(window1, gene_master):
        """
        - Deteremines closest down-stream gene by comparing the Start position of SNP and End Position of Gene.
        :param snp: Windows above threshold.
        :param gene_master: master gene list.
        :return: Closest down-stream gene and the distance from the .
        """
        winner = ['chrom', 'start', 1, 'gene']

        for genemst in gene_master[1:]:
            if "N/A" in genemst:
                continue

            else:
                try:
                    current = genemst
                    snpchrom = window1[0]
                    currchrom = current[0]
                    if snpchrom == currchrom:
                        if window1[1] > genemst[2]:
                            # print(current)
                            if int(current[2]) > int(winner[2]):
                                winner = current
                                continue
                            else:
                                continue
                        else:
                            continue
                    else:
                        continue
                except IndexError:
                    # print('IndexError Downstream',genemst)
                    continue

        return winner

    def upstreamV2(window2, gene_master):
        """
        
        :param snp: SNP's above threshold.
        :param gene_master: Master gene list
        :return: Closest Up-stream gene and the distance from it. 
        """
        winner = ['chrom', 1000000000000000000000000000, 'holder', 'gene']
        for genemst in gene_master[1:]:  #Skips first value because it indicates the length of the entire chromosome, throwing off the program.
            if 'N/A' in genemst:
                continue
            elif genemst is []:
                continue
            else:
                try:
                    current = genemst
                    if window2[0] == current[0]:
                        if float(window2[1]) < genemst[1]:
                            if current[1] < winner[1]:
                                winner = current
                                continue
                            else:
                                continue
                    else:
                        continue
                except IndexError:
                    # print('Index Error Upstream', genemst)
                    continue

        return winner

    def withingene(window3, mstr):
        """

        :param window3: Single window for analysis
        :param mstr: Master Gene list for reference
        """
        upstream4 = upstreamV2(window3, mstr)
        downstream4 = downstreamV2(window3, mstr)

        if upstream4[3] == downstream4[3]:
            return upstream4[3]
        else:
            return 'INTERGENIC_REGION'
    thresh = float(thresh1)
    snp_of_interest = fststartstop(gatherinterest(df['Z-score'], thresh))
# ----------------------------------------------------------------------------------------------------------------------
    ## Gathering all information for output file creation ##
    Chromosome_final_list = ['Chromosome', ]
    Position_final_list = ['SNP_Pos.', ]
    within_final_list = ['Within_gene', ]
    distance_downstream = ['Distance Downstream', ]
    downstream_final_list = ['Downstream Gene', ]
    dis_upstream = ['Distance Upstream', ]
    upstream_final_list = ['Upstream Gene', ]
    fst_final = ['Fst_score', ]
    z_score_final = ['Z-score', ]

    for i in snp_of_interest:
        Chromosome_final_list.append(i[0])

    for j in snp_of_interest:
        Position_final_list.append(j[1])

    for l in range(0, len(snp_of_interest)):
        snp = snp_of_interest[l]
        down = downstreamV2(snp, mastergenefile)
        disdown = (snp[1] - down[2])
        distance_downstream.append(disdown)
        downstream_final_list.append(down[3].strip('gene='))

        snp1 = snp_of_interest[l]
        up = upstreamV2(snp1, mastergenefile)
        disup = (up[1] - snp1[1])
        dis_upstream.append(disup)
        upstream_final_list.append(up[3].strip('gene='))

        snp2 = snp_of_interest[l]
        within = withingene(snp2, mastergenefile)
        within_final_list.append(within.strip('gene='))
        continue

    for jj in snp_of_interest:
        fst_final.append(jj[3])

    for kk in snp_of_interest:
        z_score_final.append(kk[2])

    itertools.zip_longest(Chromosome_final_list, Position_final_list, distance_downstream, downstream_final_list,
                          within_final_list, upstream_final_list, dis_upstream, fst_final, z_score_final)

    with open(output1, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(itertools.zip_longest(Chromosome_final_list, Position_final_list, distance_downstream, downstream_final_list, within_final_list, upstream_final_list, dis_upstream, fst_final, z_score_final))
        f.close()
    total_windows += len(Chromosome_final_list)
    return total_windows


def logfile(total_windows):
    # Indicates the amount of windows found above threshold #
    print()
    print("---------------------------------------------------------------------------------------------------")
    print()
    print("Number of windows found above threshold: ", total_windows)
    print()
    print('---------------------------------------------------------------------------------------------------')
    print()
    print("Run output file through 'Truncated_Region_Creator.py' to truncate data.")
    print()
    print('---------------------------------------------------------------------------------------------------')
    return


if __name__ =='__main__':
    try:
        mainnn = main_call(args.input, args.output, args.reference, args.threshold)
        logfile(mainnn)
    except RuntimeWarning:
        pass
