import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as stats
import math
import csv
import argparse

# ----------------------------------------------------------------------------------------------------------------------

parser = argparse.ArgumentParser(description="Takes in a .fst file, performs Z-score transformation, and provides "
                                             "interactive distribution plot of Z-score data.")
parser.add_argument('-I', '--input', required=True, help='Please indicate the pathway to input .fst file.')
parser.add_argument('-o', '--output', required=True, help="Indicate pathway to output location")
args = parser.parse_args()

# ----------------------------------------------------------------------------------------------------------------------


def maincall_output_mult(input4, output4):
    DSLDcsv = input4
    name_of_fst = (list(DSLDcsv)[-1])
    fst = DSLDcsv[str(name_of_fst)]
    fst_list = []
    for current_fst in fst:  # Obtain all Fst values, ignore 'na' values.
        if 'na' in current_fst:
            continue
        else:
            fst_list.append(current_fst[4:])

    def popmean(x):  #Calculate population mean
        total = 0
        for value in x:
            total += float(value)

        return total*1.0/len(x)

    def popstand(x):  #Calculate Population Standard Deviation
        length = len(x)
        m = popmean(x)
        total_sum = 0
        for current_length in range(length):
            total_sum += (float(x[current_length])-float(m))**2
        under_root = total_sum*1.0/length
        return math.sqrt(under_root)

    def zscores(x, mean, stand):  #Calculates Z-Score for a given Fst value
        score = ((float(x)-float(mean))/float(stand))
        return score

    def MainCall():
        score_list = []
        for g in DSLDcsv[name_of_fst]:
            if g[4:] == 'na':
                score_list.append('na')
            else:
                score_list.append(zscores(g[4:], popumean, popustand))

        return score_list

    def PlotCall():
        score_list = []
        for g in DSLDcsv[name_of_fst]:
            if g[4:] == 'na':
                continue
            else:
                score_list.append(zscores(g[4:], popumean, popustand))

        return score_list

    popumean = popmean(fst_list)
    print('popmean = ', popumean)
    popustand = popstand(fst_list)
    print('popstan = ', popustand)
# ----------------------------------------------------------------------------------------------------------------------
    ## Data organization for output ##

    Chromo_list = ['Chromosome', ]
    Chromo_pos = ['Position', ]
    FSTSCORE = ["Fst_Score", ]
    Z_score = MainCall()
    Z_scoretitle = ['Z-score', ]

    for score in DSLDcsv[name_of_fst]:
        FSTSCORE.append(score[4:])

    for name in DSLDcsv[str((list(DSLDcsv)[0]))]:
        Chromo_list.append(name)

    for pos in DSLDcsv[str((list(DSLDcsv)[1]))]:
        Chromo_pos.append(pos)
    for tit in Z_score:
        Z_scoretitle.append(tit)

    zip(Chromo_list, Chromo_pos, FSTSCORE, Z_scoretitle)


    output_new = output4

    if '.bed' not in output_new:
        raise TypeError(print("Output file does not appear to be in .bed format. "
                              "Please re-run program with corrections."))

    with open(output_new, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(zip(Chromo_list, Chromo_pos, FSTSCORE, Z_scoretitle))

# ----------------------------------------------------------------------------------------------------------------------

    zscore_plot = PlotCall()  #Next steps create distribution plot of all Z-scores.
    zscore_plot.sort()
    zmean = np.mean(zscore_plot)
    zstand = np.std(zscore_plot)
    pdf = stats.norm.pdf(zscore_plot, zmean, zstand)
    plt.title(str(name_of_fst))
    plt.plot(zscore_plot, pdf)
    plt.show()

# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------


def maincall2(input3, output3):

    names1 = ['Chromo', 'Position', '3', '4', 'dd']
    file_read = pd.read_csv(input3, sep='\t')

    for i in file_read:  # i is not used because the first line is on
        name_holder = []
        tag_hold = []
        file_list = list(file_read)
        name_holder = file_list[5:]
        for i2 in name_holder:
            tag_hold.append(i2[0:3])
            continue
        for j in tag_hold:
            temp_hold = 'Fst_' + j
            names1.append(temp_hold)
            continue
        break
    file_read = [] #Empty memory of temp table.
    new_read = pd.read_csv(input3, sep='\t', names=names1)
    file_df = pd.DataFrame(data=new_read)
    name_crop = int(len(names1))
    for i4 in range(5, name_crop):
        place = names1[i4]
        print(place)
        tag_holder1 = tag_hold[(i4-5)]  # Ignore warning, 'tag_hold' is made during data processing.
        out_fix = output3 + '/Comparison_' + str(tag_holder1[0]) + '-' + 'to' + '-' + str(tag_holder1[2]) + '_Zscores.bed'
        maincall_output_mult(file_df[['Chromo', 'Position', '3', '4', 'dd', str(place)]], out_fix)
    return


if __name__ == '__main__':
    maincall2(args.input, args.output)

