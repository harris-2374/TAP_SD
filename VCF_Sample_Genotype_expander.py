import pandas as pd
import argparse

#----------------------------------------------------------------------------------------------------------------------#
parser = argparse.ArgumentParser(description="Takes in vcf file and expands genotypes for further analysis.")
parser.add_argument('-I', '--input', type=str, required=True, help='Please indicate the pathway to input vcf file.')
parser.add_argument('-o', '--output', type=str, required=True, help="Indicate pathway to output locaiton, including "
                                                                    "what you would like to call the output.")
parser.add_argument('-p', '--ploidy', type=int, required=True, help='Indicate the ploidy of the organism')
args = parser.parse_args()

#----------------------------------------------------------------------------------------------------------------------#


def mainfunction(vcf_file, output, ploidy):

    def header_fun():
        with open(vcf_file, mode='r', encoding='utf-8') as vcf_reader:

            header_len = 0
            header1 = []

            # print(vcf_reader)
            try:
                insert_flag = 0
                for line in vcf_reader:
                    # print(line)
                    if insert_flag == 2:
                        header1.append("##Genotype_Expander - Andrew Harris" + '\n')
                    if '##' in str(line):
                        # print(line)
                        # print("here")
                        header_len += 1
                        insert_flag += 1
                        header1.append(line)
                        continue
                    else:
                        return header_len, header1
            except Exception:
                vcf_reader.close()

                return header_len, header1

    head_call = header_fun()
    # print(head_call)
    header_length = head_call[0]
    header = head_call[1]

    # print("Header length:", header_length)

#-----------------------------------------------------------------------------------------------------------------------
    def get_sample_amount(headlist):
        """
        :param headlist: VCF Headers and Sample names
        :return: The number of Samples, and a list of the sample names
        """
        sample_names = []
        for pos in headlist[9:]:
            temp = []
            temp.append(pos)

            for i6 in vcf_df[pos]:
                tally = 0
                for j in i6:
                    if j == '/':
                        tally += 1
                    else:
                        continue
                number_of_samps = ((tally + 1) / int(ploidy))
                temp.append(int(number_of_samps))
                break
            sample_names.append(temp)
        # print("Get_samp_amount Done")


        return sample_names

    def name_creator(samp_info):
        """
        :param samp_info: return of get_sample_amount, gives number of samples and the sample names.
        :return: list of list that contains the new names for the samples
        """

        names = samp_info
        names_fixed = []
        for i5 in names:
            temp = []
            for j in range(1, i5[1] + 1):
                temp1 = str(i5[0]) + "_PseudoSample_" + str(j)
                temp.append(temp1)
            names_fixed.append(temp)
        # print("Name_creat Done")
        return names_fixed

    def genotype_handler(df_column, newnames, current, ploid):
        """

        :param df_column: Singular column of genotype data from a ONE sample
        :param newnames: list of new sample names
        :return: Dictionary of single sample genotype breakdown
        """
        # print(df_column)
        dicts = []
        strt = 0
        stp = (ploid + 1)
        for j in newnames:
            try:
                if current in j:
                    temp = [j]
                    for ii in df_column:
                        temp.append(str(ii[strt:stp]))

                    dicts.append(temp)
                else:
                    continue
            except:
                continue
            strt = (stp + 1)
            stp += 4
        # print("Genotype Handler Done")
        return dicts

    def genotype_dictionary_maker(df, names1, newnames, ploid):
        dictionaries = []
        placehold = 0
        for i3 in names1:
            dictionaries.append(genotype_handler(df[i3[0]], newnames[placehold], i3[0], ploid))
            placehold += 1
            continue
        return dictionaries

    def newdataframe(geno_info):
        for pos2 in geno_info:
            for j in pos2:
                se = pd.Series(j[1:])
                vcf_df[j[0]] = se.values
        # print("Genotype Dictionary maker done")
        return

#----------------------------------------------------------------------------------------------------------------------#
    header_temp = header
    # print(header_temp)
    with open(output, 'w') as vcf:
        for iqq in header_temp:
            vcf.write(str(iqq))

    holder = str(vcf_file)

    df_read = pd.read_table(holder, skiprows=header_length, sep='\t', encoding="ISO-8859-1", low_memory=True)


    start = 0
    end = 1000
    popp = 1

    while popp > 0:
        try:
            vcf_df = pd.DataFrame(data=df_read[start:end])

            orig_names = get_sample_amount(list(vcf_df))

            new_names = name_creator(orig_names)

            geno1 = genotype_dictionary_maker(vcf_df, orig_names, new_names, ploidy)

            newdataframe(geno1)

            for i34 in orig_names:
                hold1 = str(i34[0])
                vcf_df.drop(hold1, axis=1, inplace=True)

            vcf_df.to_csv(output, sep='\t', mode='a', encoding='utf-8', index=False)
            start = end + 1
            end += 1001
            # print('Done With That One!')
            # print(start)
            # print(end)
        except:
            break
    return


if __name__ == '__main__':
    mainfunction(args.input, args.output, args.ploidy)
