import os
import errno
import vcf
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Produces two folders of organized .vcf region files. One folder "
                                             "contains gene regions with only mis-sense and HIGH mutations only, "
                                             "other has all other variants. ")
parser.add_argument('-I', '--input', required=True, help='Please indicate the pathway to annotated .vcf file. ')
parser.add_argument('-o', '--output', required=True, help="Indicate pathway to output location, DO NOT INDICATE "
                                                          "FILE NAME. This program produces two folders with files "
                                                          "named based off their respective regions.")
parser.add_argument('-r', '--region', required=True, help='Provide pathway to reference .txt Gene Region file.')
args = parser.parse_args()


def main(region, vcf_file, output):

    reg_file = pd.read_csv(region, sep='\t', dtype={"Region": int, "Chromosome": str, "Start": int, "End": int, "Genes": str, "Max_Fst": float, "Max_Zscore": float, "Max_Position": int})
    vcf_reader2 = vcf.Reader(filename=vcf_file, compressed='gz')
    collect_folder = str(output + "/OUTPUT/")
    wanted = str(collect_folder + "WANTED/")
    other = str(collect_folder + 'OTHER/')

    if not os.path.exists(os.path.dirname(wanted)):
        try:
            os.makedirs(os.path.dirname(wanted))
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise

    if not os.path.exists(os.path.dirname(wanted)):
        try:
            os.makedirs(os.path.dirname(wanted))
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    if not os.path.exists(os.path.dirname(other)):
        try:
            os.makedirs(os.path.dirname(other))
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise

    for i3 in range(0, len(reg_file["Region"])):
        chrom = reg_file['Chromosome']
        start1 = reg_file['Start']
        end1 = reg_file["End"]
        reg = reg_file["Region"]
        reg1 = str(reg[i3])

        if len(reg1) == 1:
            reg1 = '00' + reg1
        elif len(reg1) == 2:
            reg1 = "0" + reg1

        wanted_name2 = wanted + "/" + "REG" + str(reg1) + "_" + str(chrom[i3]) + '_' + str(start1[i3]) + '_' + str(end1[i3]) + '.ANN.vcf'
        other_name2 = other + "/" + "REG" + str(reg1) + "_" + str(chrom[i3]) + '_' + str(start1[i3]) + '_' + str(end1[i3]) + '.ANN.vcf'

        wanted_file = vcf.Writer(open(wanted_name2, 'w'), vcf_reader2)
        unwanted_file = vcf.Writer(open(other_name2, 'w'), vcf_reader2)
        try:
            for record in vcf_reader2.fetch(chrom[i3], start1[i3], end1[i3]):
                miss = 'missense_variant'
                comp_string = str(record.INFO['ANN'])
                if miss in comp_string:
                    wanted_file.write_record(record)
                    continue
                elif 'HIGH' in comp_string:
                    wanted_file.write_record(record)
                    continue
                else:
                    unwanted_file.write_record(record)
                    continue
        except ValueError:
            pass


def region_number_fix():
    reg_read = pd.read_csv('temp_region.txt', delimiter='\t')
    reg_df = pd. DataFrame(reg_read)
    reg_nums = reg_df['Region']
    current_number = 1
    placement = 0
    for char in reg_nums:
        if char == current_number:
            print(char, current_number)
            current_number += 1
            placement += 1
            continue
        else:
            print(char, current_number)
            reg_nums.replace(to_replace=char, value=current_number, inplace=True)
            current_number += 1
            placement += 1
            continue
    reg_df.to_csv('temp_region_final.txt', sep='\t', index=False)
    return


def na_check(orig_region):
    trip = 0
    with open(orig_region, 'r') as f:
        with open('temp_region.txt', 'w') as output_file:
            for line in f:
                if 'na' in line:
                    trip += 1
                    continue
                elif 'NA' in line:
                    trip += 1
                    continue
                elif 'nan' in line:
                    trip += 1
                    continue
                else:
                    output_file.writelines(line)
    if trip != 0: #if there are no na values, it will skip the correction of region numbers
        region_number_fix()
    return


if __name__ == '__main__':
    na_check(args.region)
    main('temp_region_final.txt', args.input, args.output)


