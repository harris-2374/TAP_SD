# TAP_SD
Toolset for the Analysis of Pooled Sequenced Data

Created by: Andrew Harris

__________________________________________________________________________________________________________________________________________

Introduction

	The purpose of this toolset is to aid in the analysis of pooled sequenced data, and to provide a set of tools that will enable the user to better organize their research of genetic variation. The toolset is divided into two parts: a pipeline for variant analysis, and several scripts that are useful outside of the pipeline. This toolset is intended to be used in conjunction with CRISP, GATK, and Popoolation. 

__________________________________________________________________________________________________________________________________________

Required Modules
- numpy
- matplotlib
- pandas
- scipy
- csv
- argparse 
- itertools
- PyVCF
__________________________________________________________________________________________________________________________________________

This toolset picks up after the creation of .bam files. For the creation of this toolset, .bam files were created using SpeedSeq, but any alignment program that creates .bam files should work. 

The pipeline is designed so the output of one program acts as part of the input of the next program.

__________________________________________________________________________________________________________________________________________
Implementation:

- This toolset is implemented in Python, utilizing: Popoolation, CRISP, GATK, and SnpEff. 
__________________________________________________________________________________________________________________________________________

Running the Pipeline:


	Fst Z-score transformer:
	  Summary:
	    - This program takes in the output .fst file, produced by Popoolation, and calculates transformed Z-scores for each Fst value. It also provides an interactive distribution plot, using Matplotlib. The distribution graph may be used to help determine the threshold to use in the next step of the pipeline. 

	  Command Line Requirements:
	    - Input: (-I) (--input) - (.fst) file created by Popoolation.
	    - Output: (-o) (--output) - Pathway to output file. 

	  Example:
	    python3 Fst_Z-score_transformer.py --input <pathway_to_fst_file> --output <pathway_to output_folder>

	  Important Notes:
	    - This program is able to handle .fst files with more than two pools of data. If there are more than two pools present, more 		than one output file will be produced. Each file will contain Fst and Z-score data for a singular comparison.  
	    - In reference to the output files produced, “1_to_2” indicates the “1:2” comparison. This rule is consistent for all output files.  

__________________________________________________________________________________________________________________________________________

	Up Down Stream Gene Identification:

	  Summary:
	    - This program determines if a given window, with a transformed Z-score above a given threshold, is within a gene. It also indicates the closest up and down stream gene, with the distance from each. This is based off the organisms reference (.gff/.gtf) file.

	  Command Line Requirements:
	     - Input: (-I) (--input) - (.bed) file created by Fst Z-score transformer.py.
	     - Output: (-o) (--output) - Pathway to output file and name you would like to call the file.
	     - Reference: (-r) (--reference) – Pathway to reference file.
	     - Threshold: (-t) (--threshold) – Minimum standard deviation to analyze. 

	  Example:
	    python3 Up Down Stream Gene Identification.py --input <pathway_to_fst_file> --output <pathway_to_output file.txt> --reference <pathway_to_reference_file> --threshold <wanted_threshold>

	  Important Notes:
	    - Ensure that the Chromosome names are identical to the names of the reference file. 
	    - In reference to the output files produced, “1_to_2” indicates the “1:2” Fst comparison. This rule is consistent for all output files created.

__________________________________________________________________________________________________________________________________________

	Truncated Region Creator:
	  Summary:
	    - Utilizing the output text file of Up_Down_Stream_Gene_Identification.py, the truncated region creator will produce a new tab delimited file that contains regions of continuous window reads of Z-score data above the indicated threshold. This file is best used in Excel as an organization tool when beginning your analysis of genes of interest. 

	  Command Line Requirements:
	    - Input: (-I) (--input) - (.txt) file created by Up_Down_Stream_Gene_Identification.py.
	    - Output: (-o) (--output) - Pathway to output file location - indicate file name, but not file format. (i.e. '.txt') 
	    - Step: (-s) (--step) – Step size used to create Fst calculations in Popoolation. (i.e. 1000)

	  Example:
	    python3 Truncated_Region_Creator.py --input <pathway_to_fst_file> --output <pathway_to_output_folder> --step <numerical_step_value>

__________________________________________________________________________________________________________________________________________

	Annotated VCF Separator:

	  Summary:
	    - Takes annotated VCF file and pulls variants in regions of interest, then separates them into mis-sense mutations and HIGH effect mutations into one file, and the rest of the variants within the given region into another file. These will be separated into “WANTED” (mis-sense and HIGH variants) and “OTHER” (all other variants) folders found in the main “OUTPUT” folder. 

	  Command Line Requirements:
	    - Input: (-I) (--input) – Pathway to annotated .vcf file.
	    - Output: (-o) (--output) - Pathway to output FOLDER. (Note: Do not indicate a file name, simply the location)
	    - Region: (-r) (--region) – Pathway to truncated (.txt) region file created by Truncated_Region_Creator.py.

	  Example:
	    python3 Annotated_Variant_Separator.py --input <pathway_to_fst_file> --output <pathway_to_output_folder> --region <pathway_to_region_file>
	    
	    Important Notes:
	    	- Ensure that the annotated variant VCF file is bgzipped, and then indexed using tabix. This is ESSENTIAL to properly run this program. Not doing this step will throw errors.

__________________________________________________________________________________________________________________________________________


Extra Scripts:

	Genotype Expander:
	  Summary:
	    - The Genotype Expander takes vcf files with pooled genotypic data and separates the data into individual pseudo-samples. This program is built to take in pooled data of any size, and any ploidy. 

	  Command Line Requirements:
	    - Input: (-I) (--input) – Pathway to annotated .vcf file.
	    - Output: (-o) (--output) - Pathway to output FOLDER. (Note: Do not indicate a file name, simply the location)
	    - Ploidy: (-p) (--ploidy) – The ploidy of the organism’s (**One value only**)

	  Example:
	      python3 Genotype_Expander.py --input <pathway_to_fst_file> --output <pathway_to_output_file> --ploidy <ploidy_of_samples>

__________________________________________________________________________________________________________________________________________

	Allele Frequency and Z-score Interactive Plot:

	  Summary:
	    -	Produces interactive graphs to compare Z-score and allele frequencies. 

	  Command Line Requirements:
	    - Input: (-I) (--input) – Indicate the pathway to input .bed file created with Fst_Z_score_transformer.py.
	    - Frequency: (-f) (--frequency) - Indicate pathway to ‘_pwc’ delta allele frequency file created by Popoolation.
	    - Reference: (-r) (--reference) – Provide pathway to reference ‘.gff’ or ‘.gtf’ file.
	    - Threshold: (-t) (--threshold) – Indicate the SAME Z-score threshold used in Up-Down Stream Gene Identification.
	    - Chromosome: (-c) (--chromosome) – Indicate which chromosome(s) to produce. (i.e. ‘NC_009144.3’ for single chromosome, ‘NC_009144.3,NC_009154.3’for a list of chromosomes, or ‘all’ for all chromosomes)

	  Example:
	    python3 AlleleFrequency_Zscore_Interactive_Plot.py --input < pathway_to_fst_file > --frequency < pathway_to_pwc_file > --threshold <Z_score_threshold_used > --chromosome < ‘all’, ‘NC_009144.3, NC_009154.3’, or ‘NC_009144.3’ >

	  Important Notes:
	    - This program loads one chromosome at a time, to continue to the next chromosome, you must close the current window. 
	    - Due to the large size of most data sets, the graphs may take a few minutes to load. 
    
    
   _________________________________________________________________________________________________________________________________________
    
    Contributors:
    	- Brian W. Davis, Ph.D.
