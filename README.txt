DNA_jets contains the code used to characterise jets.
Authors: Rosalba Garcia-Millan and Gunnar Pruessner

0. List of files:
gplist_util.h
gplist_util2.h
data_structure.h
gplist_util.c
gplist_util2.c
DNAsliding_window.c
DNA_expected.c 
DNA_panoramic.c
DNA_stencil.c
panoramic_mode.c
pseudoBAMfile.c


1. Compilation

cc -o DNAsliding_window DNAsliding_window.c gplist_util.c

cc -O3 -Wall -o BAM pseudoBAMfile.c gplist_util2.c -lm -lgsl -lgslcblas

cc -O3 -Wall -o exp DNA_expected.c gplist_util2.c -lm -lgsl -lgslcblas

cc -O3 -Wall -o Jets DNA_panoramic.c gplist_util2.c rgm_Gfit.c -lm -lgsl -lgslcblas

cc -O3 -Wall -o max panoramic_mode.c -lm -lgsl -lgslcblas

cc -O3 -Wall -o JetsSten DNA_stencil.c gplist_util2.c -lm -lgsl -lgslcblas


2. Pipeline
DNAsliding_window: Generate bin files from dat files with entries "x  y  count"

pseudoBAMfile: Generate file with total counts for each position x (in case BAM files are not provided).

DNA_expected: Generate file with the average count for a given distance between two positions (in case expected data is not provided).

DNA_panoramic: Generate panoramic curves.

panoramic_mode: Characterise orientation, and strength of each jet.

DNA_stencil: Generate stencil curves, characterise distances reached by the jet.


3. Notes
Initially, the struct data_struct called data in the header file data_structure.h needs the following variables to be assigned.
    char *ident; // name of the dataset, e.g. "WTR1"
    char *filenameBAM; //path to the BAM file (either provided or generated using pseudoBAMfile). 
    char *basenameEXP; //path to file with "the expected data" (either provided or generated using DNA_expected).
    char *basenameHIC; //path to bin file generated with DNAsliding_window
    double total_counts; //total number of counts in dataset



While running the pipeline above, the following path directories in the header file data_structure.h need to be updated.

	bed_file: contains chromosome locations to be analysed
	panoramic_file: contains panoramic curves output from DNA_panoramic.
	jet_analysis_file: contains output from panoramic_mode.
	jet_consensus_file: list of locations considered to have a jet.
	jet_reach_file: contains output from DNA_stencil.



