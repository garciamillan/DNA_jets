/*
 * This code is based on GP's raw2binned.c, which
 * is based on Francesc's freq2dens_lo.c.
 * GP ironed out some bits for better ease of use.
 * RGM has included the initial part regarding to the DNA data processing.
 * date: Thu 21 Feb 2019 19:06:12 GMT (prev version Mon 18 Feb 2019 14:48:16 GMT)
 *
 * $Header: /Users/garciamillan/Documents/uni/ImperialCL/DNAfolding/HiC_2kb/contactdomains/results10_slidingwindow/DNAsliding_window.c
 *
 *
 This program reads a list of chromosomes and their divisions (e.g. compartments or domains). Then, from the input files with the contact data, finds to which
 division belongs each point and decides whether or not it belongs to the ensemble of interest. If it doesn't, the data (count) is disregarded. If it does, then
 the data is processed in the way raw2bin.c works (see below). Now, the output does not only contain the density, the frequency is also calculated. Future step:
 choose ensemble from input argument, include errorbars.
 Compile: gcc -Wall -o dna DNAcontactdensity.c
 
 The header before this edit is: /home/ma/p/pruess/.cvsroot/BD/raw2binned.c,v 1.13 2016/12/01 22:06:49 pruess Exp $
 
 raw2bin.c: This program takes a list of numbers corresponding to frequencies of a variable, or simply avalanche sizes. They don't need to be sorted. It outputs the pdf.
 raw2bin.c: It takes in account discreteness effects as proposed by Corral et al. in [1]. It can also be used for continuous variables simply by by using a resolution of R=0

 raw2bin.c: 20 Apr 2015
 raw2bin.c: Just to clarify:
 raw2bin.c: This program takes a stream of event sizes (not of their frequencies, as it says above) and creates a nizely binned histogram from that data.
 raw2bin.c: The header prior to me writing this message is
 raw2bin.c: Header: /home/ma/p/pruess/.cvsroot/BD/raw2binned.c,v 1.9 2014/07/23 16:03:24 pruess Exp



cc -o DNAsliding_window DNAsliding_window.c gplist_util.c
./DNAsliding_window -w 5000:5000 ../RUN04_contact_domains/CD69negDPCTCFKOR1R2/CD69negDPCTCFKOR1R2ch_17_17_5kb.dat
./DNAsliding_window -w 5000:500000 ../RUN04_contact_domains/CD69negDPCTCFKOR1R2/CD69negDPCTCFKOR1R2ch_17_17_5kb.dat > visited.dat
grep "#VISIT 1 " visited.dat | sed 's/#VISIT //' > visited.dat_1
grep "#VISIT 2 " visited.dat | sed 's/#VISIT //' > visited.dat_2

grep "#VISIT 1 " visited_5000_50000.dat3 | sed 's/#VISIT //' > visited_5000_50000.dat3_1
grep "#VISIT 2 " visited_5000_50000.dat3 | sed 's/#VISIT //' > visited_5000_50000.dat3_2
grep "#VISIT 3 " visited_5000_50000.dat3 | sed 's/#VISIT //' > visited_5000_50000.dat3_3
grep "#VISIT 4 " visited_5000_50000.dat3 | sed 's/#VISIT //' > visited_5000_50000.dat3_4
grep "#VISIT 5 " visited_5000_50000.dat3 | sed 's/#VISIT //' > visited_5000_50000.dat3_5
grep "#VISIT 6 " visited_5000_50000.dat3 | sed 's/#VISIT //' > visited_5000_50000.dat3_6
grep "#VISIT 7 " visited_5000_50000.dat3 | sed 's/#VISIT //' > visited_5000_50000.dat3_7
grep "#VISIT 8 " visited_5000_50000.dat3 | sed 's/#VISIT //' > visited_5000_50000.dat3_8
grep "#VISIT 9 " visited_5000_50000.dat3 | sed 's/#VISIT //' > visited_5000_50000.dat3_9
grep "#VISIT 10 " visited_5000_50000.dat3 | sed 's/#VISIT //' > visited_5000_50000.dat3_10
grep "#VISIT 11 " visited_5000_50000.dat3 | sed 's/#VISIT //' > visited_5000_50000.dat3_11
grep "#VISIT 12 " visited_5000_50000.dat3 | sed 's/#VISIT //' > visited_5000_50000.dat3_12
grep "#VISIT 13 " visited_5000_50000.dat3 | sed 's/#VISIT //' > visited_5000_50000.dat3_13
grep "#VISIT 14 " visited_5000_50000.dat3 | sed 's/#VISIT //' > visited_5000_50000.dat3_14
grep "#VISIT 15 " visited_5000_50000.dat3 | sed 's/#VISIT //' > visited_5000_50000.dat3_15
grep "#VISIT 16 " visited_5000_50000.dat3 | sed 's/#VISIT //' > visited_5000_50000.dat3_16


*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <time.h>
//#include <sys/syslimits.h>
#define PATH_MAX 1024

#include "gplist_util.h"


int spit_out_minmax=0;
int minmax_fixed_by_hand=0;
int data_out_of_minmax=0;
int force_last_bin_correction=0;
int flag_windows=1;
int flag_compartments=0;
int flag_ranked_compartments=0;

//careful with next line, it can cause overwriting
#define MAX_LENGTH_ARRAY 1000
#define MAX_NUM_WINDOWS 5000
#define NUM_MOMS (6)
#define LAST_CHR 24
#define NBINS 21
#define n(i,g) n[((g-1)*kmax+i)]
/* Rosalba: Do not give the macro the same name as the variable. Also, bracket all arguments of the macro -- imagine a call chr_w_moms_distance(i+1,j) */
#define chr_w_moms_distance(i,j) chr_w_moms_distance[(i*NUM_MOMS+j)]

#define MALLOC(a,n) if ((a=malloc(sizeof(*a)*(n)))==NULL) { fprintf(stderr, "Not enough memory for %s, requested %i bytes, %i items of size %i. %i::%s\n", #a, (int)(sizeof(*a)*n), n, (int)sizeof(*a), errno, strerror(errno)); exit(EXIT_FAILURE); } else { printf("# Info malloc(3)ed %i bytes (%i items of %i bytes) for %s.\n", (int)(sizeof(*a)*(n)), n, (int)sizeof(*a), #a); }


struct compartment {
    int group;
    double start, end;
    char label[100];
};

struct chromosome {
    int num_comp, chr_with_input_file;
    char filename[PATH_MAX];
    struct compartment comp[MAX_LENGTH_ARRAY];
};

//find compartment of position x by divide and conquer
void find_compartment(double x, struct chromosome chr, int *ind, int *group, char label[]) {
    int indA, indB, indC;
    double A,B;
    indA=0;
    if(chr.num_comp>=0) indB=chr.num_comp; else exit(EXIT_FAILURE);
    A=chr.comp[indA].start;
    B=chr.comp[indB].end;
    if(B<A){ fprintf(stderr,"error A %lf\tB %lf\n",A,B); exit(EXIT_FAILURE); }
    if(x<A || x>B){ *ind=-1; *group=0; strncpy(label,"none",10);}
    else {
        while (indB-indA>1) {
            indC=(int)(indA+indB)/2;
            if(x<=chr.comp[indC].start) indB=indC;
            else indA=indC;
        }
        if(x<=chr.comp[indA].end){ *ind=indA; *group=chr.comp[indA].group; strncpy(label,chr.comp[indA].label,10); }
        else if(x==chr.comp[indB].start){ *ind=indB; *group=chr.comp[indB].group; strncpy(label,chr.comp[indB].label,10); }
        else{ *ind=-1; *group=0; strncpy(label,"none",10); }
    }
}

/*struct window {
    //double start, end, middle;
    long double count;
};

struct chromosome_windows {
    int num_wind, chr_with_input_file, chr_num;
    double length;
    char filename[PATH_MAX];
    struct window window[MAX_NUM_WINDOWS];
};*/

//find windows where x lies in
/*void find_windows(double x, struct chromosome_windows chr_w, double resolution, int *ind_min, int *ind_max) {
    int i_m, i_M;
    i_m=i_M=(int)x/resolution;
    if(x>=chr_w.window[i_m].start && x<=chr_w.window[i_m].end){
        while(x>=chr_w.window[i_m].start && x<=chr_w.window[i_m].end) i_m--;
        while(x>=chr_w.window[i_M].start && x<=chr_w.window[i_M].end) i_M++;
        *ind_min=i_m++;
        *ind_max=i_M--;
    } else { fprintf(stderr,"Error: error in locating point %lf\n",x); exit(EXIT_FAILURE); }
}*/

int main(int argc, char *argv[]){
    long long int valid, lines, comments, line_spoiled_by_comment, data_in_ensemble;
    char buffer[8192];
    char pattern[8192];
    //strcpy(pattern,"%*[^_]_%*[^_]_%[0-9*XYxy]_%[0-9*XYxy]%*[^0-9XYxy]");
    strcpy(pattern,"%*[^_]_%*[^_]_%*[^_]_%[0-9*XYxy]_%[0-9*XYxy]%*[^0-9XYxy]");
    char identifier1[PATH_MAX];
    char identifier2[PATH_MAX];
    struct chromosome chr[LAST_CHR+1]; //chrs X and Y will be stored in positions LAST_CHR-1 and LAST_CHR
    //struct chromosome_windows chr_w; //only one chromosome
    double chr_length[LAST_CHR+1];
    double x, y, x0, x1;
    double delta_w,Delta_w;
    int chr_w_number_windows;
    int indx, indy, g, group, groupx, groupy;
    int chr_w_chr_num, indx_max,indx_min,indy_max,indy_min;
    char chr_w_filename[PATH_MAX];
    double chr_w_length;
    long double *chr_w_count;
    long double *chr_w_moms_distance; //moments contact distance
    char labelx[]="initialise";
    char labely[]="initialise";
    char chrID[100];
    char compID[100];
    char sign[1];
    int k, chr_num, chr_num_prev;
    int i, j, nbins=5;
    FILE *fp=stdin;
    FILE *fg=stdin;
    double min=0, max=-1.;
    double smin;
    long double data, count; // count is in case there is a columns with "counts" of each "data" but not necessarily forming a histogram
    long double mom;
    int ch;
    //long long int *n;
    long double *n;
    double *s;
    //long long int N;
    long double N;
    long double resolution=1., De,dens,freq;
    int clean=0;
    int flag_identifier_set=0;
    char *p;
    long long int length_provided=-1LL;
    
    spit_out_minmax=0;
    resolution=5000.;
    min=5000; //this is the resolution
    max=20000000;
    delta_w=500000;//50000;
    Delta_w=1500000;//250000;

    setlinebuf(stdout);
    printf("# $/Users/garciamillan/Documents/uni/ImperialCL/DNAfolding/HiC_2kb/contactdomains/results10_slidingwindow/DNAsliding_window.c\n");
    printf("# Command:");
    for (i=0; i<argc; i++) printf(" %s", argv[i]);
    printf("\n");
    {time_t tm;
    tm=time(NULL);
    printf("# Info: %s", ctime(&tm));
    }
    for(i=0;i<=LAST_CHR;i++){ strcpy(chr[i].filename,""); chr[i].num_comp=-2; }
    strcpy(chr_w_filename,"");

    while ((ch = getopt(argc, argv, "C:cf:g:p:l:Mm:n:b:r:h:t:w:")) != -1)
        switch (ch) {
	    case 'C':
	        strcpy(identifier1, optarg);
		{ char *p;
		  for (p=identifier1; *p; p++) if (*p==':') break;
		  if (*p) {
		    strcpy(identifier2,p+1);
		    *p=(char)0;
		  } else {
		    fprintf(stderr, "Failed to recognise chromosome identifiers.\n");
		    exit(EXIT_FAILURE);
		  }
		}
		printf("# Info: identifier1: %s identifier2: %s\n", identifier1, identifier2);
		flag_identifier_set=1;
		break;
            case 'c':
                clean=1;
                break;
            case 'f': //compartment_file
                flag_compartments=1;
                if ((fg=fopen(optarg, "rt"))==NULL) {
                    if (strcmp(optarg, "-")==0) fg=stdin;
                    else {
                        fprintf(stderr, "Cannot open file \"%s\" (%i::%s)\n", optarg, errno, strerror(errno));
                        exit(EXIT_FAILURE);
                    }
                }
                break;
            case 'g': //ranked_compartment_file
                flag_ranked_compartments=1;
                if ((fg=fopen(optarg, "rt"))==NULL) {
                    if (strcmp(optarg, "-")==0) fg=stdin;
                    else {
                        fprintf(stderr, "Cannot open file \"%s\" (%i::%s)\n", optarg, errno, strerror(errno));
                        exit(EXIT_FAILURE);
                    }
                }
                break;
            case 'p': //pattern of raw_contact_data_file(s)
                strncpy(pattern,optarg,50);
                break;
            case 'w': //windows spacing delta_w and length Delta_w
                flag_windows=1;
                if (sscanf(optarg, "%lg:%lg", &delta_w, &Delta_w)!=2) {
                    fprintf(stderr, "Failed to scan delta_w:Delta_w, should be -w 50000:250000.\n");
                    exit(EXIT_FAILURE);
                }
                if (delta_w>Delta_w) {
                    fprintf(stderr, "Failed to scan delta_w:Delta_w, should be -m 50000:250000, delta_w<Delta_w.\n");
                    exit(EXIT_FAILURE);
                }
                break;
	    case 'l':
	    	length_provided=strtoll(optarg, NULL, 10);
		break;
            case 'm':
                if (sscanf(optarg, "%lg:%lg", &min, &max)!=2) {
                    fprintf(stderr, "Failed to scan min:max, should be -m 1.3:42.0.\n");
                    exit(EXIT_FAILURE);
                }
                if (min>max) {
                    fprintf(stderr, "Failed to scan min:max, should be -m 1.3:42.0, min<max.\n");
                    exit(EXIT_FAILURE);
                }
                minmax_fixed_by_hand=1;
                break;
            case 't':
                force_last_bin_correction=1;
                break;
            case 'M':
                spit_out_minmax=1;
                break;
            case 'n':
            case 'b':
                if ((nbins=atoi(optarg))<=0) {
                    fprintf(stderr, "nbins must be strictly positive.\n");
                    exit(EXIT_FAILURE);
                }
                break;
            case 'r':
                if ((resolution=atof(optarg))<0) {
                    fprintf(stderr, "resolution r must be non-negative.\n");
                    exit(EXIT_FAILURE);
                }
                break;
            default:
                fprintf(stderr, "Flag %c unknown.\n", ch);
            case 'h':
                fprintf(stderr, "Usage: ./dna -f filename_compartments -p pattern -M -m min:max -b nbins -r resolution file_contact1 file_contact2 file_contact3 ...\n");
                fprintf(stderr, "-f filename_compartments specified input file (with compartments). By default the program reads from stdin.\n"
                        "-p give the pattern to read the chromosomes from the filename file_contact*. Default example WT_ch_1_3.dat\n"
                        "-M switch to make the program spit out min and max of data only. This is \n"
                        "   useful when reading from stdin. Determine min and max first, then rerun\n"
                        "   with those supplied via...\n"
                        "-m min:max Optional. Defaults to 5000:20000000\n"
                        //"-m min:max Note that for binning from stdin, those _have_ to be given.\n"
                        //"   Otherwise they are determined first, after which the stream is rewound.\n"
                        "-b nbins Number of bins _per_ _decade_. Defaults to 5.\n"
                        "-r resolution gives the resolution for discretisation of bins towards small\n"
                        "   values. Must correspond to the resolution of the data. Defaults to 5000, \n"
                        "   should be set to 0. for continuous data.\n"
                        "-w delta_w:Delta_w Optional. Defaults to 50000:250000\n"
                        "-t Force correction of last bin normalization. To be used when the user knows\n"
                        "   beforehand that the data being suppliend has been truncated, so that the min\n"
                        "   and or max are not statistically natural.\n"
                        " ./dna -f CD69negDPCTCFKOR1R2_eigenvector_KR100kb.bed CD69negDPWTR1R2R3R4_ch_10_10_01.dat\n");
                return(0);
                break;
        }

    
    if(optind==argc){ //raw_contact_data_file(s)
        fprintf(stderr, "Need at least one input file with the raw data (x, y, count).\n");
        exit(EXIT_FAILURE);
    }

/* To apease cc on macOS. */
chr_num=0;

    if (flag_identifier_set) {
		/* Copied from below */
                if(!strcmp(identifier1,identifier2)){
                    if(identifier1[0]>='0' && identifier1[0]<='9') chr_w_chr_num=chr_num=atoi(identifier1);
                    else if (identifier1[0]=='X' || identifier1[0]=='x') chr_w_chr_num=chr_num=LAST_CHR-1;
                    else if (identifier1[0]=='Y' || identifier1[0]=='y') chr_w_chr_num=chr_num=LAST_CHR;
                    strcpy(chr[chr_num].filename,argv[optind]);
                    strcpy(chr_w_filename,argv[optind]);
                    /*if(identifier2[0]>='0' && identifier2[0]<='9') chr_num=atoi(identifier2);
                    else if (identifier2[0]=='X' || identifier2[0]=='x') chr_num=LAST_CHR-1;
                    else if (identifier2[0]=='Y' || identifier2[0]=='y') chr_num=LAST_CHR;*/
                } else {
                    fprintf(stderr,"The filename %s (in position %i) contains contacts between different chromosomes.\n", argv[i], i);
                    exit(EXIT_FAILURE);
                }
	} else { 
      
    //check and "load" input files with counts

    for(i=optind;i<argc;i++){
        if ((fp=fopen(argv[i], "rt"))==NULL) {
            fprintf(stderr, "Cannot open file \"%s\" (%i::%s)\n", argv[i], errno, strerror(errno));
            exit(EXIT_FAILURE);
        }
        //check pattern
        //strlen(argv[i]);
        if ( (strlen(argv[i])<sizeof(identifier1)) && (strlen(argv[i])<sizeof(identifier2)) ) {
            if (sscanf(argv[i], pattern, identifier1, identifier2)!=2) { /* pattern does not match */
                fprintf(stderr,"The filename %s (in position %i) does not match the pattern %s\n", argv[i], i, pattern);
                exit(EXIT_FAILURE);
            } else { /* pattern matches and the identifier is found. */
                fprintf(stdout,"# Info: the file %s matches the pattern, the identifiers are %s and %s\n", argv[i], identifier1, identifier2);
                if(!strcmp(identifier1,identifier2)){
                    if(identifier1[0]>='0' && identifier1[0]<='9') chr_w_chr_num=chr_num=atoi(identifier1);
                    else if (identifier1[0]=='X' || identifier1[0]=='x') chr_w_chr_num=chr_num=LAST_CHR-1;
                    else if (identifier1[0]=='Y' || identifier1[0]=='y') chr_w_chr_num=chr_num=LAST_CHR;
                    strcpy(chr[chr_num].filename,argv[i]);
                    strcpy(chr_w_filename,argv[i]);
                    /*if(identifier2[0]>='0' && identifier2[0]<='9') chr_num=atoi(identifier2);
                    else if (identifier2[0]=='X' || identifier2[0]=='x') chr_num=LAST_CHR-1;
                    else if (identifier2[0]=='Y' || identifier2[0]=='y') chr_num=LAST_CHR;*/
                } else {
                    fprintf(stderr,"The filename %s (in position %i) contains contacts between different chromosomes.\n", argv[i], i);
                    exit(EXIT_FAILURE);
                }
            }
        }
        fclose(fp);
    }

    } /* else part of if (flag_identifier_set) */
    
    if ((flag_identifier_set) && (argc>optind+1)) fprintf(stderr, "Too many input files, only %s will be taken into account.\n", chr_w_filename);
    if(flag_windows && (argc>optind+1)) fprintf(stderr, "Too many input files, only %s will be taken into account.\n", chr_w_filename);

    if(flag_windows) { //generate windows
        fprintf(stderr,"# Info: windows are spaced by delta_w = %lf (bp) and have length Delta_w = %lf (bp) .\n", delta_w, Delta_w);
        //length
        // mm9 reference genome
        //Chromosome_size=[('chr1',197195432),('chr2',181748087),('chr3',159599783),('chr4',155630120),('chr5',152537259),('chr6',149517037),('chr7',152524553),('chr8',131738871),('chr9',124076172),('chr10',129993255),('chr11',121843856),('chr12',121257530),('chr13',120284312),('chr14',125194864),('chr15',103494974),('chr16',98319150),('chr17',95272651),('chr18',90772031),('chr19',61342430),('chrX',166650296),('chrY',15902555)]
        chr_length[1]=197195432;
        chr_length[2]=181748087;
        chr_length[3]=159599783;
        chr_length[4]=155630120;
        chr_length[5]=152537259;
        chr_length[6]=149517037;
        chr_length[7]=152524553;
        chr_length[8]=131738871;
        chr_length[9]=124076172;
        chr_length[10]=129993255;
        chr_length[11]=121843856;
        chr_length[12]=121257530;
        chr_length[13]=120284312;
        chr_length[14]=125194864;
        chr_length[15]=103494974;
        chr_length[16]=98319150;
        chr_length[17]=95272651;
        chr_length[18]=90772031;
        chr_length[19]=61342430;
        chr_length[LAST_CHR-1]=166650296;
        chr_length[LAST_CHR]=15902555;
        
/* pruess 27 May 2021
 * This code is SIGSEGVing when called via
 * ../RUN37_DNAsliding_window_debug/DNAsliding_window -w 5000:5000 -C 1:1 /home/clustor2/ma/g/gunsim/DNAfolding/HICdata_2kb_20181026/RUN36_HeLa/GSM_HeLa_G1/GSM4106789_HeLa_G1_HiRes_controlKRobs_ch_1_1_5kb.dat
 * gdb says it happens in line 533:
 * Program received signal SIGSEGV, Segmentation fault.
0x000055555555731a in main (argc=6, argv=0x7fffffffe6c8) at DNAsliding_window.c:533
533	                                    chr_w_moms_distance(i,j)+=mom*count;
 * I suspect that this is because the genome is longer than expected.
 * The only place chr_w_chr_num is ever used is here. But chr_w_length enters 
 * into all sorts of things, including gplist.
 * 
 * So, I will add a flag to tell the code to find out the length all by itself.
 * 
 * Ok. This needs to be done using -l 0
 *
	Testing: ../RUN37_DNAsliding_window_debug/DNAsliding_window -w 5000:5000 -w 5000:50000 -C 16:16 -l 0 ../RUN04_contact_domains/CD69negDPCTCFKOR1R2/CD69negDPCTCFKOR1R2_ch_16_16_5kb.dat > tmp3.dat
	...
	Info: Chromosome length provided was 0, is 0.000000 will be 98315000.000000
        chr_length[16]=98319150;
	Sounds good to me. I am a bit worried by an off-by-one.
 *
 */


/* You can force the length to be determined by the code by saying -l 0 */

        if (length_provided>=0) chr_w_length=length_provided;
	else if ((chr_w_chr_num<=LAST_CHR) && (chr_w_chr_num>=0)) chr_w_length=chr_length[chr_w_chr_num];

	gplist_init(chr_w_length, Delta_w, chr_w_filename);

	printf("# Info: Chromosome length provided was %lli, is %f will be %f\n", length_provided, chr_w_length, gplist_chr_w_length);
	chr_w_length=gplist_chr_w_length;


        /* Rosalba: I don't think the next line is correct -- below you calculate the max index as floor(x/delta_w); */
        //chr_w_number_windows = ceil((chr_w_length-Delta_w)/delta_w);
        chr_w_number_windows = floor(chr_w_length/delta_w);

	printf("# Info: chr_w_length=%g, delta_w=%g, chr_w_number_windows=%i\n", chr_w_length, delta_w, chr_w_number_windows);

        MALLOC(chr_w_count, chr_w_number_windows);
        MALLOC(chr_w_moms_distance, chr_w_number_windows*NUM_MOMS);

        /* DEPRECATED
	 * The final 1 in gplist_init is for "no coarse graining" -- which means that Delta_w, 
	 * the window size, is identical to the resolution, so that no edges connect the same
	 * two (coarse grained) locations; but in gplist_reorganise, I am now able to cope
	 * with repeated targets, I just merge them there and then. */





        for(i=0;i<chr_w_number_windows;i++){ chr_w_count[i]= 0.; for(j=0;j<NUM_MOMS;j++){ chr_w_moms_distance(i,j)=0.; } }
    } else if(flag_compartments){
        fprintf(stderr,"# Info: reading compartments from provided file.\n");
        chr_num_prev=k=0;
        while (fgets(buffer, sizeof(buffer)-1, fg)!=NULL) {
            if(buffer[0]=='c'){
                if (sscanf(buffer,"%*[^0-9*XYxy]%s %lf %lf %s %*s %s",chrID,&x0,&x1,compID,sign)) {
                    if(chrID[0]>='0' && chrID[0]<='9') chr_num=atoi(chrID);
                    else if(chrID[0]=='X' || chrID[0]=='x') chr_num =LAST_CHR-1;
                    else if(chrID[0]=='Y' || chrID[0]=='y') chr_num =LAST_CHR;
                    if(chr_num!=chr_num_prev){ k=0; }
                    chr[chr_num].comp[k].start=x0;
                    chr[chr_num].comp[k].end=x1;
                    if(strcmp(sign,"+")==0) chr[chr_num].comp[k].group=+1;
                    else chr[chr_num].comp[k].group=-1;
                    chr_num_prev=chr_num;
                    chr[chr_num].num_comp++;
                    k++;
                }
            }
        }
    } else if(flag_ranked_compartments){
        fprintf(stderr,"# Info: reading ranked compartments from provided file.\n");
        chr_num_prev=k=0;
        while (fgets(buffer, sizeof(buffer)-1, fg)!=NULL) {
            if(buffer[0]=='c'){
                if (sscanf(buffer, "%*[^0-9*XYxy]%[^,],%lf,%lf,%d,%*s",chrID,&x0,&x1,&group)) {
                    if(chrID[0]>='0' && chrID[0]<='9') chr_num=atoi(chrID);
                    else if(chrID[0]=='X' || chrID[0]=='x') chr_num =LAST_CHR-1;
                    else if(chrID[0]=='Y' || chrID[0]=='y') chr_num =LAST_CHR;
                    if(chr_num!=chr_num_prev){ k=0; }
                    chr[chr_num].comp[k].start=x0;
                    chr[chr_num].comp[k].end=x1;
                    chr[chr_num].comp[k].group=group;
                    chr_num_prev=chr_num;
                    chr[chr_num].num_comp++;
                    k++;
                }
            }
        }
    } else { fprintf(stderr,"# Info: no provided file with compartments.\n"); }
    

    
    if(flag_windows){
        //this part is to determine for a sliding window of fixed size along the DNA the _total_ _number_ of contacts seen to be made from there
        //also first and second moments of contact distance
        
        printf("# Info: parameters delta_w, Delta_w: %g %g\n", delta_w, Delta_w);
        lines=comments=valid=data_in_ensemble=line_spoiled_by_comment=0LL;
        
        // counting starts here
        for(chr_num=0;(chr_num<=LAST_CHR); chr_num++){
            if ((fp=fopen(chr[chr_num].filename, "rt"))!=NULL){
                printf("# Info reading data from file %s\n", chr[chr_num].filename);
                while (fgets(buffer, sizeof(buffer)-1, fp)!=NULL) {
                    buffer[sizeof(buffer)-1]=0;
                    if ((++lines %100000)==0) printf("# Info %lli lines read at binning.\n", lines);
                    if (buffer[0]=='#') {comments++; continue;}
                    else {
                        if (clean==1) {
                            for (p=buffer+1; *p; p++) { if (*p=='#') break;}
                            if (*p=='#') {
                                line_spoiled_by_comment++;
                                printf("# Spoiled line: [");
                                for (p=buffer; *p; p++) {
                                    if (*p!='\n') fputc(*p, stdout);
                                    else printf("\n# ");
                                }
                                printf("]\n");
                                continue;
                            }
                        }
                        if ( (sscanf(buffer,"%lf %lf %Lf",&x,&y,&count)==3) ? (!isnan(count)) : 0 ) {
                            valid++;
                            data=fabs(x-y);
			    /* Rosalba: Lots of boilerplate code here. */
                            indx_min = ceil((x-Delta_w)/delta_w);
                            indx_max = floor(x/delta_w);
                            indy_min = ceil((y-Delta_w)/delta_w);
                            indy_max = floor(y/delta_w);
			    //printf("Indeces: %i %i %i %i of %g %g\n", indx_min, indx_max, indy_min, indy_max, x, y);

/* BEGINNING OF NEW MATERIAL */
			    gplist_add_edge(x/Delta_w, y/Delta_w, count);
			    if ( ((int)(x/Delta_w)) != ((int)(y/Delta_w)) )
			      gplist_add_edge(y/Delta_w, x/Delta_w, count);

/* END OF NEW MATERIAL */
                            for(i=indx_min;i<=indx_max;i++) {
				if (i<0) {
				/* pruess 27 May 2021: This caused another SIGSEGV. For some reason the first entry in this
				 * HeLa file
				 * ../RUN37_DNAsliding_window_debug/DNAsliding_window -w 5000:5000 -C 17:17 -l 0 /home/clustor2/ma/g/gunsim/DNAfolding/HICdata_2kb_20181026/RUN36_HeLa/GSM_HeLa_G1/GSM4106790_HeLa_G1_HiRes_SA1-siRNAKRobs_ch_17_17_5kb.dat
				 * contains zeros. The ceil((x-Delta_w)/delta_w) and the ceil((y-Delta_w)/delta_w) above
				 * then produce -1. This is irrelevant as far as the production of the bin files is concerned.
				 *
				 * 28 Feb 2022
				 * I have now switched off the use of
				 * chr_w_moms_distance(i,j) completely. When we
				 * run 
				 * ./DNAsliding_window -w 1:1 -C 1:1 -l 0 /home/clustor2/ma/g/gunsim/DNAfolding/kntagiantas/sim_bin/control.dat
				 * the code produced lots of errors and eventually complained about scan failed, when the 
				 * buffer looked like it got mangeled. 
				 *
				 */
				  /* pruess 28 Feb 2022: wrong order of arguments fixed. */
				  printf("# Warning: Not addressing %i of %i in line %i of code line %lli of data (Indeces: %i %i %i %i of %g %g)\n", i*NUM_MOMS, chr_w_number_windows*NUM_MOMS, __LINE__, lines, indx_min, indx_max, indy_min, indy_max, x, y);
				  continue;
				}
				if (i<chr_w_number_windows)
                                  chr_w_count[i]+=count;
                                mom=1;
                                //chr_w_moms_distance(i,0)++;
                                for(j=1;j<NUM_MOMS;j++){
                                    mom*=data;
                                    //chr_w_moms_distance(i,j)+=mom*count;
                                }
                            }
                            for(i=indy_min;i<=indy_max;i++) {
				if (i<0) {
				  /* pruess 28 Feb 2022: wrong order of arguments fixed. */
				  printf("# Warning: Not addressing %i of %i in line %i of code line %lli of data (Indeces: %i %i %i %i of %g %g)\n", i*NUM_MOMS, chr_w_number_windows*NUM_MOMS, __LINE__, lines, indx_min, indx_max, indy_min, indy_max, x, y);
				  continue;
				}
				if (i<chr_w_number_windows)
                                  chr_w_count[i]+=count;
                                mom=1;
                                //chr_w_moms_distance(i,0)++;
                                for(j=1;j<NUM_MOMS;j++){
                                    mom*=data;
                                    //chr_w_moms_distance(i,j)+=mom*count;
                                }
                            }
                        } else {
                            printf("# Scan failed (%i) in line %lli: [", __LINE__, lines);
                            for (p=buffer; *p; p++) {
                                if (*p!='\n') fputc(*p, stdout);
                                else printf("\n# ");
                            }
                            printf("]\n");
                        }
                    }
                }
                fclose(fp);
            }
        }
        
        printf("# Info: Read %lli lines, %lli comments, %lli valid, %lli spoiled; %lli unaccounted.\n",
               lines, comments, valid, line_spoiled_by_comment, lines-valid-comments-line_spoiled_by_comment);
        
        if (comments+valid!=lines) {
            fprintf(stderr, "WARNING: Some invalid lines were not comments.\n");
            fprintf(stdout, "# WARNING: Some invalid lines were not comments.\n");
        }

	//printf("vsiiting\n");
	
	gplist_reorganise();
	{ char command[4096];

	#warning "Should check length of string written."
	for (i=0; i<argc; i++) {
	  if (i) sprintf(command+strlen(command), " \"%s\"", argv[i]);
	  else strcpy(command, argv[0]);
	}
	gplist_prep_header(GPLIST_HDRFLG_SORTED, command);
	}
	printf("# Info: Writing the file (which might take a moment).\n");
	gplist_write_sorted("gplist_data.bin");
	printf("# Info: Done! So long and thanks for all the fish.\n");
	//gplist_print_edges_sorted();
	exit(0);
	{
	int visited;
        for (visited=1; visited<20; visited++) 
	  if (gplist_traversal(visited)) {/*exit(0);*/}
	}
	exit(0);
        
        for(N=0, i=0;i<chr_w_number_windows;i++) N+=chr_w_count[i];
        
        //print results
        for(i=0;i<chr_w_number_windows;i++){
            printf("CHR%d\t%d\t%.16g\t%.16LG\t%.16LG\t%.16LG",chr_w_chr_num,i, (i*delta_w+0.5*Delta_w),chr_w_count[i],chr_w_count[i]/(long double)Delta_w,N/chr_w_number_windows);
            for(j=0;j<NUM_MOMS;j++){ printf("\t%.16LG",chr_w_moms_distance(i,j)); }
            printf("\n");
        }
    } else{
    
    /*.. here starts the code that bins the raw data into a nice density function ..
     * The following couple of environments are not useful for our setting, it will be
     * avoided with the new default parameters (to avoid scanning through massive files).
     */
if ((spit_out_minmax==0) && (fp==stdin) && (min>max)) {
  fprintf(stderr, "You cannot bin a stdin stream without stating min and max.\n");
  exit(EXIT_FAILURE);
}

#define DD fprintf(stderr, "%s::%i\n", __FILE__, __LINE__)

valid=-1LL;
lines=comments=line_spoiled_by_comment=0LL;
if ((spit_out_minmax==1) || (min>max)) {
  valid=0LL;
  lines=0LL;
  while (fgets(buffer, sizeof(buffer)-1, fp)!=NULL) {
    if ((++lines %1000000)==0) printf("# Info %lli lines read at minmax.\n", lines);
    buffer[sizeof(buffer)-1]=0;
    if (buffer[0]=='#') {comments++; continue;}
    else {
      if (clean==1) {
	for (p=buffer+1; *p; p++) { if (*p=='#') break;}
	if (*p=='#') {
	  line_spoiled_by_comment++; 
	  printf("# Spoiled line: [");
	  for (p=buffer; *p; p++) {
	    if (*p!='\n') fputc(*p, stdout);
	    else printf("\n# ");
	  }
	  printf("]\n");
	  continue;
	}
      }
      if (sscanf(buffer,"%Lf %Lf",&data,&count)==2) {
        valid++;
	if (data<=0) printf("# Info data=%Lg in line %lli, data excluded.\n", data, lines);
	else if (min>max) min=max=data;
	else if(data<min) min=data; //here is where min and max are defined
	else if(data>max) max=data;
      } else {
          printf("# Scan failed (%i) in line: [", __LINE__);
	for (p=buffer; *p; p++) {
	  if (*p!='\n') fputc(*p, stdout);
	  else printf("\n# ");
	}
	printf("]\n");
      }
    }
  }
rewind(fp);
}

/* This is in the case that valid lines should have been found, but have not. */
if (valid==0LL) {
  printf("# No valid lines. Exiting.\n");
  exit(0);
} //environments that we want to switch off end here (to avoid scanning through the massive files)

printf("# Info: Read %lli lines, %lli comments, %lli valid, %lli spoiled; %lli unaccounted.\n", 
  lines, comments, valid, line_spoiled_by_comment, lines-valid-comments-line_spoiled_by_comment);
  

if (spit_out_minmax==1) {
  printf("Min: %10.20g\nMax: %10.20g\n-m %10.20g:%10.20g\n", min, max, min, max);
  return(0);
}
  
printf("# Info: parameters nbins, min, max: %i %g %g\n", nbins, min, max);

if ((min<0) || (max<0)) {
  fprintf(stderr, "min<0 or max<0, %g %g\n", min, max);
  exit(EXIT_FAILURE);
}

//log binning base
double b=pow(10,1/((double)nbins));			
smin=pow(b,(int)(log(min)/log(b)))/sqrt(sqrt(b));
int kmax=(int)(log(max/smin)/log(b))+1;		

printf("# Info: parameters b, smin, kmax: %g %g %i\n", b, smin, kmax);


MALLOC(n, kmax*NBINS);
MALLOC(s, kmax+1);

for(i=0;i<kmax;i++){
  s[i]=pow(b,(double)i)*smin;
  for(g=1;g<=NBINS;g++) n(i,g)=0LL;
}
s[kmax]=pow(b,(double)kmax)*smin;

    lines=comments=valid=data_in_ensemble=line_spoiled_by_comment=0LL;

    // counting starts here
    for(chr_num=0;(chr_num<=LAST_CHR); chr_num++){
      if ((fp=fopen(chr[chr_num].filename, "rt"))!=NULL){
        printf("# Info reading data from file %s\n", chr[chr_num].filename);
        while (fgets(buffer, sizeof(buffer)-1, fp)!=NULL) {
          buffer[sizeof(buffer)-1]=0;
          if ((++lines %100000)==0) printf("# Info %lli lines read at binning.\n", lines);
          if (buffer[0]=='#') {comments++; continue;}
          else {
            if (clean==1) {
              for (p=buffer+1; *p; p++) { if (*p=='#') break;}
              if (*p=='#') {
                line_spoiled_by_comment++;
                printf("# Spoiled line: [");
                for (p=buffer; *p; p++) {
                  if (*p!='\n') fputc(*p, stdout);
                  else printf("\n# ");
                }
                printf("]\n");
                continue;
              }
            }
            //if (sscanf(buffer,"%lf %lf %Lf",&x,&y,&count)==3 && (!isnan(count))) {
            if ( (sscanf(buffer,"%lf %lf %Lf",&x,&y,&count)==3) ? (!isnan(count)) : 0 ) {
              valid++;
                if(flag_compartments || flag_ranked_compartments){
                    find_compartment(x, chr[chr_num], &indx, &groupx, labelx);
                    find_compartment(y, chr[chr_num], &indy, &groupy, labely);
                }
              data=fabs(x-y);
              // if the points satisfy the constraints, data is accounted
                if( (flag_compartments || flag_ranked_compartments) ? ((data>=min) && (data<=max) && (groupx>0 || groupy>0)) : ((data>=min) && (data<=max)) ){
                data_in_ensemble++;
                n((int)(log(data/smin)/log(b)),groupx)+=count;   // counting
                n((int)(log(data/smin)/log(b)),groupy)+=count;   // counting
              } else data_out_of_minmax=1;
            } else {
              printf("# Scan failed (%i) in line: [", __LINE__);
              for (p=buffer; *p; p++) {
                if (*p!='\n') fputc(*p, stdout);
                else printf("\n# ");
              }
            printf("]\n");
          }
        }
      }
      fclose(fp);
    }
  }

  printf("# Info: Read %lli lines, %lli comments, %lli valid, %lli spoiled; %lli unaccounted.\n",
  lines, comments, valid, line_spoiled_by_comment, lines-valid-comments-line_spoiled_by_comment);
  printf("# Info: Out of %lli valid lines, %lli were included in the ensemble considered.\n",
         valid, data_in_ensemble);
    
  if (comments+valid!=lines) {
    fprintf(stderr, "WARNING: Some invalid lines were not comments.\n");
    fprintf(stdout, "# WARNING: Some invalid lines were not comments.\n");
  }


//print results for each group of bins
for(g=1;g<NBINS;g++){

  // get total
  for(N=0, i=0;i<kmax;i++) N+=n(i,g);
  if (N==0) N=-1;
  fprintf(stdout, "\n# Info: total number of counts is %Lf\n",N);
	
// CAUTION!! 	modify s[0] and s[kmax] _only_if_ min and max where supplied by hand
//		and if there was data outside the stated min:max. In this case, we
//		are actually filtering data, and so the min
//		and max are not "statistically natural". so we need special (smaller)
//		last bins.
// 		do it as welf if the -t flag is passed

    //RGM I don't get this comment (nor the two enviroments below)
    
if( force_last_bin_correction){

	fprintf(stderr,	"WARNING: option -t caused first and last bin to be\n"
			"         corrected (they are smaller). This assumes that\n"
			"         the data has been truncated or filtered.\n");
	fprintf(stdout,	"# WARNING: option -t caused first and last bin to be\n"
			"#          corrected (they are smaller). This assumes that\n"
			"#	    the data has been truncated or filtered.\n");
		
	for(i=0;i<=kmax;i++) if(n(i,g)!=0){
		s[i]=min;
		break;
	}
	for(i=kmax-1;i>=0;i--) if(n(i,g)!=0){
		s[i+1]=max;
		break;
	}

}


if( data_out_of_minmax ){

	fprintf(stdout,	"# WARNING: some data was outside the range (possibly provided\n"
			"#          via -m) as a result, first and last bin were corrected\n"
			"#          corrected (they became smaller)\n");
		
	fprintf(stderr,	"WARNING: some data was outside the range (possibly provided\n"
			"         via -m) as a result, first and last bin were\n"
			"         corrected (they became smaller)\n");
		
	for(i=0;i<=kmax;i++) if(n(i,g)!=0){
		s[i]=min;
		break;
	}
	for(i=kmax-1;i>=0;i--) if(n(i,g)!=0){
		s[i+1]=max;
		break;
	}
}

// normalizei and print
for(i=0; i<kmax-1; i++){
  if(resolution!=0){
  /* resolution really is the resolution of the data. resolution is 1 if it's discrete, it's 0 if it's continuous.
   * The next line checks how big the bins are _given_ a certain resolution.
   */
	De=resolution*(ceil(s[i+1]/resolution)-ceil(s[i]/resolution));
  }
  else De=s[i+1]-s[i];	
  
  /* I think there should be a 
   * if (De==0.) continue;
   * here. Used to be just
   * dens=(double)n[i]/(De*(double)N);
   */

  if (De==0.) {
    freq=dens=-1.;
    if (n(i,g)>0) {
      printf("# %Lf entries ignored as De==0.\n", n(i,g));
    }
  } else { dens=n(i,g)/(De*N); freq=n(i,g)/(De); }

  if(resolution!=0){
	data=resolution*sqrt(ceil(s[i]/resolution)*ceil(s[i+1]/resolution-1));
  }
  else data=sqrt(s[i+1]*s[i]);

//		fprintf(stderr,"[%f,%f) <-- %lld, De=%.4G \t %.16G\t%.16G\n",s[i],s[i+1],n[i],De,data,dens);
/* This is a bit dangerous: If we are under the resolution limit, then we may have
 * data in the bin, but we are not spitting it out. */
  if(De!=0) printf("GROUP%d\t%.16LG\t%.16LG\t%.16LG\t%.16LG\n",g,data,dens,freq,N);
}		

// this is only last bin!!! its different: [a,b] instead of [a,b)	
for(i=kmax-1;i<kmax;i++){
  if(resolution!=0) De=resolution*(floor(s[i+1]/resolution)+1.-ceil(s[i]/resolution));
  else De=s[i+1]-s[i];

  dens=n(i,g)/(De*N);
  freq=n(i,g)/(De);

  if(resolution!=0) data=resolution*sqrt(ceil(s[i]/resolution)*ceil(s[i+1]/resolution));
  else data=sqrt(s[i+1]*s[i]);

  if(De!=0) printf("GROUP%d\t%.16LG\t%.16LG\t%.16LG\t%.16LG\n",g,data,dens,freq,N);
}
}
}
        
return 0;
}



