/*
 * $Header: /Users/garciamillan/Documents/uni/ImperialCL/DNAfolding/HiC_2kb/contactdomains/results34_generalised_code/generalised_code/DNA_stencil.c
 * COMPILE: cc -O3 -Wall -o JetsSten DNA_stencil.c gplist_util2.c -lm -lgsl -lgslcblas
 * EXECUTE: ./JetsSten > test.dat
 * EXECUTE: ./JetsSten -r 500000 -R 5000 -n 17 > test.dat
 * date: Thu 14 Oct 2021 11:34:21 BST
 *
 * To generate the binary files, run
 * ./DNAsliding_window -w 5000:5000 -C 17:17 ../RUN04_contact_domains/CD69negDPCTCFKOR1R2/CD69negDPCTCFKOR1R2_ch_17_17_5kb.dat
 * /bin/mv -f gplist_data.bin CD69negDPCTCFKOR1R2_ch_17_17_5kb_coarsened5kb.bin
 * gzip CD69negDPCTCFKOR1R2_ch_17_17_5kb_coarsened5kb.bin
 */



#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <time.h>
#include "gplist_util2.h"
#include "data_structure.h"

int read_depth(struct data_struct *dt);
int read_hic(const char *filename, struct gplist_header_struct **data, int chr);
int read_expected(const char *basenameEXP, double **data, int chr);
int protractor_grid(double rad, const double resolution, int n, struct point_stencil **point);
int stencil_condition(struct condition_struct *condition, int n);
int match_parents_daughters(struct condition_struct *condition);
int Gaussian_filter(struct condition_struct *condition, int num_sectors);
int generate_points_stencil(struct data_struct *dt, int k, int n, struct point_stencil *point,int *num_points);
int stencil_background(struct data_struct *dt, int k, struct point_stencil *point, int num_points, int n, const double resolution, double *background);
int background_condition(struct condition_struct *condition, struct data_struct *dt, int i, int n);
int subtraction_condition(struct condition_struct *condition, int n);
int Gaussian_filter_subtraction(struct condition_struct *condition, int num_sectors);
int Gaussian_filter_background(struct condition_struct *condition, int num_sectors);
int jet_depth(struct condition_struct *condition, int num_sectors, double threshold, double resolution);

double angle(double x, double y){
    double theta;
    if(y>0) {
     if(x>0) theta=atan(y/x);
     else if(x==0) theta=M_PI_2;
     else theta=M_PI+atan(y/x);}
    else theta= M_PI+atan(y/x);
    return theta;
}

int sector(double x, double y, int n){
    if ((x==0) && (y==0)) return -1;
    double theta;
    theta=angle(x,y);
    if (floor((1.25-theta/M_PI))==1) return n-1;
    return floor(n*(1.25-theta/M_PI)); //this way s increases as x increases
}

int stencil_grid(double theta, int n, struct point_stencil **point);
double boundary_bottom (double x, double theta);
double boundary_top (double x, double theta);
double boundary_left (double x, double theta);
double boundary_right (double x, double theta);
double slab(double x, double y, double theta);
int stencil_curve(struct data_struct *dt, int k, struct point_stencil *point, int num_points, int n, const double resolution);


int main(int argc, char *argv[]){
    setlinebuf(stdout);
    double rad;
    int n,num_jet_candidates, num_points;
    char buffer[8192];
    double x0, x1;
    double theta;
    int k, chr_num, chr_num_read, chr_num_prev;
    int i;
    FILE *fg=stdin;
    int ch;
    int replicate;
    int num_replicates=(int)(sizeof(data)/sizeof(*data));
    int num_conditions=(int)(sizeof(condition)/sizeof(*condition));
    struct point_stencil *point=NULL;

    n=NUM_SLABS;

    printf("# $Header$\n");
    printf("# Command:");
    for (i=0; i<argc; i++) printf(" %s", argv[i]);
    printf("\n");
    { time_t tm;
      tm=time(NULL);
      printf("# Info: %s", ctime(&tm));
    }
    
    { char hostname[128];
        gethostname(hostname, sizeof(hostname)-1);
        hostname[sizeof(hostname)-1]=(char)0;
        printf("# Info: Hostname: %s\n", hostname);
    }
    
    { char cwd[1024];
        cwd[0]=0;
        if(getcwd(cwd, sizeof(cwd)-1)!=NULL){
        cwd[sizeof(cwd)-1]=(char)0;
        printf("# Info: Directory: %s\n", cwd);}
    }
    
    printf("# Info: PID: %i\n", (int)getpid());
    

    
    while ((ch = getopt(argc, argv, "r:n:h")) != -1)
            switch (ch) {
                case 'r':
                    if ((rad=atof(optarg))<0) {
                        fprintf(stderr, "radius r must be non-negative.\n");
                        exit(EXIT_FAILURE);
                    }
                    break;
                case 'n':
                    if ((n=atoi(optarg))<=0) {
                        fprintf(stderr, "number of sectors  must be positive.\n");
                        exit(EXIT_FAILURE);
                    }
                    break;
                default:
                    fprintf(stderr, "Flag %c unknown.\n", ch);
                case 'h':
                    fprintf(stderr, "Usage: ./dna -f filename_compartments -p pattern -M -m min:max -b nbins -r resolution file_contact1 file_contact2 file_contact3 ...\n");
                    fprintf(stderr, "-f filename_compartments specified input file (with compartments). By default the program reads from stdin.\n"
                            "-r radius of the circle sector (default 5e5bp)\n"
                            "-R resolution of the lattice (default 5e3bp)\n"
                            "-n number of sectors we consider in half a circle (default 9)\n"
                            "-c chromosome number (by default, n=1).\n"
                            " ./dna -c 2 -f CD69negDPCTCFKOR1R2_eigenvector_KR100kb.bed -r 500000 -R 5000 -n 9\n");
                    return(0);
                    break;
            }
    
    {//two reads to identify jet candidates from the .bed file
    if ((fg=fopen(jet_analysis_file, "rt"))==NULL) {
        fprintf(stderr, "Can't open file \"%s\" (%i::%s)\n", jet_analysis_file, errno, strerror(errno));
        exit(EXIT_FAILURE);
    }

    printf("# Info: Reading compartments from jet_analysis_file %s\n",jet_analysis_file);
        
        chr_num_read=chr_num_prev=0;
        num_jet_candidates=0;
        
    //the first read is to find how many chromosomes and jet candidates there are
        int maxIndexOfChromosomes=0;
    printf("# Info: First read.\n");
    while (fgets(buffer, sizeof(buffer)-1, fg)!=NULL) {
        if(buffer[0]=='P'){
            if (sscanf(buffer,"%*s %*s %d", &chr_num_read)) {
                    maxIndexOfChromosomes=MAX(maxIndexOfChromosomes, chr_num_read);
                    num_jet_candidates++;
            }
        }
    }

        printf("# Info: Total number of jet candidates %d\n"
               "# Info: Total number of chromosomes %d\n",num_jet_candidates,maxIndexOfChromosomes);
        if(num_jet_candidates>=MAX_INDEX_OF_POSITIONS){
            printf("# Warning: Curtailing num_jet_candidates = %i to %i\n", num_jet_candidates, MAX_INDEX_OF_POSITIONS);
            num_jet_candidates=MAX_INDEX_OF_POSITIONS;
        }
            
        for(replicate=0;replicate<num_replicates;replicate++) {
            data[replicate].num_positions=num_jet_candidates;
            MALLOC(data[replicate].position, data[replicate].num_positions);
//#warning "Never free data[replicate].position, because it has been shifted."
        }
    rewind(fg);

        char label[100];
        //this is hard-coded again. may lead to issues
        double angle0, angle1, angle2, strength0, strength1, strength2;
    k=0;
        
        for (i=0;i<num_conditions;i++){
        condition[i].num_positions=condition[i].data[0]->num_positions;
        MALLOC(condition[i].position,condition[i].num_positions);
        }
        

    printf("# Info: Second read.\n");
    replicate=0;
    while (fgets(buffer, sizeof(buffer)-1, fg)!=NULL && k<data[0].num_positions) {
        if(buffer[0]=='P'){
            if (sscanf(buffer,"%s %d %d %lf %lf %lf %lf %lf %lf %lf %lf", label,&k,&chr_num_read,&x0,&x1,&angle0,&angle1,&angle2,&strength0,&strength1,&strength2) && k<data[0].num_positions) {
                
                for(replicate=0;replicate<num_replicates;replicate++) {
                    data[replicate].position[k].chromosome=chr_num_read;
                    data[replicate].position[k].start=x0*1000.;
                    data[replicate].position[k].end=x1*1000.;
                    strcpy(data[replicate].position[k].label, label);
                }
                condition[0].position[k].angle=angle0;
                condition[1].position[k].angle=angle1;
                condition[2].position[k].angle=angle2;
                condition[0].position[k].strength=strength0;
                condition[1].position[k].strength=strength1;
                condition[2].position[k].strength=strength2;
                for (i=0;i<num_conditions;i++) {
                    condition[i].position[k].flag_has_stencil=0;
                    strcpy(condition[i].position[k].label, label);
                }

            }
        }
    }
    fclose(fg);

    printf("# Info: Finished reading jet candidates.\n");

    }

    { // read list of consensus jets and flick the flag_has_stencil
        if ((fg=fopen(jet_consensus_file, "rt"))==NULL) {
            fprintf(stderr, "Can't open file \"%s\" (%i::%s)\n", jet_analysis_file, errno, strerror(errno));
            exit(EXIT_FAILURE);
        }
        
        printf("# Info: Consensus jets:\n");

        char positionID[100];
        
        while (fgets(buffer, sizeof(buffer)-1, fg)!=NULL) {
            if(buffer[0]=='P'){
                if (sscanf(buffer,"%*[^0-9*XYxy]%s", positionID)) {
                    if(positionID[0]>='0' && positionID[0]<='9') {
                        k=atoi(positionID)-1;
                        for (i=0;i<num_conditions;i++) condition[i].position[k].flag_has_stencil=1;
                        printf("%s\n",condition[0].position[k].label);
                    }
                }
            }
        }
    }
    
    {//load bamfiles (to find white lines)
        printf("# Info: Reading BAM_depth files to find white lines.\n");
        for(replicate=0;replicate<num_replicates;replicate++) read_depth(&(data[replicate]));
        printf("# Info: Finished reading BAM_depth files.\n");
    }

    { //generate stencils
        MALLOC(point,STENCIL_MAXPOINTS);
        for (i=0;i<num_conditions;i++) match_parents_daughters(&(condition[i]));
    printf("# Info: Generating stencil curves.\n");
    chr_num_prev=0;
    int replicate_prev=-1;
        printf("# Info: Generating grid of points in rectangle (stencil).\n");
    for(replicate=0;replicate<num_replicates;replicate++) {
        printf("# Info: Working on replicate %s\n", data[replicate].ident);
        for(k=0; k<data[replicate].num_positions; k++) if(data[replicate].parent->position[k].flag_has_stencil){
            printf("# Info: Working on candidate %s, in chromosome %d, replicate %s\n",
                   data[replicate].position[k].label,data[replicate].position[k].chromosome,data[replicate].ident);
	printf("# Woking %s : %d [ %lf , %lf ]\n",data[replicate].position[k].label,data[replicate].position[k].chromosome,data[replicate].position[k].start, data[replicate].position[k].end);
            
            chr_num=data[replicate].position[k].chromosome;
            if(chr_num!=chr_num_prev || replicate_prev!=replicate) {
                gplist_free(data[replicate].hic);
                free(data[replicate].expected);
                
                //load HiC data with KR normalisation (observed)
                read_hic(data[replicate].basenameHIC, &(data[replicate].hic), chr_num);
            
                //load HiC expected values
                read_expected(data[replicate].basenameEXP, &(data[replicate].expected), chr_num);
                chr_num_prev=chr_num;
                replicate_prev=replicate;
            }
                     
            if(!isnan(data[replicate].parent->position[k].angle)) generate_points_stencil(&(data[replicate]),k,n,point,&num_points);
            else data[replicate].parent->position[k].flag_has_stencil=0;
            
            if(data[replicate].parent->position[k].flag_has_stencil) stencil_curve(&(data[replicate]),k,point,num_points,n,resolution);
            
            if(flag_background && replicate>=num_replicates-2) if(data[replicate].parent->position[k].flag_has_stencil) {
                MALLOC(data[replicate].position[k].background_stencil, num_conditions);
                for(i=0;i<num_conditions;i++) {
                    MALLOC(data[replicate].position[k].background_stencil[i].stencil, n);
                    theta=condition[i].position[k].angle * M_PI/180.;
                    num_points = stencil_grid(theta,n,&point);
                    printf("# Info: The background stencil grid of %s %s has %d points.\n",condition[i].position[k].label,condition[i].ident,num_points);
                    
                    if((num_points)<n) {
                        fprintf(stderr,"# Error: The background stencil grid of %s %s has not enough points:  %d < %d .\n",condition[i].position[k].label,condition[i].ident,num_points,n);
                        exit(EXIT_FAILURE);
                    }

                    MALLOC(data[replicate].position[k].background_stencil[i].stencil,n);
                    stencil_background(&(data[replicate]),k,point,num_points,n,resolution,data[replicate].position[k].background_stencil[i].stencil);
            }
        }
    }
    }}
    
    
    {//group replicates by condition, calculate mean and std stencil
        printf("# Info: Calculating stencil curve mean and variance for each condition.\n");
        for (i=0;i<num_conditions;i++) stencil_condition(&(condition[i]),n);
    }
    
    {//Gaussian filter on stencil
        for (i=0;i<num_conditions;i++) Gaussian_filter(&(condition[i]),n);
    }
    
    {//generate subtraction stencils (WT or CTCFko minus background given by Dko)
        //ignore whitelines for the moment
        if(flag_background) {
            //calculate background mean and std for each condition
            for (i=0;i<num_conditions;i++) background_condition(&(condition[i]),data,i,n);

            //subtract background
            for (i=0;i<num_conditions;i++) subtraction_condition(&(condition[i]), n);
            
            //Gaussian filter on subtraction stencil
            for (i=0;i<num_conditions;i++) Gaussian_filter_subtraction(&(condition[i]),n);
            
            //Gaussian filter on background
            for (i=0;i<num_conditions;i++) Gaussian_filter_background(&(condition[i]),n);
        }
    }
    
    
    {//Find jet depth (= "distance" on HiC map where the stencil curve minus background falls below a certain threshold)
        //and find projection of "jet ends" on chromatin
        double threshold=SIGNAL_THRESHOLD;
        for (i=0;i<num_conditions;i++) jet_depth(&(condition[i]),n,threshold,resolution);
    }

    double label_slab[n];
    int j;
    
    for(i=0;i<n;i++) label_slab[i]= (i+0.5) * ((double)STENCIL_LENGTH) * (resolution/1000.) / ((double)NUM_SLABS);

        printf("# Info: Printing output.\n");
    
    {//print stencil curves
        printf("STENCIL #label  k  chromosome  start  end  index  position_stencil");
        for (i=0;i<num_conditions;i++) printf("  %s_mean  %s_std  %s_smooth",condition[i].ident, condition[i].ident,condition[i].ident);
        if(flag_background) for (i=0;i<num_conditions;i++) printf("  %s_subtraction_mean  %s_subtraction_std  %s_subtraction_smooth",condition[i].ident, condition[i].ident,condition[i].ident);
        if(flag_background) for (i=0;i<num_conditions;i++) printf("  %s_background_mean  %s_background_std  %s_background_smooth",condition[i].ident, condition[i].ident,condition[i].ident);
        printf("\n");

        for(k=0; k<data[0].num_positions; k++) for(j=0;j<n;j++) if(data[0].parent->position[k].flag_has_stencil) {
            printf("STENCIL %s  %d  %d  %.16G  %.16G  %d  %.16G",
                   data[0].position[k].label, k, data[0].position[k].chromosome, data[0].position[k].start, data[0].position[k].end,  j, label_slab[j]);
            for (i=0;i<num_conditions;i++) {
                if(condition[i].position[k].flag_has_stencil) printf("  %.16G  %.16G  %.16G", condition[i].position[k].stencil_mean[j], condition[i].position[k].stencil_std[j], condition[i].position[k].stencil_smooth[j]);
                else printf("  nan  nan  nan");
            }
            if(flag_background) {
                for (i=0;i<num_conditions;i++) {
                if(condition[i].position[k].flag_has_stencil) printf("  %.16G  %.16G  %.16G",
                                                                     condition[i].position[k].stencil_subtraction_mean[j], condition[i].position[k].stencil_subtraction_std[j], condition[i].position[k].stencil_subtraction_smooth[j]);
                else printf("  nan  nan  nan");
            }
                for (i=0;i<num_conditions;i++) {
                    if(condition[i].position[k].flag_has_stencil) printf("  %.16G  %.16G  %.16G", condition[i].position[k].stencil_background_mean[j], condition[i].position[k].stencil_background_std[j], condition[i].position[k].stencil_background_smooth[j]);
                    else printf("  nan  nan  nan");
                }
            }
            printf("\n");
        }
    }
    
    {//print jet depth and projection
        for(k=0; k<data[0].num_positions; k++) if(data[0].parent->position[k].flag_has_stencil) {
            printf("#PROJECTION  %s  %d  %.16G  %.16G",
                   data[0].position[k].label, data[0].position[k].chromosome, (data[0].position[k].start)/1000., (data[0].position[k].end)/1000.);
            for (i=0;i<num_conditions;i++) if(condition[i].position[k].flag_has_stencil) printf("  %.16G  %.16G  %.16G  %.16G", round((condition[i].position[k].upstream_reach)/5)*5, round((condition[i].position[k].upstream_reach_std)/5)*5, round((condition[i].position[k].downstream_reach)/5)*5, round((condition[i].position[k].downstream_reach_std)/5)*5);
            else printf("  nan  nan  nan  nan");
            printf("\n");
        }
        
        double midpoint;
        printf("D1D2 #label  k  chromosome  start  end");
        for (i=0;i<num_conditions;i++) printf("  %s_d1  %s_d1_std  %s_d2  %s_d2_std",condition[i].ident, condition[i].ident,condition[i].ident,condition[i].ident);
        printf("\n");

        for(k=0; k<data[0].num_positions; k++) if(data[0].parent->position[k].flag_has_stencil) {
            printf("#D1D2  %s  %d  %d  %.16G  %.16G",
                   data[0].position[k].label, k+1, data[0].position[k].chromosome, (data[0].position[k].start)/1000., (data[0].position[k].end)/1000.);
            midpoint = 0.5*((data[0].position[k].start)/1000. + (data[0].position[k].end)/1000.);
            for (i=0;i<num_conditions;i++) if(condition[i].position[k].flag_has_stencil) printf("  %.16G  %.16G  %.16G  %.16G", - round((condition[i].position[k].upstream_reach)/5)*5 + midpoint, round((condition[i].position[k].upstream_reach_std)/5)*5, round((condition[i].position[k].downstream_reach)/5)*5 - midpoint, round((condition[i].position[k].downstream_reach_std)/5)*5);
            else printf("  nan  nan  nan  nan");
            printf("\n");
        }

    }
    
    {//print GNUPLOT commands
    char filePath[PATH_MAX],fileName[PATH_MAX];
    if (readlink("/proc/self/fd/1", filePath, sizeof(filePath)-1)!=-1) strncpy(fileName, filePath + 73, PATH_MAX - 73);
                
        printf("#GNUPLOT00 FILE = '%s'\n"
            "#GNUPLOT00 a=-10\n"
            "#GNUPLOT00 b=150\n"
            "#GNUPLOT00 set yrange [a:b]\n"
            "#GNUPLOT00 set xlabel 'Distance along stencil (kb)'\n"
            "#GNUPLOT00 set ylabel 'Strength'\n"
            "#GNUPLOT00 set key out horizontal bottom\n"
            "#GNUPLOT00 set zeroaxis\n"
            "#GNUPLOT00\n"
               "#GNUPLOT00\n", fileName);
        
        for(k=0; k<data[0].num_positions; k++) if(data[0].parent->position[k].flag_has_stencil) {
            printf("#GNUPLOT02_%d set term post eps color solid\n"
            "#GNUPLOT02_%d set out '%s_HeLa_stencil.eps'\n",
                   k+1,k+1,data[0].position[k].label);
            printf("#GNUPLOT%d set title 'HeLa %s: chr%d [ %.16G , %.16G ]'\n ", k+1,data[0].position[k].label,data[0].position[k].chromosome,data[0].position[k].start,data[0].position[k].end);
            printf("#GNUPLOT%d plot ",k+1);
            for (i=0;i<num_conditions;i++){
                printf(" '<grep %s[^0-9] '.FILE.'' u %s lc '%s' t \"%s\","
                       " '<grep %s[^0-9] '.FILE.'' u %s lc '%s' t \"\"", data[0].position[k].label, condition[i].plotting2,condition[i].colour,condition[i].ident,data[0].position[k].label, condition[i].plotting3,condition[i].colour);
                if(i<num_conditions-1) printf(",");
            }
            printf("\n");

        }

    }
    
    if(flag_background){//print GNUPLOT commands for subtraction and division stencils
    char filePath[PATH_MAX],fileName[PATH_MAX];
    if (readlink("/proc/self/fd/1", filePath, sizeof(filePath)-1)!=-1) strncpy(fileName, filePath + 73, PATH_MAX - 73);
                
        printf("#GNUPLOT_SUB00 FILE = '%s'\n"
            "#GNUPLOT_SUB00 a=-10\n"
            "#GNUPLOT_SUB00 b=100\n"
            "#GNUPLOT_SUB00 set yrange [a:b]\n"
            "#GNUPLOT_SUB00 set xlabel 'Distance along stencil (kb)'\n"
            "#GNUPLOT_SUB00 set ylabel 'Strength - Dko'\n"
            "#GNUPLOT_SUB00 set key out horizontal bottom\n"
            "#GNUPLOT_SUB00 set zeroaxis\n"
            "#GNUPLOT_SUB00\n"
               "#GNUPLOT_SUB00\n", fileName);
        
        for(k=0; k<data[0].num_positions; k++) if(data[0].parent->position[k].flag_has_stencil) {
            printf("#GNUPLOT_SUB%d set term post eps color solid\n"
            "#GNUPLOT_SUB%d set out '%s_subtraction_gnu.eps'\n",
                   k+1,k+1,data[0].position[k].label);
            printf("#GNUPLOT_SUB%d set title '%s: chr %d [ %.16G , %.16G ]'\n",
                   k+1,data[0].position[k].label, data[0].position[k].chromosome,data[0].position[k].start, data[0].position[k].end);
            printf("#GNUPLOT_SUB%d plot '<grep %s[^0-9] '.FILE.'' u 7:17:18 w e pt 4 lc 'black' t \"Wild-type\","
                   " '<grep %s[^0-9] '.FILE.'' u 8:20 w l lw 3 lc 'gray' t \"\","
                   " '<grep %s[^0-9] '.FILE.'' u 8:21:22 w e pt 6 lc 3 t \"CTCF ko\","
                   " '<grep %s[^0-9] '.FILE.'' u 8:23 w l lw 3 lc 3 t \"\","
                   " '<grep %s[^0-9] '.FILE.'' u 8:24:25 w e pt 8 lc 'pink' t \"RAD21 CTCF ko\","
                   " '<grep %s[^0-9] '.FILE.'' u 8:26 w l lw 3 lc 'pink' t \"\"\n",
                   k+1,data[0].position[k].label,data[0].position[k].label,data[0].position[k].label,data[0].position[k].label,data[0].position[k].label,data[0].position[k].label);
        }

        
        printf("#GNUPLOT_BACK00 FILE = '%s'\n"
            "#GNUPLOT_BACK00 a=-10\n"
            "#GNUPLOT_BACK00 b=100\n"
            "#GNUPLOT_BACK00 set yrange [a:b]\n"
            "#GNUPLOT_BACK00 set xlabel 'Distance along stencil (kb)'\n"
            "#GNUPLOT_BACK00 set ylabel 'background (Dko)'\n"
            "#GNUPLOT_BACK00 set key out horizontal bottom\n"
            "#GNUPLOT_BACK00 set zeroaxis\n"
            "#GNUPLOT_BACK00\n"
               "#GNUPLOT_BACK00\n", fileName);
        
        for(k=0; k<data[0].num_positions; k++) if(data[0].parent->position[k].flag_has_stencil) {
            printf("#GNUPLOT_BACK%d set term post eps color solid\n"
            "#GNUPLOT_BACK%d set out '%s_background_gnu.eps'\n",
                   k+1,k+1,data[0].position[k].label);
            printf("#GNUPLOT_BACK%d set title '%s: chr %d [ %.16G , %.16G ]'\n",
                   k+1,data[0].position[k].label, data[0].position[k].chromosome,data[0].position[k].start, data[0].position[k].end);
            printf("#GNUPLOT_BACK%d plot '<grep %s[^0-9] '.FILE.'' u 8:27:28 w e pt 4 lc 'black' t \"Wild-type\","
                   " '<grep %s[^0-9] '.FILE.'' u 8:27 w l lw 3 lc 'gray' t \"\","
                   " '<grep %s[^0-9] '.FILE.'' u 8:28:29 w e pt 6 lc 3 t \"CTCF ko\","
                   " '<grep %s[^0-9] '.FILE.'' u 8:30 w l lw 3 lc 3 t \"\","
                   " '<grep %s[^0-9] '.FILE.'' u 8:31:32 w e pt 8 lc 'pink' t \"RAD21 CTCF ko\","
                   " '<grep %s[^0-9] '.FILE.'' u 8:33 w l lw 3 lc 'pink' t \"\"\n",
                   k+1,data[0].position[k].label,data[0].position[k].label,data[0].position[k].label,data[0].position[k].label,data[0].position[k].label,data[0].position[k].label);
        }

    }
    
    
printf("# Info: c'est fini.\n");

return 0;
}



int read_depth(struct data_struct *dt)
{
    FILE *fg;
    char *p;
    char buffer[8192];
    int chr;
    int locus,depth,index;
    char chrID[MAX_CHROMOSOME_NAME];

    /* Determine number of chromosomes and each chromosome length */
    
    for (chr=0; chr<=MAX_INDEX_OF_CHROMOSOMES; chr++) {
        dt->BAM[chr].maxIndex=-1;
        dt->BAM[chr].chromosomeName[0]=(char)0;
    }
    if ((fg=fopen(dt->filenameBAM, "rt"))==NULL) {
        fprintf(stderr, "Can't open file \"%s\" (%i::%s)\n", dt->filenameBAM, errno, strerror(errno));
        exit(EXIT_FAILURE);
    } else printf("#Info: reading file \"%s\"\n", dt->filenameBAM);

    dt->maxIndexOfChromosomes=-1;
    while (fgets(buffer, sizeof(buffer)-1, fg)!=NULL) {
            if(buffer[0]=='c'){
                chr=0;
                if (sscanf(buffer,"chr%s %d %d",chrID,&locus,&depth)) {
                    if(chrID[0]>='0' && chrID[0]<='9') chr=strtol(chrID, &p, 10);
                    else {chr=-1; continue;}
                    /* Ignoring X,Y chromosomes. */
                    //else if(chrID[0]=='X' || chrID[0]=='x') chr =LAST_CHR-1;
                    //else if(chrID[0]=='Y' || chrID[0]=='y') chr =LAST_CHR;
                    //strings like chr17_random will be ignored
                    if (*p!=0)
                        if (!isspace(*p)) {chr=-1; continue;}
                    
                    index=(int)(locus/resolution-CHROMOSOME_OFFSET);
                    /* At this stage, we keep track only of max of chr and max of index. */
                    dt->maxIndexOfChromosomes=MAX(dt->maxIndexOfChromosomes, chr);
                    if ((chr>=0) && (chr<=MAX_INDEX_OF_CHROMOSOMES)) {
                        if (dt->BAM[chr].maxIndex==-1) strcpy(dt->BAM[chr].chromosomeName, chrID);
                        dt->BAM[chr].maxIndex=MAX(dt->BAM[chr].maxIndex,index);
                    }
                }
            }
        }
    
    if (dt->maxIndexOfChromosomes>=MAX_INDEX_OF_CHROMOSOMES) {
        printf("# Warning: Curtailing dt->maxIndexOfChromosomes=%i to %i\n", dt->maxIndexOfChromosomes, MAX_INDEX_OF_CHROMOSOMES);
        dt->maxIndexOfChromosomes=MAX_INDEX_OF_CHROMOSOMES;
    }
    
    
    /* Now we know the maximum index of the chromosomes and the maximum indeces of those ... */
    for (chr=0; chr<=dt->maxIndexOfChromosomes; chr++) {
        printf("# Info: Chromosome %i is called [%s] has maximum index %i\n", chr, dt->BAM[chr].chromosomeName, dt->BAM[chr].maxIndex);
        if (dt->BAM[chr].maxIndex>0) {
            MALLOC(dt->BAM[chr].depth, dt->BAM[chr].maxIndex+1);
            for(index=0;index<=dt->BAM[chr].maxIndex;index++) dt->BAM[chr].depth[index]=0.;
        } else {
            dt->BAM[chr].depth=NULL;
            strcat(dt->BAM[chr].chromosomeName, "--EMPTY--");
        }
    }
    
    rewind(fg);
    
#warning "Boilerplate code."
    
    while (fgets(buffer, sizeof(buffer)-1, fg)!=NULL) {
            if(buffer[0]=='c'){
                chr=0;
                if (sscanf(buffer,"chr%s %d %d",chrID,&locus,&depth)) {
                    if(chrID[0]>='0' && chrID[0]<='9') chr=strtol(chrID, &p, 10);
                    else {chr=-1; continue;}
                    /* Ignoring X,Y chromosomes. */
                    //else if(chrID[0]=='X' || chrID[0]=='x') chr =LAST_CHR-1;
                    //else if(chrID[0]=='Y' || chrID[0]=='y') chr =LAST_CHR;
                    //strings like chr17_random will be ignored
                    if (*p!=0)
                        if (!isspace(*p)) {chr=-1; continue;}
                    
                    index=(int)(locus/resolution-CHROMOSOME_OFFSET);
                    /* At this stage, we keep track only of max of chr and max of index. */
                    if ((chr>=0) && (chr<=dt->maxIndexOfChromosomes)) {
                        if ((index>=0) && (index<=dt->BAM[chr].maxIndex)) dt->BAM[chr].depth[index]+=depth;
                    }
                }
            }
        }

    /* End of copied code */
    fclose(fg);

    return(0);
}


int read_hic(const char *basenameHIC, struct gplist_header_struct **data, int chr){
    if ((*data)!=NULL) {
        printf("# Warning: Pointer to data is not NULL. Does this need to be free(3)d?\n");
    }

    char filename[PATH_MAX];
    sprintf(filename, "%s%d_%d_5kb.bin.gz", basenameHIC,chr,chr);
    
    printf("# Info: Reading HiC from file %s\n", filename);
    if (((*data)=gplist_read(filename))==NULL) {
        fprintf(stderr, "# Error: Failed to read %s. (%i::%s)\n", filename, errno, strerror(errno));
        exit(EXIT_FAILURE);
    }

    return 0;
}

int read_expected(const char *basenameEXP, double **data, int chr){
    FILE *fexpected;
    char buffer[8192];
    int i;
    double expected;
    char filename[PATH_MAX];

    if ((*data)!=NULL) {
        printf("# Warning: Pointer to data is not NULL. Does this need to be free(3)d?\n");
    }
    
    MALLOC((*data),LENGHT_EXP);
    sprintf(filename, "%s%d_%d_5kb.dat", basenameEXP,chr,chr);

    
    if ((fexpected=fopen(filename, "rt"))==NULL) {
        fprintf(stderr, "Cannot open file \"%s\" (%i::%s)\n", filename, errno, strerror(errno));
        exit(EXIT_FAILURE);
    }
    
    printf("# Info: Reading expected from file %s, chromosome %d.\n", filename, chr);
    i=0;
    while(fgets(buffer, sizeof(buffer)-1, fexpected)!=NULL) if(sscanf(buffer,"%lf",&expected)) {
        (*data)[i++]=expected;
            if ((i>=LENGHT_EXP) || (i<0)) {
                    fprintf(stderr, "i=%i out of bounds [0,%i]\n", i, LENGHT_EXP);
                    exit(EXIT_FAILURE);
            }
        }
    fclose(fexpected);
    
    return 0;
}


double boundary_bottom (double x, double theta) {
    double m = tan(0.75*M_PI-theta);
    return -x/m-0.5;
}
double boundary_top (double x, double theta) {
    double m = tan(0.75*M_PI-theta);
    return -x/m -STENCIL_LENGTH*sqrt(1.+m*m)/m+0.5;
}
double boundary_left (double x, double theta) {
    double m = tan(0.75*M_PI-theta);
    return m*x +0.5*STENCIL_WIDTH*m*sqrt(1.+1./(m*m))-0.5;
}
double boundary_right (double x, double theta) {
    double m = tan(0.75*M_PI-theta);
    return m*x -0.5*STENCIL_WIDTH*m*sqrt(1.+1./(m*m))+0.5;
}
double slab(double x, double y, double theta){
    double m = tan(0.75*M_PI-theta);
    return floor(-m*((double)NUM_SLABS)*(x/m+y)/(((double) STENCIL_LENGTH)*sqrt(1.+m*m)));
}


int stencil_grid(double theta, int n, struct point_stencil **point){
    //expects to see theta in radians
    int i=0;
    double x,y;
    
    for(y=STENCIL_LENGTH; y>0 && i<=STENCIL_MAXPOINTS; y--){ //y>-0.5*STENCIL_WIDTH
        for(x=-STENCIL_LENGTH; x<0 && i<=STENCIL_MAXPOINTS; x++){ //x<0.5*STENCIL_WIDTH
            if(y<boundary_top(x,theta) && y>boundary_left(x,theta) && y<boundary_right(x,theta)) {
                
                //y>boundary_bottom(x,theta) &&
                if(i<STENCIL_MAXPOINTS){
                    (*point)[i].x=x;
                    (*point)[i].y=y;
                    (*point)[i].slab=slab(x,y,theta);
                }
                if ((*point)[i].slab>=0 && (*point)[i].slab<n) i++;
            }
        }
    }
    if(i==STENCIL_MAXPOINTS) printf("# Warn: stencil reached STENCIL_MAXPOINTS %d points (max is %d).\n",i,STENCIL_MAXPOINTS);
    printf("# Info: stencil has %d points (max is %d).\n",i,STENCIL_MAXPOINTS);

    return i;
}

int generate_points_stencil(struct data_struct *dt, int k, int n, struct point_stencil *point,int *num_points){
    double theta = dt->parent->position[k].angle * M_PI/180.;

    *num_points = stencil_grid(theta,n,&point);
    dt->parent->position[k].stencil_num_points = *num_points;
    if((*num_points)<n) {
        fprintf(stderr,"# Error: The stencil grid of %s %s has not enough points:  %d < %d .\n",dt->parent->position[k].label,dt->parent->ident,*num_points,n);
        exit(EXIT_FAILURE);
    }
    printf("# Info: The stencil grid of %s %s has %d points.\n",dt->parent->position[k].label,dt->parent->ident,*num_points);
    return 0;
}

int stencil_curve(struct data_struct *dt, int k, struct point_stencil *point, int num_points, int n, const double resolution){
    //based on the function panoramic_curve
    int i,sector_number,chr_num,epos;
    double count;
    int x_bp,y_bp;
    
    if (dt->hic==NULL) { fprintf(stderr, "Error: dt->hic==NULL unexpected in line %i of %s\n", __LINE__, __FILE__); exit(EXIT_FAILURE); }
    if (dt->expected==NULL) { fprintf(stderr, "Error: dt->expected==NULL unexpected in line %i of %s\n", __LINE__, __FILE__); exit(EXIT_FAILURE); }
    if (point==NULL) { fprintf(stderr, "Pointer 'point' is NULL\n"); exit(EXIT_FAILURE); }
    if ((k>= dt->num_positions) || (k<0)) { fprintf(stderr, "k=%i out of bounds [0,%i]\n", k, dt->num_positions-1); exit(EXIT_FAILURE); }
    
    MALLOC(dt->position[k].stencil, n);
    MALLOC(dt->position[k].whitelines, n);
    
    for(sector_number=0;sector_number<n;sector_number++) dt->position[k].stencil[sector_number]=dt->position[k].whitelines[sector_number]=0.; //initialisation

    for(i=0;i<num_points;i++) {
        x_bp=floor(point[i].x+(dt->position[k].start)/resolution);
        y_bp=floor(point[i].y+(dt->position[k].end)/resolution);
        epos=abs(x_bp-y_bp);

        sector_number=(int)point[i].slab;
        if ((sector_number>=n) || (sector_number<0)) { fprintf(stderr, "sector_number = %i out of bounds [0,%i]\n", sector_number, n); exit(EXIT_FAILURE); }
        if ((epos<0) || (epos>=LENGHT_EXP)) { fprintf(stderr, "Error: epos=%i out of bounds 0..%i\n", epos, LENGHT_EXP);  exit(EXIT_FAILURE); }
        
        count = gplist_fetch_count(dt->hic,x_bp,y_bp) - dt->expected[epos];
        dt->position[k].stencil[sector_number]+= count;

        chr_num=dt->position[k].chromosome;
        
        if (((int)(x_bp-CHROMOSOME_OFFSET) > dt->BAM[chr_num].maxIndex) || ((int)(x_bp-CHROMOSOME_OFFSET)<0)) dt->position[k].whitelines[sector_number]++;
        else if(dt->BAM[chr_num].depth[(int)(x_bp-CHROMOSOME_OFFSET)]==0) dt->position[k].whitelines[sector_number]++;
        
        if (((int)(y_bp-CHROMOSOME_OFFSET) > dt->BAM[chr_num].maxIndex) || ((int)(y_bp-CHROMOSOME_OFFSET)<0)) dt->position[k].whitelines[sector_number]++;
        else if(dt->BAM[chr_num].depth[(int)(y_bp-CHROMOSOME_OFFSET)]==0) dt->position[k].whitelines[sector_number]++;
        }
    
    //normalise by total counts
    for(sector_number=0;sector_number<n;sector_number++) dt->position[k].stencil[sector_number] *= PREF/dt->total_counts;
        
    return 0;
}


int stencil_background(struct data_struct *dt, int k, struct point_stencil *point, int num_points, int n, const double resolution, double *background){
    //based on the function stencil_curve
    int i,sector_number,epos;
    double count;
    int x_bp,y_bp;
    
    if (dt->hic==NULL) { fprintf(stderr, "Error: dt->hic==NULL unexpected in line %i of %s\n", __LINE__, __FILE__); exit(EXIT_FAILURE); }
    if (dt->expected==NULL) { fprintf(stderr, "Error: dt->expected==NULL unexpected in line %i of %s\n", __LINE__, __FILE__); exit(EXIT_FAILURE); }
    if (point==NULL) { fprintf(stderr, "Pointer 'point' is NULL\n"); exit(EXIT_FAILURE); }
    if ((k>= dt->num_positions) || (k<0)) { fprintf(stderr, "k=%i out of bounds [0,%i]\n", k, dt->num_positions-1); exit(EXIT_FAILURE); }
        
    for(sector_number=0;sector_number<n;sector_number++) background[sector_number]=0.; //initialisation

    for(i=0;i<num_points;i++) {
        x_bp=floor(point[i].x+(dt->position[k].start)/resolution);
        y_bp=floor(point[i].y+(dt->position[k].end)/resolution);
        epos=abs(x_bp-y_bp);

        sector_number=(int)point[i].slab;
        if ((sector_number>=n) || (sector_number<0)) { fprintf(stderr, "sector_number = %i out of bounds [0,%i]\n", sector_number, n); exit(EXIT_FAILURE); }
        if ((epos<0) || (epos>=LENGHT_EXP)) { fprintf(stderr, "Error: epos=%i out of bounds 0..%i\n", epos, LENGHT_EXP);  exit(EXIT_FAILURE); }
        
        count = gplist_fetch_count(dt->hic,x_bp,y_bp) - dt->expected[epos];
        background[sector_number]+= count;

        }
    
    //normalise by total counts
    for(sector_number=0;sector_number<n;sector_number++) background[sector_number] *= PREF/dt->total_counts;
        
    return 0;
}

int match_parents_daughters(struct condition_struct *condition){
    int replicate;
    for(replicate=0; replicate<condition->num_replicates; replicate++){
        condition->data[replicate]->parent= condition;
        printf("#Info %s is daughter of %s\n",condition->data[replicate]->ident,condition->data[replicate]->parent->ident);
    }
    
    return 0;
}

    
int stencil_condition(struct condition_struct *condition, int n){
    int replicate,k,j,num_points;
    double m0;
    double tolerance_whitelines;
    
    for(k=0; k<condition->num_positions; k++) if(condition->position[k].flag_has_stencil){
        MALLOC(condition->position[k].stencil_mean,n);
        MALLOC(condition->position[k].stencil_std,n);
        num_points=condition->position[k].stencil_num_points;
        tolerance_whitelines=TOLWHITELINES*num_points/((double)n);

        for(j=0;j<n;j++) {
            condition->position[k].stencil_mean[j] = condition->position[k].stencil_std[j] = 0.; //initialise
            
            //calculate mean
            m0=0.;
            for(replicate=0; replicate<condition->num_replicates; replicate++){
                if(condition->data[replicate]->position[k].whitelines[j] < tolerance_whitelines){
                    condition->position[k].stencil_mean[j]+= condition->data[replicate]->position[k].stencil[j];
                    m0++;
                    }
            }
            if(m0<1) condition->position[k].stencil_mean[j]= nan("1");
            else condition->position[k].stencil_mean[j]*= 1./m0;

            //calculate std
            if(m0>1){
                for(replicate=0;replicate<condition->num_replicates;replicate++){
                    if(condition->data[replicate]->position[k].whitelines[j] < tolerance_whitelines) condition->position[k].stencil_std[j]+= pow(condition->data[replicate]->position[k].stencil[j] - condition->position[k].stencil_mean[j],2);
                }
                condition->position[k].stencil_std[j] = sqrt(condition->position[k].stencil_std[j]/(m0*(m0-1.)));
            }
    }   }
    
    return 0;
}

int background_condition(struct condition_struct *condition, struct data_struct *dt, int i, int n){
    int replicate,k,j;//num_points;
    double m0;
    int num_replicates=(int)(sizeof(data)/sizeof(*data));
    //I will ignore white lines here
    //double tolerance_whitelines;
    
    for(k=0; k<condition->num_positions; k++) if(condition->position[k].flag_has_stencil){
        MALLOC(condition->position[k].stencil_background_mean,n);
        MALLOC(condition->position[k].stencil_background_std,n);

        for(j=0;j<n;j++) {
            condition->position[k].stencil_background_mean[j] = condition->position[k].stencil_background_std[j] = 0.; //initialise
            
            //calculate mean
            m0=0.;
            for(replicate=6; replicate<num_replicates; replicate++){
                //if(condition->data[replicate]->position[k].whitelines[j] < tolerance_whitelines){
                    condition->position[k].stencil_background_mean[j]+= data[replicate].position[k].background_stencil[i].stencil[j];
                    m0++;
            }
            if(m0<1) condition->position[k].stencil_background_mean[j]= nan("1");
            else condition->position[k].stencil_background_mean[j]*= 1./m0;

            //calculate std
            if(m0>1){
                for(replicate=6;replicate<num_replicates;replicate++){
                        condition->position[k].stencil_background_std[j]+= pow(data[replicate].position[k].background_stencil[i].stencil[j] - condition->position[k].stencil_background_mean[j],2);
                }
                condition->position[k].stencil_background_std[j] = sqrt(condition->position[k].stencil_background_std[j]/(m0*(m0-1.)));
            }
        }}
    return 0;
}

int subtraction_condition(struct condition_struct *condition, int n){
    int k,j;
    //the subtraction stencil effectively inherits the whitelines from stencil_mean
    
    for(k=0; k<condition->num_positions; k++) if(condition->position[k].flag_has_stencil){
        MALLOC(condition->position[k].stencil_subtraction_mean,n);
        MALLOC(condition->position[k].stencil_subtraction_std,n);

        for(j=0;j<n;j++) {
            condition->position[k].stencil_subtraction_mean[j] = condition->position[k].stencil_mean[j] - condition->position[k].stencil_background_mean[j];
            condition->position[k].stencil_subtraction_std[j] = sqrt(condition->position[k].stencil_std[j]*condition->position[k].stencil_std[j] + condition->position[k].stencil_background_std[j]*condition->position[k].stencil_background_std[j]);
        }}
    
    return 0;
}

double Gaussian(double x, double h){
    return exp(-0.5*x*x/(h*h))/(h*M_SQRT2*sqrt(M_PI));
}

int Gaussian_filter(struct condition_struct *condition, int num_sectors) {
    int i,j,k;
    double bandwidth=BANDWIDTH;

    for(k=0;k<condition->num_positions;k++) if(condition->position[k].flag_has_stencil){
        MALLOC(condition->position[k].stencil_smooth,num_sectors);
        for(i=0;i<num_sectors;i++) {
            //smoothing through Gaussian kernel
            condition->position[k].stencil_smooth[i]=0.; //initialisation
            for(j=0;j<num_sectors;j++) {
                if(!isnan(condition->position[k].stencil_mean[j])) condition->position[k].stencil_smooth[i] += condition->position[k].stencil_mean[j]*Gaussian((double)(i-j),bandwidth);
            }
        }
    }
    return 0;
}

int Gaussian_filter_subtraction(struct condition_struct *condition, int num_sectors) {
    int i,j,k;
    double bandwidth=BANDWIDTH;

    for(k=0;k<condition->num_positions;k++) if(condition->position[k].flag_has_stencil){
        MALLOC(condition->position[k].stencil_subtraction_smooth,num_sectors);
        for(i=0;i<num_sectors;i++) {
            //smoothing through Gaussian kernel
            condition->position[k].stencil_subtraction_smooth[i]=0.; //initialisation
            for(j=0;j<num_sectors;j++) {
                if(!isnan(condition->position[k].stencil_subtraction_mean[j])) condition->position[k].stencil_subtraction_smooth[i] += condition->position[k].stencil_subtraction_mean[j]*Gaussian((double)(i-j),bandwidth);
            }
        }
    }
    return 0;
}

int Gaussian_filter_background(struct condition_struct *condition, int num_sectors) {
    int i,j,k;
    double bandwidth=BANDWIDTH;

    for(k=0;k<condition->num_positions;k++) if(condition->position[k].flag_has_stencil){
        MALLOC(condition->position[k].stencil_background_smooth,num_sectors);
        for(i=0;i<num_sectors;i++) {
            //smoothing through Gaussian kernel
            condition->position[k].stencil_background_smooth[i]=0.; //initialisation
            for(j=0;j<num_sectors;j++) {
                if(!isnan(condition->position[k].stencil_background_mean[j])) condition->position[k].stencil_background_smooth[i] += condition->position[k].stencil_background_mean[j]*Gaussian((double)(i-j),bandwidth);
            }
        }
    }
    return 0;
}


int project_to_chromatin(struct position_struct *position, int i_thrs, double resolution) {
    double theta,m,x,y,x_std,y_std,dummy;
    
    theta = position->angle * M_PI/180.;
    
    m = tan(0.75*M_PI-theta);
    x = 1./sqrt(1.+m*m) *(-i_thrs * STENCIL_LENGTH / NUM_SLABS );
    x_std = 1./sqrt(1.+m*m) *(1. * STENCIL_LENGTH / NUM_SLABS );
    
    y = m*x;
    y_std = m*x_std;
    
    //make sure that x is smaller than y
    if(x>y) {
        dummy=x; x=y; y=dummy;
        dummy=x_std; x_std=y_std; y_std=dummy;
    }

    position->upstream_reach= (x*resolution+position->start)/1000.;
    position->downstream_reach= (y*resolution+position->end)/1000.;
    position->upstream_reach_std= x_std*resolution/1000.;
    position->downstream_reach_std= -y_std*resolution/1000.;
    

    return 0;
}


int jet_depth(struct condition_struct *condition, int num_sectors, double threshold, double resolution) {
    int i,i_max,i_thrs,k;

    for(k=0;k<condition->num_positions;k++) if(condition->position[k].flag_has_stencil){
        //find first local maximum
        i_max=0;
        i=1;
        while (condition->position[k].stencil_subtraction_smooth[i] > condition->position[k].stencil_subtraction_smooth[i_max] && i<num_sectors){
            i_max=i;
            i++;
        }
        //find last instance before the signal decays below the threshold
        i=i_max;
        i_thrs=0;
        while (condition->position[k].stencil_subtraction_smooth[i] > threshold && i<num_sectors){
            i_thrs=i;
            i++;
        }
        //project i_thrs to find upstream and dowsnstream reach
        condition->position[k].start= condition->data[0]->position[k].start;
        condition->position[k].end= condition->data[0]->position[k].end;
        project_to_chromatin(&(condition->position[k]),i_thrs,resolution);

        printf("#i_thrs  Pos%d  %s  %d  %.16G  %.16G\n",k+1, condition->ident, i_thrs, (double) (i_thrs+0.5) * STENCIL_LENGTH * 5. / (double) NUM_SLABS, condition->position[k].stencil_subtraction_smooth[i_thrs]);
    }
    return 0;
}
