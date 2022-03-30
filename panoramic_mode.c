/*
 * $Header: /Users/garciamillan/Documents/uni/ImperialCL/DNAfolding/HiC_2kb/contactdomains/results31_mode_curvature/mode_curvature.c$
 * COMPILE: cc -O3 -Wall -o max panoramic_mode.c rgm_Gfit2.c -lm -lgsl -lgslcblas
 * EXECUTE: ./max > test.dat
 * date: Sat 20 Mar 2021 13:58:42 GMT
 *
 * To generate the panoramic file, run
 * ./Jets
 * from /Users/garciamillan/Documents/uni/ImperialCL/DNAfolding/HiC_2kb/contactdomains/results31_mode_curvature/DNA_jets_replicates2.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <time.h>
//#include "gplist_util2.h"
#include "rgm_Gfit2.h"
#include "data_structure.h"

#ifndef PATH_MAX
#define PATH_MAX (1024)
#endif

#ifndef NUM_CANDIDATES
#define NUM_CANDIDATES (277)
#endif

#define MALLOC(a,n) if ((a=malloc(sizeof(*a)*(n)))==NULL) { fprintf(stderr, "Not enough memory for %s, requested %i bytes, %i items of size %i. %i::%s\n", #a, (int)(sizeof(*a)*n), n, (int)sizeof(*a), errno, strerror(errno)); exit(EXIT_FAILURE); } else { printf("# Info malloc(3)ed %i bytes (%i items of %i bytes) for %s.\n", (int)(sizeof(*a)*(n)), n, (int)sizeof(*a), #a); }



#define point(i,j) point[i*4+j]
#define LAST_CHR (24)
#define UPPER_BOUND_CHR (19)
#define TOLWHITELINES (0.1)
#define MAX_LENGTH_BAM (38841)
#define LENGHT_EXP (40000)
//#define LBSTRENGHT_STRONG (3.)//(4.)//(3.)
#define LBSTRENGHT_MEDIUM (1.)//(3.)//(1.)
#define LB_STD (0.5)
#define LENGTH_STAMP (2000.)
#define UPPER_BOUND_ANGLE_FIT (0.33333*M_PI)
#define UPPER_BOUND_ANGLE_JET (0.25*M_PI)
#define JET_OPENING (45.)
#define BANDWIDTH (1.)
#define xFITTING_RANGE (8)//(20.)
#define FITTING_RANGE (6)//(20.)
#define SHORT_FIT_RANGE_IND (4)

#define DD fprintf(stderr, "%s::%i\n", __FILE__, __LINE__)

//Total contacts in each replicate
#define TOTALCONTACTS_WTR1 (244438962.)
#define TOTALCONTACTS_WTR2 (223193977.)
#define TOTALCONTACTS_WTR3 (417830413.)
#define TOTALCONTACTS_WTR4 (427257794.)
#define TOTALCONTACTS_DKOR1 (158131390.)
#define TOTALCONTACTS_DKOR2 (175281770.)
#define TOTALCONTACTS_CTCFKOR1 (424715008.)
#define TOTALCONTACTS_CTCFKOR2 (325929677.)

#define PREF (1.e8)
#define NORMALISATION_WTR1 (PREF/TOTALCONTACTS_WTR1)
#define NORMALISATION_WTR2 (PREF/TOTALCONTACTS_WTR2)
#define NORMALISATION_WTR3 (PREF/TOTALCONTACTS_WTR3)
#define NORMALISATION_WTR4 (PREF/TOTALCONTACTS_WTR4)
#define NORMALISATION_CTCFKOR1 (PREF/TOTALCONTACTS_CTCFKOR1)
#define NORMALISATION_CTCFKOR2 (PREF/TOTALCONTACTS_CTCFKOR2)
#define NORMALISATION_DKOR1 (PREF/TOTALCONTACTS_DKOR1)
#define NORMALISATION_DKOR2 (PREF/TOTALCONTACTS_DKOR2)

#define PREFSTRENGTH (1.92e-6)


int Gaussian_filter(struct condition_struct *condition, int num_sectors);
int find_modes(struct condition_struct *condition, int num_sectors);
int fit_modes(struct condition_struct *condition, int num_sectors, double *label_angle);
int strength(struct condition_struct *condition, int num_sectors);
int integrate_panoramic(struct condition_struct *condition,int num_sectors);

int flag_read_from_bed_file=1;

struct compartment {
    int chromosome;
    double start, end;
    double *panoramicWTmean;
    double *panoramicCTCFKOmean;
    double *panoramicDKOmean;
    double *panoramicWTstd;
    double *panoramicCTCFKOstd;
    double *panoramicDKOstd;
    double *panoramicWTsmooth;
    double *panoramicCTCFKOsmooth;
    double *panoramicDKOsmooth;
    double *curvature_WT;
    double *curvature_CTCFKO;
    double *curvature_DKO;
    int modeWT_ind;
    int modeCTCFKO_ind;
    int modeDKO_ind;
    double jetstrengthWT;
    double jetstrengthCTCFKO;
    double jetstrengthDKO;
    double *offset;
    double normWT;
    double normCTCFKO;
    double normDKO;
    double est_A_WT; // these are calculated by fitting a Gaussian
    double est_sigmasq_WT;
    double std_A_WT;
    double std_sigmasq_WT;
    double est_A_CTCFKO;
    double est_sigmasq_CTCFKO;
    double std_A_CTCFKO;
    double std_sigmasq_CTCFKO;
    double est_A_DKO;
    double est_sigmasq_DKO;
    double std_A_DKO;
    double std_sigmasq_DKO;
};

void print_jet_classificationGNU(double strength){
    if(isnan(strength)) printf(" \\\\textcolor{black}{NO}");
    else if(strength>=LBSTRENGHT_STRONG) printf(" \\\\textcolor{Green}{STRONG}");
    else if(strength>=LBSTRENGHT_MEDIUM) printf(" \\\\textcolor{YellowOrange}{MEDIUM}");
    else if(strength<LBSTRENGHT_MEDIUM) printf(" \\\\textcolor{red}{WEAK}");
}

void print_jet_classificationSIMPLE(double strength){
    if(isnan(strength)) printf(", NO");
    else if(strength>=LBSTRENGHT_STRONG) printf(", STRONG");
    else if(strength>=LBSTRENGHT_MEDIUM) printf(", MEDIUM");
    else if(strength<LBSTRENGHT_MEDIUM) printf(", WEAK");
}

void curvature(int i, double *f, double *second_derivative, double dtheta){
    double a0,a1,a2,a3;
    a0=-49./36.;
    a1=3./4.;
    a2=-3./40.;
    a3=1./180.;
    second_derivative[i] = a0*f[i]+a1*(f[i-1]+f[i+1])+a2*(f[i-2]+f[i+2])+a3*(f[i-3]+f[i+3]);
    second_derivative[i] *= 2./(dtheta*dtheta);
}

double Gaussian(double x, double h){
    return exp(-0.5*x*x/(h*h))/(h*M_SQRT2*sqrt(M_PI));
}




struct compartment *comp;
double *label_angle;
double *xdata_fit;
double *ydata_fit;
double *weights;

int main(int argc, char *argv[]){
    setlinebuf(stdout);
    int num_positions;// = NUM_CANDIDATES;
    //struct compartment comp[NUM_CANDIDATES];
    int i,j,k,chr,start,end,num_sectors, i1,rad;
    double WTmean, WTstd, CTCFKOmean, CTCFKOstd, DKOmean, DKOstd;
    int ch;
    char buffer[8192];
    FILE *fg=stdin;
    //int num_replicates=(int)(sizeof(data)/sizeof(*data));
    int num_conditions=(int)(sizeof(condition)/sizeof(*condition));


    while ((ch = getopt(argc, argv, "h")) != -1)
            switch (ch) {
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
        
    //char filename[PATH_MAX];
    //strcpy(filename, "results/CD69negDPKR_jets_170977_PANORAMIC.dat");

    {//read input file
        
        if ((fg=fopen(panoramic_file, "rt"))==NULL) {
            fprintf(stderr, "Can't open file \"%s\" (%i::%s)\n", panoramic_file, errno, strerror(errno));
            exit(EXIT_FAILURE);
        }

        printf("# Info: Reading panoramics from %s\n",panoramic_file);

        printf("# Info: First read.\n"); //find number of jet candidates and sectors in panoramics

        while (fgets(buffer, sizeof(buffer)-1, fg)!=NULL) {
            if(buffer[0]=='P') sscanf(buffer,"%*s Pos%d %*s %*s %*s %d", &num_positions,&num_sectors);
        }
        printf("# Info: Total number of jet candidates %d\n",num_positions);
        printf("# Info: Total number of sectors in panoramic %d\n",++num_sectors);
        
        if(num_sectors!=NUM_SECTORS) { printf("Error: number of sectors %d != NUM_SECTORS %d\n",num_sectors,NUM_SECTORS); exit(EXIT_FAILURE); }
        
        for (i=0;i<num_conditions;i++) {
            condition[i].num_positions=num_positions;
            MALLOC(condition[i].position,condition[i].num_positions);
            for (k=0; k<condition[i].num_positions; k++){
                MALLOC(condition[i].position[k].panoramic_mean,num_sectors);
                MALLOC(condition[i].position[k].panoramic_std,num_sectors);
            }
        }
        
        printf("# Info: Second read.\n"); //read data
        rewind(fg);

        int index;
        double angle;
        MALLOC(label_angle,num_sectors);

        while (fgets(buffer, sizeof(buffer)-1, fg)!=NULL) {
            if(buffer[0]=='P') if(sscanf(buffer,"%*s Pos%d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf", &k,&chr,&start,&end,&index,&angle,&WTmean,&WTstd,&CTCFKOmean,&CTCFKOstd,&DKOmean,&DKOstd)){
                //I'd like to improve this (generalise it to any number of conditions)
                //printf("%s",buffer);
                if(k>=1 && k<=num_positions && index>=0 && index<num_sectors) {
                    label_angle[index]=angle;
                    for (i=0;i<num_conditions;i++){
                        condition[i].position[k-1].chromosome=chr;
                        condition[i].position[k-1].start=start;
                        condition[i].position[k-1].end=end;
                        sprintf(condition[i].position[k-1].label, "Pos%d", k);
                    }
                    condition[0].position[k-1].panoramic_mean[index]=WTmean;
                    condition[0].position[k-1].panoramic_std[index]=WTstd;
                    condition[1].position[k-1].panoramic_mean[index]=CTCFKOmean;
                    condition[1].position[k-1].panoramic_std[index]=CTCFKOstd;
                    condition[2].position[k-1].panoramic_mean[index]=DKOmean;
                    condition[2].position[k-1].panoramic_std[index]=DKOstd;
                }
            }
        }
    }

    {//Gaussian filter on panoramics
        for (i=0;i<num_conditions;i++) Gaussian_filter(&(condition[i]),num_sectors);
    }
     
    {//find modes
        for (i=0;i<num_conditions;i++) find_modes(&(condition[i]),num_sectors);
    }
    
    {//calculate offsets
        for (i=0;i<num_conditions;i++) for(k=0;k<condition[i].num_positions;k++){
            condition[i].position[k].offset=0.;
            for(j=0;j<num_sectors;j++) condition[i].position[k].offset = MIN(condition[i].position[k].offset,condition[i].position[k].panoramic_mean[j]);
        }
    }
    
    /*{//fit Gaussian around mode
        for (i=0;i<num_conditions;i++) fit_modes(&(condition[i]),num_sectors,label_angle);
    }*/

    {//jet strength
        for (i=0;i<num_conditions;i++) strength(&(condition[i]),num_sectors);
    }
    
    /*{//panoramic integral
        for (i=0;i<num_conditions;i++) integrate_panoramic(&(condition[i]),num_sectors);
    }*/
    
    {//print ANALYSIS_JETS
        printf("#ANALYSIS_JETS #label, k, chromosome,  start,  end");
        for (i=0;i<num_conditions;i++) printf(",  mu_%s",condition[i].ident);
        for (i=0;i<num_conditions;i++) printf(",  strength_%s",condition[i].ident);
        for (i=0;i<num_conditions;i++) printf(",  jet_%s",condition[i].ident);
        printf(",  radius protractor\n");
        for(k=0;k<num_positions;k++){
            printf("#ANALYSIS_JETS %s, %d, %d, %.16G, %.16G",
                   condition[0].position[k].label,k,condition[0].position[k].chromosome,(condition[0].position[k].start)/1000.,(condition[0].position[k].end)/1000.);
            for (i=0;i<num_conditions;i++) { if(condition[i].position[k].mode>0) printf(", %.16G",label_angle[condition[i].position[k].mode]); else printf(", nan"); }
            for (i=0;i<num_conditions;i++) printf(", %.16G",condition[i].position[k].strength);
            //for (i=0;i<num_conditions;i++) printf(", %.16G",condition[i].position[k].panoramic_smooth[condition[i].position[k].mode]);
            //for (i=0;i<num_conditions;i++) printf(", %.16G",condition[i].position[k].integral_panoramic);
            for (i=0;i<num_conditions;i++) print_jet_classificationSIMPLE(condition[i].position[k].strength);
            printf(",  %.16G\n",RADIUS_PROTRACTOR);
        }
    }
    
    
    
    double pan;
    int ind;
    i1 = ceil(0.5*((-JET_OPENING/90.+1.)*(double)num_sectors-1));
    rad=RADIUS_PROTRACTOR;
    {//print GNUPLOT commands
        //char filePath[PATH_MAX];,fileName[PATH_MAX];
    //if (readlink("/proc/self/fd/1", filePath, sizeof(filePath)-1)!=-1) strncpy(fileName, filePath + 73, PATH_MAX - 73);

        printf("#GNUPLOT00 FILE = '%s'\n"
            "#GNUPLOT00 set xrange [-90:90]\n"
            "#GNUPLOT00 a=-200\n"
            "#GNUPLOT00 b=300\n"
            "#GNUPLOT00 set yrange [a:b]\n"
            "#GNUPLOT00 set arrow from 45.,a to 45.,b nohead dt '.'\n"
            "#GNUPLOT00 set arrow from -45.,a to -45.,b nohead dt '.'\n"
            "#GNUPLOT00 set xlabel 'Angle (degrees)'\n"
            "#GNUPLOT00 set ylabel 'Strength'\n"
            "#GNUPLOT00 set key right bmargin\n"
            "#GNUPLOT00 set zeroaxis\n"
            "#GNUPLOT00\n"
               "#GNUPLOT00\n", panoramic_file);

        int flag_print;
        for(k=0; k<num_positions; k++) {
            flag_print=0;
            for (i=0;i<num_conditions;i++) if(!isnan(condition[i].position[k].est_A*condition[i].position[k].norm)) flag_print++;
            if (flag_print>0) {
                printf("#GNUPLOT01_%d set term 'svg'\n"
                "#GNUPLOT01_%d set out '%s_gnu_%dMB.svg'\n"
                "#GNUPLOT02_%d set term post eps color solid\n"
                "#GNUPLOT02_%d set out '%s_gnu_%dMB.eps'\n",
                       k+1,k+1,condition[0].position[k].label,rad/1000000,k+1,k+1,condition[0].position[k].label,rad/1000000);
                
                /*printf("#GNUPLOT%d plot '<grep %s[^0-9] '.FILE.'' u 7:8:9 w e pt 4 lc 'black' t \"Wild-type\","
                       " '<grep %s[^0-9] '.FILE.'' u 7:10:11 w e pt 6 lc 'blue' t \"CTCF ko\","
                       " '<grep %s[^0-9] '.FILE.'' u 7:12:13 w e pt 8 lc 'red' t \"RAD21 CTCF ko\"\n",
                       k+1,condition[0].position[k].label,condition[0].position[k].label,condition[0].position[k].label);*/
                
                printf("#GNUPLOT%d plot",k+1);
                for (i=0;i<num_conditions;i++) {
                    printf(" '<grep %s[^0-9] '.FILE.'' u %s lc '%s' t \"%s\"",condition[0].position[k].label,condition[i].plotting,condition[i].colour,condition[i].ident);
                    if(i<num_conditions-1) printf(",");
                }
                printf("\n#GNUPLOT%d\n",k+1);
            /*for (i=0;i<num_conditions;i++) {
                ind=condition[i].position[k].mode;
                pan=condition[i].position[k].panoramic_mean[ind];
                if(ind>=i1 && !isnan(pan)) printf("#GNUPLOT%d set arrow %d nohead from %.1f,0 to %.1f,%.1f dt '-' lc '%s'\n",k,i+3,label_angle[ind],label_angle[ind],pan,condition[i].colour);
                else printf("#GNUPLOT%d unset arrow %d\n",i+3);
            }
            printf("#GNUPLOT%d plot");
            for (i=0;i<num_conditions;i++) {
                printf(" '<grep %s[^0-9] '.FILE.'' u %s lc '%s' t \"%s",k,condition[0].position[k].label,condition[i].plotting,condition[i].colour,condition[i].ident);
                ind=condition[i].position[k].mode;
                if(ind>=i1) {
                    printf(" $\\\\mu=%.1f$ $s=%.2G$",label_angle[ind],condition[i].position[k].strength);
                    print_jet_classificationGNU(condition[i].position[k].strength);
                }
                printf("\"");
                if(i<num_conditions-1) printf(",");
            }
            for (i=0;i<num_conditions;i++) {
                if(!isnan(condition[i].position[k].est_A*condition[i].position[k].norm)){
                    ind=condition[i].position[k].mode;
                    printf(", G(x,%.16G,%.16G,%.16G,%.16G,%.16G,%.16G) lc '%s' w l t ''",
                           condition[i].position[k].est_A*condition[i].position[k].norm, label_angle[ind], condition[i].position[k].est_sigmasq, condition[i].position[k].offset,label_angle[condition[i].position[k].min_i],label_angle[condition[i].position[k].max_i],condition[i].colour);
                }}
                printf("\n#GNUPLOT\n");*/

        }
            
        }
        printf("\n#GNUPLOT\n");
    }
    
    
    if(0){//print GNUPLOT commands
        //char filePath[PATH_MAX];,fileName[PATH_MAX];
    //rad=RADIUS_PROTRACTOR;
    //if (readlink("/proc/self/fd/1", filePath, sizeof(filePath)-1)!=-1) strncpy(fileName, filePath + 73, PATH_MAX - 73);

        printf("#GNUPLOT set term epslatex standalone color size 15,%g header \"\\\\usepackage{color}\\\\usepackage[dvipsnames]{xcolor}\"\n"
            "#GNUPLOT set out 'CD69negDPKR_jets_gnu_.tex'\n"
            "#GNUPLOT FILE = '%s'\n"
            "#GNUPLOT set format x \"$%%g^\\\\circ$\"\n"
            "#GNUPLOT set xrange [-90:90]\n"
            "#GNUPLOT a=-150\n"
            "#GNUPLOT b=375\n"
            "#GNUPLOT set yrange [a:b]\n"
            "#GNUPLOT set arrow 1 from 45.,a to 45.,b nohead dt '.'\n"
            "#GNUPLOT set arrow 2 from -45.,a to -45.,b nohead dt '.'\n"
            "#GNUPLOT set zeroaxis\n"
            "#GNUPLOT set key bottom\n"
            //"#GNUPLOT max(x,y) = x<y?y:x\n"
            //"#GNUPLOT min(x,y) = x<y?x:y\n"
            "#GNUPLOT g(x,A,mu,sigsq) = A*exp(-0.5*(x-mu)**2/sigsq)\n"
            "#GNUPLOT G(x,A,mu,sigsq,C,a,b) = x<a?1/0: x<b? C+g(x,A,mu,sigsq): 1/0\n"
            "#GNUPLOT set multiplot layout %g,3 title \"Sum of counts per sector minus expected divided"
            " by total counts. KR norm. Number of sectors $n=%d$. Radius $r=%d$bp.\"\n#GNUPLOT\n"
            "#GNUPLOT INF = 0\n#GNUPLOT\n",  ceil(num_positions/3.),panoramic_file,ceil(num_positions/3.)*3,num_sectors,rad);

        int flag_print;
        for(k=0; k<num_positions; k++) {
            flag_print=0;
            for (i=0;i<num_conditions;i++) if(!isnan(condition[i].position[k].est_A*condition[i].position[k].norm)) flag_print++;
            if (flag_print>0) {
            printf("#GNUPLOT set title '%s: chr %d [ %.16G , %.16G ]' offset 0,-1\n",
                   condition[0].position[k].label,condition[0].position[k].chromosome,condition[0].position[k].start, condition[0].position[k].end);
            
            for (i=0;i<num_conditions;i++) {
                ind=condition[i].position[k].mode;
                pan=condition[i].position[k].panoramic_mean[ind];
                if(ind>=i1 && !isnan(pan)) printf("#GNUPLOT set arrow %d nohead from %.1f,0 to %.1f,%.1f dt '-' lc '%s'\n",i+3,label_angle[ind],label_angle[ind],pan,condition[i].colour);
                else printf("#GNUPLOT unset arrow %d\n",i+3);
            }

            printf("#GNUPLOT plot");
            for (i=0;i<num_conditions;i++) {
                printf(" '<grep %s[^0-9] '.FILE.'' u %s lc '%s' t \"%s",condition[0].position[k].label,condition[i].plotting,condition[i].colour,condition[i].ident);
                ind=condition[i].position[k].mode;
                if(ind>=i1) {
                    printf(" $\\\\mu=%.1f$ $s=%.2G$",label_angle[ind],condition[i].position[k].strength);
                    print_jet_classificationGNU(condition[i].position[k].strength);
                }
                printf("\"");
                if(i<num_conditions-1) printf(",");
            }
            for (i=0;i<num_conditions;i++) {
                if(!isnan(condition[i].position[k].est_A*condition[i].position[k].norm)){
                    ind=condition[i].position[k].mode;
                    printf(", G(x,%.16G,%.16G,%.16G,%.16G,%.16G,%.16G) lc '%s' w l t ''",
                           condition[i].position[k].est_A*condition[i].position[k].norm, label_angle[ind], condition[i].position[k].est_sigmasq, condition[i].position[k].offset,label_angle[condition[i].position[k].min_i],label_angle[condition[i].position[k].max_i],condition[i].colour);
                }}
                printf("\n#GNUPLOT\n");
        }
            
        }
        printf("\n#GNUPLOT\n");
    }

    printf("# Info: c'est fini.\n");

    return 0;
    }



int Gaussian_filter(struct condition_struct *condition, int num_sectors) {
    int i,j,k;
    double bandwidth=BANDWIDTH;

    for(k=0;k<condition->num_positions;k++){
        MALLOC(condition->position[k].panoramic_smooth,num_sectors);
        for(i=0;i<num_sectors;i++) {
            //smoothing through Gaussian kernel
            condition->position[k].panoramic_smooth[i]=0.; //initialisation
            for(j=0;j<num_sectors;j++) {
                if(!isnan(condition->position[k].panoramic_mean[j])) condition->position[k].panoramic_smooth[i] += condition->position[k].panoramic_mean[j]*Gaussian((double)(i-j),bandwidth);
            }
        }
    }
    return 0;
}

#define CHOOSE_MAX(ind,i) { prev=condition->position[k].panoramic_smooth[i-1]; current=condition->position[k].panoramic_smooth[i]; next=condition->position[k].panoramic_smooth[i+1]; if(current>prev && current>next) ind=i; }

int find_modes(struct condition_struct *condition, int num_sectors){
    double prev, current,next;
    int i,k;
    int i1 = ceil(0.5*((-JET_OPENING/90.+1.)*(double)num_sectors-1));
    int i2 = floor(0.5*((JET_OPENING/90.+1.)*(double)num_sectors-1));

    for(k=0;k<condition->num_positions;k++){
        int ind=0;
        //scan for local maximum between i1 and i2
        for(i=i1+1;i<i2;i++) {
            if(ind<i1) CHOOSE_MAX(ind,i)
            else if(condition->position[k].panoramic_smooth[i] > condition->position[k].panoramic_smooth[ind]) CHOOSE_MAX(ind,i)
            //printf("p0 ind %d\n",ind);
        }
        
        //the edges can be mode only if they are local maxima
        if(condition->position[k].panoramic_smooth[i1]>condition->position[k].panoramic_smooth[ind]) CHOOSE_MAX(ind,i1)
        if(condition->position[k].panoramic_smooth[i2]>condition->position[k].panoramic_smooth[ind]) CHOOSE_MAX(ind,i2)
        //printf("Pos%d %d %d  %g %g %g  %g %g %g\n",k+1,i,i2,prev,current,next,condition->position[k].panoramic_smooth[i-1],condition->position[k].panoramic_smooth[i],condition->position[k].panoramic_smooth[i+1]);
        //if(k==17) { int j; for(j=50;j<70;j++) printf("%s %d %.16G %.16G\n", condition->position[k].label,j,condition->position[k].panoramic_mean[j],condition->position[k].panoramic_smooth[j]); }
        
        
        //if(k==19) exit(EXIT_FAILURE);
        //if there's no local maximum, ind==0
        condition->position[k].mode=ind;
        }
    return 0;
}

int fit_modes(struct condition_struct *condition, int num_sectors, double *label_angle){
    int i,k,num_sectors_in_fit,j;
    int i1 = ceil(0.5*((-JET_OPENING/90.+1.)*(double)num_sectors-1));
    int i2 = floor(0.5*((JET_OPENING/90.+1.)*(double)num_sectors-1));
    int lb_i=i1-SHORT_FIT_RANGE_IND;
    int ub_i=i2+SHORT_FIT_RANGE_IND;

    printf("#Info: fitting jet candidates in %s\n",condition->ident);

    for(k=0;k<condition->num_positions;k++){
        condition->position[k].norm=0.;
        if(condition->position[k].mode>0){
            condition->position[k].min_i=MAX(lb_i,condition->position[k].mode-FITTING_RANGE);
            //printf("%d min %d %d %d %.16G",k,lb_i,condition->position[k].mode-FITTING_RANGE,condition->position[k].min_i,label_angle[condition->position[k].min_i]);
            condition->position[k].max_i=MIN(ub_i,condition->position[k].mode+FITTING_RANGE);
            //printf(" max %d %d %d %.16G\n",ub_i,condition->position[k].mode+FITTING_RANGE,condition->position[k].max_i,label_angle[condition->position[k].max_i]);

            //calculate norm
            num_sectors_in_fit=0;
            for(i=condition->position[k].min_i;i<=condition->position[k].max_i;i++) if(condition->position[k].panoramic_std[i]>LB_STD){
                condition->position[k].norm+=condition->position[k].panoramic_smooth[i] - condition->position[k].offset;
                num_sectors_in_fit++;
            }
            if(num_sectors_in_fit>2){
                MALLOC(xdata_fit,num_sectors_in_fit);
                MALLOC(ydata_fit,num_sectors_in_fit);
                MALLOC(weights,num_sectors_in_fit);

                //prepare data for fit
                j=0;
                for(i=condition->position[k].min_i;i<=condition->position[k].max_i;i++) if(condition->position[k].panoramic_std[i]>LB_STD){
                    xdata_fit[j]= label_angle[i];
                    ydata_fit[j]= (condition->position[k].panoramic_smooth[i] - condition->position[k].offset)/(condition->position[k].norm);
                    weights[j++]= pow(condition->position[k].panoramic_std[i]/(condition->position[k].norm),-2.);

                }
                //fit
                rgm_gsl_Gauss_fit(xdata_fit, ydata_fit, weights, num_sectors_in_fit, &(condition->position[k].est_A), label_angle[condition->position[k].mode], &(condition->position[k].est_sigmasq), &(condition->position[k].std_A), &(condition->position[k].std_sigmasq));
            
                free(xdata_fit);
                free(ydata_fit);
                free(weights);
            } else { condition->position[k].est_A=condition->position[k].est_sigmasq=condition->position[k].std_A=condition->position[k].std_sigmasq=nan("1"); }
        } else { condition->position[k].est_A=condition->position[k].est_sigmasq=condition->position[k].std_A=condition->position[k].std_sigmasq=nan("1"); }
    }
    return 0;
}

int old_strength(struct condition_struct *condition,int num_sectors){
    int i,k;
    int i1 = ceil(0.5*((-JET_OPENING/90.+1.)*(double)num_sectors-1));
    int i2 = floor(0.5*((JET_OPENING/90.+1.)*(double)num_sectors-1));
    double excess,count;
    
    for(k=0;k<condition->num_positions;k++){
    if(condition->position[k].mode>0){
        //excess
        excess = 0.;
        for(i=i1;i<=i2;i++) {
            count=condition->position[k].panoramic_mean[i];
            if (count>0) excess+=count;
        }
        //strength
        condition->position[k].strength= PREFSTRENGTH * condition->position[k].norm * excess * condition->position[k].panoramic_smooth[condition->position[k].mode] / condition->position[k].est_sigmasq ;
    }}
    return 0;
}

int strength(struct condition_struct *condition,int num_sectors){
    int k;
    for(k=0;k<condition->num_positions;k++) if(condition->position[k].mode>0 && condition->position[k].mode< num_sectors){
        condition->position[k].strength = 0.01*condition->position[k].panoramic_smooth[condition->position[k].mode];
    }
    return 0;
}

int integrate_panoramic(struct condition_struct *condition,int num_sectors){
    double d_theta = 180./(double)num_sectors;
    int i,k;
    int i1 = ceil(0.5*((-JET_OPENING/90.+1.)*(double)num_sectors-1));
    int i2 = floor(0.5*((JET_OPENING/90.+1.)*(double)num_sectors-1));
    //double num_points=251594.; //hard-coded!!
    for(k=0;k<condition->num_positions;k++){
        if(condition->position[k].mode>0){
            condition->position[k].integral_panoramic=0.;
            for(i=i1;i<=i2;i++){
                condition->position[k].integral_panoramic += condition->position[k].panoramic_smooth[i];
            }
            condition->position[k].integral_panoramic *= d_theta; //*2./num_points;
    }}
    return 0;
}
