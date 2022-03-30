/*
 * $Header: /Users/garciamillan/Documents/uni/ImperialCL/DNAfolding/HiC_2kb/contactdomains/results36_simulatedJets/DNA_expected.c$
 * COMPILE: cc -O3 -Wall -o exp DNA_expected.c gplist_util2.c -lm -lgsl -lgslcblas
 * EXECUTE: ./exp -a CD69negDPCTCFKOR1_chr1.dat -b /home/clustor2/ma/g/gunsim/DNAfolding/HICdata_2kb_20181026/RUN18_jettiness/CD69negDPbinfiles/CD69negDPDKOR1KRobs_ch_1_1_5kb.bin.gz
 * date: Mon 14 Feb 2022 20:22:12 GMT
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

double *BAMdepth;
int flag_BAMfile=0;

int read_depth(int maxind, FILE *fg, double *BAM_depth);
int read_hic(const char *basenameHIC, struct gplist_header_struct **data, int chr);
int expected_data(struct gplist_header_struct *hic, int flag_BAMfile, double *BAM_depth, const double resolution);

int main(int argc, char *argv[]){
    setlinebuf(stdout);
    char HICfilename[8192];
    char BAMfilename[8192];
    FILE *fg=stdin;
    int ch;
    struct gplist_header_struct *hic_file=NULL;
    
    int i,maxIndex;

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
    
    while ((ch = getopt(argc, argv, "a:b:h")) != -1)
            switch (ch) {
                case 'a':
                    strcpy(BAMfilename, optarg);
                    flag_BAMfile=1;
                  if ((fg=fopen(BAMfilename, "rt"))==NULL) {
                      fprintf(stderr, "Can't open file \"%s\" (%i::%s)\n", BAMfilename, errno, strerror(errno));
                      exit(EXIT_FAILURE);
                  } else printf("#Info: success BAM file \"%s\"\n", BAMfilename);
                  break;
                case 'b':
                    strcpy(HICfilename, optarg);
                    printf("# Info: Reading HiC from file %s\n", HICfilename);
                    if ((hic_file=gplist_read(HICfilename))==NULL) {
                        fprintf(stderr, "# Error: Failed to read %s. (%i::%s)\n", HICfilename, errno, strerror(errno));
                        exit(EXIT_FAILURE);
                    } else printf("#Info: success HiC file \"%s\"\n", HICfilename);
                  break;
                case 'h':
                    fprintf(stderr, "Usage: ./exp -a hic_filename  -b BAM_file (depth of one chromosome)\n");
                    return(0);
                    break;
                default:
                    fprintf(stderr, "Flag %c unknown.\n", ch);
                    break;
            }
    

    printf("# Info: Generating \"expected\" data from \"observed\".\n");
    
    maxIndex = (int) (hic_file->chr_w_length/resolution);
    if(flag_BAMfile) { MALLOC(BAMdepth,maxIndex-CHROMOSOME_OFFSET); }

    if(flag_BAMfile) { read_depth(maxIndex, fg, BAMdepth); }
    
    //printf("%.16G  %.16G\n",BAMdepth[0], BAMdepth[1]);

    expected_data(hic_file, flag_BAMfile, BAMdepth, resolution);
    
    
printf("# Info: c'est fini.\n");

return 0;
}



int read_depth(int maxind, FILE *fg, double *BAM_depth)
{
    //FILE *fg;
    char buffer[8192];
    int locus,depth,index,chrID;

    //DD;
    while (fgets(buffer, sizeof(buffer)-1, fg)!=NULL) {
        //printf("%s",buffer);
        if (sscanf(buffer,"chr%d %d %d",&chrID,&locus,&depth)) {
            index=(int)(locus/resolution-CHROMOSOME_OFFSET);
            if ((index>=0) && (index<maxind-CHROMOSOME_OFFSET)) BAM_depth[index]+=depth;
        }
    }
    fclose(fg);
    
    //for(index=0;index<50;index++) printf("%d  %.16G\n",index,BAM_depth[index]);
    
    printf("# Info: Finished reading BAM file.\n");

    return(0);
}



int expected_data(struct gplist_header_struct *hic, int flag_BAMfile, double *BAM_depth, const double resolution){
    int dist, x_bpr, y_bpr, maxIndex;
    double mom0, mom1;
    
    if (hic==NULL) { fprintf(stderr, "Error: hic==NULL unexpected in line %i of %s\n", __LINE__, __FILE__); exit(EXIT_FAILURE); }
    
    printf("mom1/mom0  dist*resolution  dist  mom0  mom1\n");

    maxIndex = (int) (hic->chr_w_length/resolution);
    for(dist=0;dist<maxIndex;dist++){
        mom0 = mom1 = 0.;
        for(x_bpr=0; x_bpr<(maxIndex-CHROMOSOME_OFFSET) && (x_bpr+dist)<(maxIndex-CHROMOSOME_OFFSET);x_bpr++){
            y_bpr = x_bpr+dist;
             if(flag_BAMfile==0){
                mom0 ++;
                mom1 += gplist_fetch_count(hic,x_bpr,y_bpr);
             } else if (BAM_depth[x_bpr]>0 && BAM_depth[y_bpr]>0) {
                 mom0 ++;
                 mom1 += gplist_fetch_count(hic,x_bpr,y_bpr);
             }
        }
        if (mom0>0) printf("%.16G  %.16G  %d  %.16G  %.16G\n", mom1/mom0, dist*resolution, dist, mom0, mom1);

    }
    
    return 0;
}
