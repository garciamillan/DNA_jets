#define MAX_INDEX_OF_CHROMOSOMES (99)
#define MAX_CHROMOSOME_NAME (255)
#define MAX_INDEX_OF_POSITIONS (500)
#define PREF (1.e8)
#define NUMBER_CONDITIONS (3)
#define MAX_NUMBER_REPLICATES (10)
#define xxRADIUS_PROTRACTOR (800000.)
#define RADIUS_PROTRACTOR (2000000.)
#define xRADIUS_PROTRACTOR (1000000.)
#define NUM_SECTORS (80)
#define xNUM_SECTORS (60)
#define LAST_CHR (24)
#define TOLWHITELINES (0.1)
#define MAX_LENGTH_BAM (38841)
#define LENGHT_EXP (40000)
#define LBSTRENGHT_STRONG (1.5)
#define LBSTRENGHT_MEDIUM (1.)
#define LB_STD (0.5)
#define LENGTH_STAMP (2000.)
#define UPPER_BOUND_ANGLE_FIT (0.33333*M_PI)
#define UPPER_BOUND_ANGLE_JET (0.25*M_PI)
#define BANDWIDTH (1.)
#define JET_OPENING (45.)
#define SHORT_FIT_RANGE_IND (4)
#define FITTING_RANGE_DEG (27.) //degrees //(6)//(20.)
#define PREFSTRENGTH (1.92e-6)
#define SIGNAL_THRESHOLD (5.) //(15.)
#define MAX_LENGTH_PROJ (1000)
#define CHROMOSOME_OFFSET (600) // HeLa (350) //thymo

#define STENCIL_LENGTH (700) //has to be an integer (lattice points) (600 ->3000kb)
#define STENCIL_WIDTH (20) //has to be an integer (lattice points) (40 ->200kb)
#define STENCIL_MAXPOINTS ((int) ((fmax(STENCIL_LENGTH,STENCIL_WIDTH)+1)*(fmax(STENCIL_LENGTH,STENCIL_WIDTH)+1)))
#define NUM_SLABS (STENCIL_LENGTH/10)


#ifndef PATH_MAX
#define PATH_MAX (1024)
#endif

#define MALLOC(a,n) if ((a=malloc(sizeof(*a)*(n)))==NULL) { fprintf(stderr, "Not enough memory for %s, requested %i bytes, %i items of size %i. %i::%s\n", #a, (int)(sizeof(*a)*n), n, (int)sizeof(*a), errno, strerror(errno)); exit(EXIT_FAILURE); } else { printf("# Info malloc(3)ed %i bytes (%i items of %i bytes) for %s.\n", (int)(sizeof(*a)*(n)), n, (int)sizeof(*a), #a); }

#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )
#define MIN(a,b)  ( ((a)<(b)) ? (a) : (b) )

#define DD fprintf(stderr, "%s::%i\n", __FILE__, __LINE__)

const double resolution=5000.;
int flag_background = 1;

char bed_file[]="/home/clustor2/ma/g/gunsim/DNAfolding/HICdata_2kb_20181026/RUN36_HeLa/new_jet_locations_tcell_cursor.txt";

char jet_analysis_file[]="/home/clustor2/ma/g/gunsim/DNAfolding/HICdata_2kb_20181026/RUN36_HeLa/results/CD69negDPKR_jet_tweaked_merged_431973_431974_431975_431976_chopped.dat";

char jet_consensus_file[]="/home/clustor2/ma/g/gunsim/DNAfolding/HICdata_2kb_20181026/RUN36_HeLa/CD69consensus_eyeballed.dat";

char jet_reach_file[]="/home/clustor2/ma/g/gunsim/DNAfolding/HICdata_2kb_20181026/RUN39_stencil_projections/results/CD69negDPKR_jets_423494_PROJECTION.dat";

char panoramic_file[]="/home/clustor2/ma/g/gunsim/DNAfolding/HICdata_2kb_20181026/RUN36_HeLa/results/CD69negDPKR_jet_tweaked_r1MB_431973_PANORAMIC.dat"; //radius = 1MB, NUM_SECTORS = 60

double radii[]={ 400000., 500000., 600000., 700000., 800000., 900000., 1000000., 1100000., 1200000., 1300000., 1400000., 1500000., 1600000., 1700000., 1800000., 1900000., 2000000., 2100000., 2200000., 2300000., 2400000., 2500000., 2600000., 2700000., 2800000., 2900000., 3000000.};

char *consensus_jets[]={ "Pos27", "Pos34", "Pos42", "Pos57", "Pos59", "Pos67", "Pos83", "Pos90", "Pos91", "Pos99", "Pos139", "Pos152", "Pos177", "Pos179", "Pos185", "Pos188", "Pos195", "Pos221", "Pos234", "Pos248", "Pos250", "Pos254", "Pos259", "Pos271", "Pos274"};

struct BAM_struct {
    double *depth;
    int maxIndex;
    char chromosomeName[MAX_CHROMOSOME_NAME+1];
} ;


/* data_struct[EXPERIMENT = WTR1, WTR2, WTR3, CTCFKOR1, ...]
 * data_struct[WTR1].BAM[chromosome_number]=..
 * data_struct[WTR1].position[77].start;
 * */

struct background {
    //contains the background signal from Dko
    double *stencil;
    
};

struct position_struct {
    //panoramic
    int chromosome;
    double start, end;
    char label[100];
    double *panoramic;
    double *panoramic_mean;
    double *panoramic_std;
    int *whitelines;
    double *panoramic_smooth;
    int mode;
    int min_i;
    int max_i;
    double offset;
    double norm;
    double est_A;
    double est_sigmasq;
    double std_A;
    double std_sigmasq;
    double strength;
    double integral_panoramic;

    //stencil
    int flag_has_stencil;
    int stencil_num_points;
    double angle;
    double *stencil;
    double *stencil_mean;
    double *stencil_std;
    double *stencil_smooth;
    
    //jet depth
    double upstream_reach;
    double downstream_reach;
    double upstream_reach_std;
    double downstream_reach_std;

    //stencil background (given by Dko)
    struct background *background_stencil; //only malloced for Dko
    double *stencil_background_mean;
    double *stencil_background_std;
    double *stencil_background_smooth;

    //stencil subtraction (background given by Dko)
    double *stencil_subtraction_mean;
    double *stencil_subtraction_std;
    double *stencil_subtraction_smooth;

    //jet projection
    int length_up;
    int length_down;
    double *proj_upstream;
    double *proj_downstream;
    double *proj_mean_upstream;
    double *proj_mean_downstream;
    double *proj_std_upstream;
    double *proj_std_downstream;
    double *proj_smooth_upstream;
    double *proj_smooth_downstream;
    double *position_upstream;
    double *position_downstream;

};

struct data_struct {
    // Input (as supplied in initialisation)
    char *ident;
    char *filenameBAM;
    char *basenameEXP; //will need to add "%d_%d_5kb.dat" where %d is the chromosome number
    char *basenameHIC; // add "%d_%d_5kb.bin.gz"
    double total_counts;
    struct condition_struct *parent;

    // Output (as read from data files)
    int test_chr_num;
    int maxIndexOfChromosomes; // This is determined in read_depth, it counts the number of valid chromosomes, with X, Y and unrecognosed chromosome names ignored.
    struct BAM_struct BAM[MAX_INDEX_OF_CHROMOSOMES+1];
    double *expected;
    struct gplist_header_struct *hic;
    int num_positions;
    
    // Output (as calculated)
    struct position_struct *position;
    
} ;


struct data_struct data[]={
    {"WTR1",
        "/home/clustor2/ma/g/gunsim/DNAfolding/HICdata_2kb_20181026/RUN22/bamfiles_depth/CD69negDPWTR1.csv",
        "/home/clustor2/ma/g/gunsim/DNAfolding/HICdata_2kb_20181026/RUN08_Ccode/CD69negDP_replicates/CD69negDPWTR1KRexp_ch_",
        "/home/clustor2/ma/g/gunsim/DNAfolding/HICdata_2kb_20181026/RUN18_jettiness/CD69negDPbinfiles/CD69negDPWTR1KRobs_ch_",
        244438962.
    },
    {"WTR2",
        "/home/clustor2/ma/g/gunsim/DNAfolding/HICdata_2kb_20181026/RUN22/bamfiles_depth/CD69negDPWTR2.csv",
        "/home/clustor2/ma/g/gunsim/DNAfolding/HICdata_2kb_20181026/RUN08_Ccode/CD69negDP_replicates/CD69negDPWTR2KRexp_ch_",
        "/home/clustor2/ma/g/gunsim/DNAfolding/HICdata_2kb_20181026/RUN18_jettiness/CD69negDPbinfiles/CD69negDPWTR2KRobs_ch_",
        223193977.
    },
    {"WTR3",
        "/home/clustor2/ma/g/gunsim/DNAfolding/HICdata_2kb_20181026/RUN22/bamfiles_depth/CD69negDPWTR3.csv",
        "/home/clustor2/ma/g/gunsim/DNAfolding/HICdata_2kb_20181026/RUN08_Ccode/CD69negDP_replicates/CD69negDPWTR3KRexp_ch_",
        "/home/clustor2/ma/g/gunsim/DNAfolding/HICdata_2kb_20181026/RUN18_jettiness/CD69negDPbinfiles/CD69negDPWTR3KRobs_ch_",
        417830413.
    },
    {"WTR4",
        "/home/clustor2/ma/g/gunsim/DNAfolding/HICdata_2kb_20181026/RUN22/bamfiles_depth/CD69negDPWTR4.csv",
        "/home/clustor2/ma/g/gunsim/DNAfolding/HICdata_2kb_20181026/RUN08_Ccode/CD69negDP_replicates/CD69negDPWTR4KRexp_ch_",
        "/home/clustor2/ma/g/gunsim/DNAfolding/HICdata_2kb_20181026/RUN18_jettiness/CD69negDPbinfiles/CD69negDPWTR4KRobs_ch_",
        427257794.
    },
    {"CTCFKOR1",
        "/home/clustor2/ma/g/gunsim/DNAfolding/HICdata_2kb_20181026/RUN22/bamfiles_depth/CD69negDPCTCFKOR1.csv",
        "/home/clustor2/ma/g/gunsim/DNAfolding/HICdata_2kb_20181026/RUN08_Ccode/CD69negDP_replicates/CD69negDPCTCFKOR1KRexp_ch_",
        "/home/clustor2/ma/g/gunsim/DNAfolding/HICdata_2kb_20181026/RUN18_jettiness/CD69negDPbinfiles/CD69negDPCTCFKOR1KRobs_ch_",
        424715008.
    },
    {"CTCFKOR2",
        "/home/clustor2/ma/g/gunsim/DNAfolding/HICdata_2kb_20181026/RUN22/bamfiles_depth/CD69negDPCTCFKOR2.csv",
        "/home/clustor2/ma/g/gunsim/DNAfolding/HICdata_2kb_20181026/RUN08_Ccode/CD69negDP_replicates/CD69negDPCTCFKOR2KRexp_ch_",
        "/home/clustor2/ma/g/gunsim/DNAfolding/HICdata_2kb_20181026/RUN18_jettiness/CD69negDPbinfiles/CD69negDPCTCFKOR2KRobs_ch_",
        325929677.
    },
    {"DKOR1",
        "/home/clustor2/ma/g/gunsim/DNAfolding/HICdata_2kb_20181026/RUN22/bamfiles_depth/CD69negDPDKO1.csv",
        "/home/clustor2/ma/g/gunsim/DNAfolding/HICdata_2kb_20181026/RUN08_Ccode/CD69negDP_replicates/CD69negDPDKOR1KRexp_ch_",
        "/home/clustor2/ma/g/gunsim/DNAfolding/HICdata_2kb_20181026/RUN18_jettiness/CD69negDPbinfiles/CD69negDPDKOR1KRobs_ch_",
        158131390.
    },
    {"DKOR2",
        "/home/clustor2/ma/g/gunsim/DNAfolding/HICdata_2kb_20181026/RUN22/bamfiles_depth/CD69negDPDKO2.csv",
        "/home/clustor2/ma/g/gunsim/DNAfolding/HICdata_2kb_20181026/RUN08_Ccode/CD69negDP_replicates/CD69negDPDKOR2KRexp_ch_",
        "/home/clustor2/ma/g/gunsim/DNAfolding/HICdata_2kb_20181026/RUN18_jettiness/CD69negDPbinfiles/CD69negDPDKOR2KRobs_ch_",
        175281770.
    }
};


struct condition_struct{
    //Input about name and structure of replicates
    char *ident;
    char *ident_short;
    int num_replicates;
    struct data_struct *data[MAX_NUMBER_REPLICATES];
    char *colour;
    char *plotting;
    char *plotting2;
    char *plotting3;
    char *plotting4;

    //Output as calculated
    int num_positions;
    struct position_struct *position;
};

struct condition_struct condition[]={
    {"Wild-type",
        "WT",
        4,
        { &(data[0]),&(data[1]),&(data[2]),&(data[3]) },
        "black",
        "7:8:9 w e pt 4",
        "6:10 w lp pt 4",
        "6:13 w lp pt 4",
        "7:8:9 w e pt 4 lc 'black'"
    },
    {"CTCF ko",
        "CTCFKO",
        2,
        { &(data[4]),&(data[5]) },
        "blue",
        "7:10:11 w e pt 6",
        "6:11 w lp pt 6",
        "6:14 w lp pt 6",
        "7:8:9 w e pt 6 lc 'blue'"
    },
    {"RAD21 CTCF ko",
        "DKO",
        2,
        { &(data[6]),&(data[7]) },
        "red",
        "7:12:13 w e pt 8",
        "6:12 w lp pt 8",
        "6:15 w lp pt 8",
        "7:8:9 w e pt 8 lc 'red'"
    }
};

struct point_protractor{
    double x, y, distance;
    int sector;
};

struct point_stencil{
    double x, y;
    int slab;
};
