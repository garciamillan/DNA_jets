#ifndef GPLIST_UTIL_H
#define GPLIST_UTIL_H
int gplist_add_edge(int source, int target, double count);
int gplist_print_edges(void);
int gplist_print_edges_sorted(void);
//int gplist_init(double chr_w_length, double Delta_w, int no_coarse_graining);
int gplist_init(double chr_w_length, double Delta_w, char *filename);
int gplist_traversal(int visited);
int gplist_write(char *filename);
int gplist_read(char *filename);
int gplist_write_sorted(char *filename);
int gplist_read_sorted(char *filename);
int gplist_reorganise(void);
int gplist_corr(void);
int gplist_gyration(void);
int gplist_fetch_edge_sequential(int source, int target);
int gplist_fetch_edge_sorted(int source, int target);
double gplist_fetch_count_sorted(int source, int target);
int gplist_prep_header(long int flags, char *command);

#define GPLIST_AVERAGE_EDGES_PER_NODE (2+1500)
#define GPLIST_NN_LINK_COUNT (0.)

struct gplist_edge_struct {
int target; /* target address within the chromosome. */
int next; /* pointer in the linked list */
double count;
int visited;
double corr;
//char unused[64];
long int flags;
char spare[32];
};

typedef struct gplist_edge_struct gplist_edge_strct;

struct gplist_location_struct {
int link_sequential;
int link_sorted;
int visited;
int num_neighbours;
double total_count;
double corr;
long int flags;
char spare[32];
};

typedef struct gplist_location_struct gplist_location_strct;


#define GPLIST_MAGIC "Binary GPLIST Dump"
#define GPLIST_VERSION "V1.0"


#define GPLIST_HDRFLG_SORTED   (1<<0)
#define GPLIST_HDRFLG_SEQUENTIAL (1<<1)
#define GPLIST_HDRFLG_NO_COARSE_GRAINING (1<<2)
#define GPLIST_HDRFLG_COARSE_GRAINING (1<<3)

//#define GPLIST_WITHIN_BOUNDS(a) (((a)>=0) && ((a)<gplist_link_num))
#define GPLIST_WITHIN_BOUNDS_BP(bp)        GPLIST_CHECK_BOUNDS_LOCATION(bp/gplist_header.Delta_w)
#define GPLIST_WITHIN_BOUNDS_LOCATION(loc) (((loc)>=0) && ((loc)<gplist_link_num))


int gplist_no_coarse_graining;


struct gplist_header_struct {
char magic[64];
char version[16];
long int flags;
int link_num;
int max_num_neighbours;
int edge_pool_count;
double chr_w_length;
double Delta_w;
char generating_command[256];
int header_size;
int edge_size;
int location_size;
char spare[1024-64-16-sizeof(long int)-sizeof(int)-sizeof(int)-sizeof(int)-sizeof(double)-sizeof(double)-256-sizeof(int)-sizeof(int)-sizeof(int)];
};

#ifndef GPLIST_UTIL_C
#warning "I think all of this can be removed -- no global vars, have a single header pointer that contains all of this."
extern struct gplist_header_struct gplist_header;
extern gplist_location_strct *gplist_location;
extern gplist_edge_strct *gplist_edge_pool_sorted;
extern int gplist_link_num;
extern double gplist_chr_w_length;
#endif

#endif
