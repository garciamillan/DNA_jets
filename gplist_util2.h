#ifndef GPLIST_UTIL2_H
#define GPLIST_UTIL2_H
struct gplist_header_struct *gplist_read(char *filename);
int gplist_fetch_edge(struct gplist_header_struct *hdr, int source, int target);
double gplist_fetch_count(struct gplist_header_struct *hdr, int source, int target);
int gplist_free(struct gplist_header_struct *hdr);

#define GPLIST_END_OF_LIST (-1)


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
int DEPRECATED_link_sequential;
int link;
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
#define GPLIST_WITHIN_BOUNDS_BP(hdr,bp)        GPLIST_CHECK_BOUNDS_LOCATION(hdr,bp/hdr->Delta_w)
#define GPLIST_WITHIN_BOUNDS_LOCATION(hdr,loc) (((loc)>=0) && ((loc)<hdr->link_num))


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
gplist_location_strct *location;
gplist_edge_strct     *edge_pool;
char spare[1024-64-16-sizeof(long int)-sizeof(int)-sizeof(int)-sizeof(int)-sizeof(double)-sizeof(double)-256-sizeof(int)-sizeof(int)-sizeof(int)-sizeof(gplist_location_strct*)-sizeof(gplist_edge_strct*)];
};

#ifndef GPLIST_UTIL2_C
#endif

#endif
