
enum RTreeMode {
    RTREE_DISK = 0x1,
    RTREE_MEM  = 0x2
};


/*!
 * @brief A node within an R-tree holding the index for data.
 *
 * Note that there are two types of nodes: leaf and twig. Leaf 
 * nodes point to where data actually is. Twig nodes point to 
 * additional index nodes, which may or may not be leaves. Each 
 * of these nodes has additional children, which may span multiple 
 * chromosomes/contigs.
 *
 * With the start/end position, these positions refer specifically 
 * to the chromosomes specified in chr_idx_start/chr_idx_end. Any 
 * chromosomes between these are completely spanned by a given child.
 */
typedef struct RTreeNode RTreeNode_t;
struct RTreeNode {
    uint8_t mode; /* Is this a in-memory or on-disk node */
    uint8_t is_leaf; /* Is this node a leaf */
    uint16_t n_children; /* The number of children of this node, all lists have this length */
    uint32_t *base_start_l; /* A list of the starting chromosome coordinates of each child */
    uint32_t *base_end_l; /* A list of the ending chromosome coordinates of each child */
    uint64_t *record_offset_l; /* A list of the record offset. If the mode is on-disk this will be the offset from the begining of the file */
    union {
        uint64_t *size_l;  /* Leave only: the size of the record block */
        RTreeNode_t **child_l; /* Twigs only */
    } x;
};

/*!
 * A header and index that points to an R-tree that in turn points to data blocks.
 */
typedef struct RTree {
    uint8_t mode;            /* Is this a in-memory or on-disk node */
    uint32_t block_size;     /* The maximum number of children a node can have */ 
    uint32_t base_start;     /* The first position on chr_idx_start with a record */
    uint32_t base_end;       /* The last position on chrIdxEnd with an entry */
    uint64_t index_offset;     /* The offset of the index */
    RTreeNode_t *root;          /* Pointer to the root node */
} RTree_t;

typedef struct {
    uint64_t *record_offset;    /* The offset to the on-disk start of the first record */ 
    uint64_t *index_offset;     /* The offset to the on-disk start of the index. */        
} bgtfZommHdr_t;