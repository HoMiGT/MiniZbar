#pragma once
#include <stdlib.h>
#include <iterator>
#include <string>
#include <assert.h>

// 宏声明以及类型声明定义
namespace minizbar {
    typedef long refcnt_t;
    typedef int qr_point[2];


    #define RECYCLE_BUCKETS (5)
    #define NUM_SYMS (20)
    #define DECODE_WINDOW (16)
    #define NUM_SCN_CFGS (ZBAR_CFG_Y_DENSITY - ZBAR_CFG_X_DENSITY + 1)
    #define ISAAC_SZ_LOG (8)
    #define ISAAC_SZ          (1<<ISAAC_SZ_LOG)

    #define movedelta(dx, dy) do {                  \
        x += (dx);                              \
        y += (dy);                              \
        p += (dx) + ((uintptr_t)(dy) * w);       \
    } while(0);
    
    #define  BUFFER_MIN (0x20)
    #define ZBAR_SCANNER_THRESH_MIN (4)
    #define QR_FINDER_SUBPREC (2)
    #define ZBAR_FIXED (5)
    #define QR_PPOLY (0x1D)
    #define ISAAC_MASK (0xFFFFFFFFU)
    #define ROUND (1 << (ZBAR_FIXED - 1))
    # define QR_FIXED(v, rnd) ((((v) << 1) + (rnd)) << (QR_FINDER_SUBPREC - 1))
    #define ISAAC_SEED_SZ_MAX (ISAAC_SZ<<2)

    #define CFG(iscn, cfg) ((iscn)->configs[(cfg) - ZBAR_CFG_X_DENSITY])

    #define RECYCLE_BUCKETS 5

    #define zbar_fourcc(a, b, c, d)                 \
        ((unsigned long)(a) |                   \
         ((unsigned long)(b) << 8) |            \
         ((unsigned long)(c) << 16) |           \
         ((unsigned long)(d) << 24))

    #define fourcc zbar_fourcc

    # define STAT(x) iscn->stat_##x++

    # define zassert(condition, retval, format, ...) do {   \
        if(!(condition))                                \
            return(retval);                             \
    } while(0)


    #ifndef ZBAR_SCANNER_EWMA_WEIGHT
    # define ZBAR_SCANNER_EWMA_WEIGHT .78
    #endif

    #define EWMA_WEIGHT ((unsigned)((ZBAR_SCANNER_EWMA_WEIGHT              \
                                 * (1 << (ZBAR_FIXED + 1)) + 1) / 2))

    #define ZBAR_SCANNER_THRESH_FADE 8

    #ifndef ZBAR_SCANNER_THRESH_INIT_WEIGHT
    # define ZBAR_SCANNER_THRESH_INIT_WEIGHT .44
    #endif
    #define THRESH_INIT ((unsigned)((ZBAR_SCANNER_THRESH_INIT_WEIGHT       \
                                 * (1 << (ZBAR_FIXED + 1)) + 1) / 2))


    # define ASSERT_POS assert(p == data + x + y * (intptr_t)w)

    #define CACHE_PROXIMITY 1000

    #define CACHE_HYSTERESIS  2000 /* ms */
    #define CACHE_TIMEOUT     (CACHE_HYSTERESIS * 2) /* ms */

    #define QR_MAXI(_a,_b)      ((_a)-((_a)-(_b)&-((_b)>(_a))))
    #define QR_MINI(_a,_b)      ((_a)+((_b)-(_a)&-((_b)<(_a))))

    #define QR_INT_BITS    ((int)sizeof(int)*CHAR_BIT)

    #define QR_SIGNMASK(_x)     (-((_x)<0))
    #define QR_FLIPSIGNI(_a,_b) ((_a)+QR_SIGNMASK(_b)^QR_SIGNMASK(_b))
    #define QR_DIVROUND(_x,_y)  (((_x)+QR_FLIPSIGNI(_y>>1,_x))/(_y))

    #define QR_LARGE_VERSION_SLACK (3)

    #define QR_ALIGN_SUBPREC (2)

    #define QR_FIXMUL(_a,_b,_r,_s) ((int)((_a)*(long long)(_b)+(_r)>>(_s)))

    #define QR_CLAMPI(_a,_b,_c) (QR_MAXI(_a,QR_MINI(_b,_c)))

    #define QR_EXTMUL(_a,_b,_r)    ((_a)*(long long)(_b)+(_r))
    #define QR_HOM_BITS (14)
    #define QR_SMALL_VERSION_SLACK (1)

    #define QR_SORT2I(_a,_b) \
  do{ \
    int t__; \
    t__=QR_MINI(_a,_b)^(_a); \
    (_a)^=t__; \
    (_b)^=t__; \
  } \
  while(0)

    #define QR_ILOG0(_v) (!!((_v)&0x2))
    #define QR_ILOG1(_v) (((_v)&0xC)?2+QR_ILOG0((_v)>>2):QR_ILOG0(_v))
    #define QR_ILOG2(_v) (((_v)&0xF0)?4+QR_ILOG1((_v)>>4):QR_ILOG1(_v))
    #define QR_ILOG3(_v) (((_v)&0xFF00)?8+QR_ILOG2((_v)>>8):QR_ILOG2(_v))
    #define QR_ILOG4(_v) (((_v)&0xFFFF0000)?16+QR_ILOG3((_v)>>16):QR_ILOG3(_v))
    #define QR_ILOG(_v) ((int)QR_ILOG4((unsigned)(_v)))

    #define QR_INT_BITS    ((int)sizeof(int)*CHAR_BIT)
    #define QR_INT_LOGBITS (QR_ILOG(QR_INT_BITS))

    #define QR_M0 (0)

    #define QR_MODE_HAS_DATA(_mode) (!((_mode)&(_mode)-1))

    #define QR_SWAP2I(_a,_b) \
  do{ \
    int t__; \
    t__=(_a); \
    (_a)=(_b); \
    (_b)=t__; \
  } \
  while(0)

    #define MOD(mod) (1 << (mod))


    

}  // 宏声明以及类型声明定义结束


// 结构体声明 以及 enum的声明与定义
namespace minizbar {

    struct zbar_scanner_s;
    typedef struct zbar_scanner_s zbar_scanner_t;

    struct zbar_image_scanner_s;
    typedef struct zbar_image_scanner_s zbar_image_scanner_t;

    struct zbar_image_s;
    typedef struct zbar_image_s zbar_image_t;

    typedef void (zbar_image_data_handler_t)(zbar_image_t* image, const void* userdata);
    typedef void (zbar_image_cleanup_handler_t)(zbar_image_t* image);

    struct zbar_symbol_s;
    typedef struct zbar_symbol_s zbar_symbol_t;

    struct zbar_symbol_set_s;
    typedef struct zbar_symbol_set_s zbar_symbol_set_t;

    typedef enum zbar_symbol_type_e {
        ZBAR_NONE = 0, 
        ZBAR_PARTIAL = 1,  
        ZBAR_EAN2 = 2,  
        ZBAR_EAN5 = 5,  
        ZBAR_EAN8 = 8,  
        ZBAR_UPCE = 9, 
        ZBAR_ISBN10 = 10,  
        ZBAR_UPCA = 12, 
        ZBAR_EAN13 = 13,  
        ZBAR_ISBN13 = 14,  
        ZBAR_COMPOSITE = 15, 
        ZBAR_I25 = 25,  
        ZBAR_DATABAR = 34, 
        ZBAR_DATABAR_EXP = 35,  
        ZBAR_CODABAR = 38, 
        ZBAR_CODE39 = 39,  
        ZBAR_PDF417 = 57,  
        ZBAR_QRCODE = 64,  
        ZBAR_CODE93 = 93, 
        ZBAR_CODE128 = 128,  
        ZBAR_SYMBOL = 0x00ff,
        ZBAR_ADDON2 = 0x0200,
        ZBAR_ADDON5 = 0x0500,
        ZBAR_ADDON = 0x0700,
    } zbar_symbol_type_t;

    struct point_s;
    typedef struct point_s point_t;

    typedef enum zbar_orientation_e {
        ZBAR_ORIENT_UNKNOWN = -1,  
        ZBAR_ORIENT_UP,            
        ZBAR_ORIENT_RIGHT,         
        ZBAR_ORIENT_DOWN,          
        ZBAR_ORIENT_LEFT,           
    } zbar_orientation_t;

    typedef enum zbar_config_e {
        ZBAR_CFG_ENABLE = 0,        /**< enable symbology/feature */
        ZBAR_CFG_ADD_CHECK,         /**< enable check digit when optional */
        ZBAR_CFG_EMIT_CHECK,        /**< return check digit when present */
        ZBAR_CFG_ASCII,             /**< enable full ASCII character set */
        ZBAR_CFG_NUM,               /**< number of boolean decoder configs */

        ZBAR_CFG_MIN_LEN = 0x20,    /**< minimum data length for valid decode */
        ZBAR_CFG_MAX_LEN,           /**< maximum data length for valid decode */

        ZBAR_CFG_UNCERTAINTY = 0x40,/**< required video consistency frames */

        ZBAR_CFG_POSITION = 0x80,   /**< enable scanner to collect position data */

        ZBAR_CFG_X_DENSITY = 0x100, /**< image scanner vertical scan density */
        ZBAR_CFG_Y_DENSITY,         /**< image scanner horizontal scan density */
    } zbar_config_t;

    struct zbar_decoder_s;
    typedef struct zbar_decoder_s zbar_decoder_t;

    struct qr_finder_s;
    typedef struct qr_finder_s qr_finder_t;

    struct recycle_bucket_s;
    typedef struct recycle_bucket_s recycle_bucket_t;

    typedef enum zbar_modifier_e {
        ZBAR_MOD_GS1 = 0,
        ZBAR_MOD_AIM,
        ZBAR_MOD_NUM,
    } zbar_modifier_t;

    struct zbar_dcoder_s;
    typedef struct zbar_decoder_s zbar_decoder_t;

    struct zbar_format_def_s;
    typedef struct zbar_format_def_s zbar_format_def_t;

    enum zbar_format_group_e;
    typedef enum zbar_format_group_e zbar_format_group_t;

    enum zbar_color_e;
    typedef enum zbar_color_e zbar_color_t;

} // 结构体声明与枚举体的声明与定义结束

// 结构体的定义
namespace minizbar {
    struct zbar_scanner_s {
        zbar_decoder_t* decoder; 
        unsigned y1_min_thresh;
        unsigned x;            
        int y0[4];            
        int y1_sign;           
        unsigned y1_thresh;    
        unsigned cur_edge;     
        unsigned last_edge;    
        unsigned width;       
    };

    struct point_s { int x, y; };


    struct zbar_symbol_s {
        zbar_symbol_type_t type;    /* symbol type */
        unsigned int configs;       /* symbology boolean config bitmask */
        unsigned int modifiers;     /* symbology modifier bitmask */
        unsigned int data_alloc;    /* allocation size of data */
        unsigned int datalen;       /* length of binary symbol data */
        char* data;                 /* symbol data */

        unsigned pts_alloc;         /* allocation size of pts */
        unsigned npts;              /* number of points in location polygon */
        point_t* pts;               /* list of points in location polygon */
        zbar_orientation_t orient;  /* coarse orientation */

        refcnt_t refcnt;            /* reference count */
        zbar_symbol_t* next;        /* linked list of results (or siblings) */
        zbar_symbol_set_t* syms;    /* components of composite result */
        unsigned long time;         /* relative symbol capture time */
        int cache_count;            /* cache state */
        int quality;                /* relative symbol reliability metric */
    };


    struct zbar_symbol_set_s {
        refcnt_t refcnt;
        int nsyms;                  /* number of filtered symbols */
        zbar_symbol_t* head;        /* first of decoded symbol results */
        zbar_symbol_t* tail;        /* last of unfiltered symbol results */
    };

    struct zbar_image_s {
        uint32_t format;          
        unsigned width, height;    
        const void* data;         
        unsigned long datalen;     
        unsigned crop_x, crop_y;  
        unsigned crop_w, crop_h;
        void* userdata;           
        zbar_image_cleanup_handler_t* cleanup;
        refcnt_t refcnt;          
        // zbar_video_t* src;         
        int srcidx;                
        zbar_image_t* next;       
        unsigned seq;             
        zbar_symbol_set_t* syms;   
    };

    typedef struct rs_gf256 {
        unsigned char log[256];
        unsigned char exp[511];
    }rs_gf256;

    typedef struct isaac_ctx {
        unsigned n;
        unsigned r[ISAAC_SZ];
        unsigned m[ISAAC_SZ];
        unsigned a;
        unsigned b;
        unsigned c;
    }isaac_ctx;

    typedef struct qr_finder_line {
        qr_point pos;
        int len;
        int boffs;
        int eoffs;
    } qr_finder_line;

    typedef struct qr_finder_lines {
        qr_finder_line* lines;
        int nlines, clines;
    } qr_finder_lines;

    typedef struct qr_reader {
        rs_gf256 gf;
        isaac_ctx isaac;
        qr_finder_lines finder_lines[2];
    }qr_reader;

    struct recycle_bucket_s {
        int nsyms;
        zbar_symbol_t* head;
    };

    struct zbar_image_scanner_s {
        zbar_scanner_t* scn;  
        zbar_decoder_t* dcode;   
        qr_reader* qr;   
        const void* userdata;  
        zbar_image_data_handler_t* handler;
        unsigned long time;     
        zbar_image_t* img;       
        int dx, dy, du, umin, v;   
        zbar_symbol_set_t* syms;    
        recycle_bucket_t recycle[RECYCLE_BUCKETS];
        int enable_cache;          
        zbar_symbol_t* cache; 
        unsigned config;          
        unsigned ean_config;
        int configs[NUM_SCN_CFGS];  
        int sym_configs[1][NUM_SYMS]; 

        int stat_syms_new;
        int stat_iscn_syms_inuse, stat_iscn_syms_recycle;
        int stat_img_syms_inuse, stat_img_syms_recycle;
        int stat_sym_new;
        int stat_sym_recycle[RECYCLE_BUCKETS];
    };

    

    struct qr_finder_s {
        unsigned s5;
        qr_finder_line line;
        unsigned config;
    };

    typedef void (zbar_decoder_handler_t)(zbar_decoder_t* decoder);
    struct zbar_decoder_s {
        unsigned char idx;                 
        unsigned w[DECODE_WINDOW];      
        zbar_symbol_type_t type;           
        zbar_symbol_type_t lock;           
        unsigned modifiers;             
        int direction;  
        unsigned s6;
        unsigned buf_alloc;          
        unsigned buflen;                  
        unsigned char* buf;              
        void* userdata;                    
        zbar_decoder_handler_t* handler;    
        qr_finder_t qrf;                    
    };

    enum zbar_format_group_e {
        ZBAR_FMT_GRAY,
        ZBAR_FMT_YUV_PLANAR,
        ZBAR_FMT_YUV_PACKED,
        ZBAR_FMT_RGB_PACKED,
        ZBAR_FMT_YUV_NV,
        ZBAR_FMT_JPEG,
        ZBAR_FMT_NUM
    };

    struct zbar_format_def_s {
        uint32_t format;
        zbar_format_group_t group;
        union {
            uint8_t gen[4];
            struct {
                uint8_t bpp;
                uint8_t red, green, blue;
            }rgb;
            struct {
                uint8_t xsub2, ysub2;
                uint8_t packorder;
            } yuv;
            uint32_t cmp;
        } p;
    };

    enum zbar_color_e {
        ZBAR_SPACE = 0,
        ZBAR_BAR = 1,
    };


    typedef struct qr_finder_edge_pt {
        qr_point pos;
        int      edge;
        int      extent;
    } qr_finder_edge_pt;

    typedef struct qr_finder_center {
        qr_point           pos;
        qr_finder_edge_pt* edge_pts;
        int                nedge_pts;
    } qr_finder_center;


    typedef struct qr_finder_cluster {
        qr_finder_line** lines;
        int              nlines;
    } qr_finder_cluster;


    typedef enum qr_mode {
        QR_MODE_NUM = 1,
        QR_MODE_ALNUM,
        QR_MODE_STRUCT,
        QR_MODE_BYTE,
        QR_MODE_FNC1_1ST,
        QR_MODE_ECI = 7,
        QR_MODE_KANJI,
        QR_MODE_FNC1_2ND
    }qr_mode;

    typedef struct qr_code_data_entry {
        qr_mode mode;
        union {
            struct {
                unsigned char* buf;
                int            len;
            }data;
            unsigned eci;
            int      ai;
            struct {
                unsigned char sa_index;
                unsigned char sa_size;
                unsigned char sa_parity;
            }sa;
        }payload;
    } qr_code_data_entry;

    typedef struct qr_code_data {
        qr_code_data_entry* entries;
        int                 nentries;
        unsigned char       version;
        unsigned char       ecc_level;
        unsigned char       sa_index;
        unsigned char       sa_size;
        unsigned char       sa_parity;
        unsigned char       self_parity;
        qr_point            bbox[4];
    }qr_code_data;

    typedef struct qr_code_data_list {
        qr_code_data* qrdata;
        int           nqrdata;
        int           cqrdata;
    }  qr_code_data_list;


    struct qr_aff {
        int fwd[2][2];
        int inv[2][2];
        int x0;
        int y0;
        int res;
        int ires;
    };

    struct qr_hom {
        int fwd[3][2];
        int inv[3][2];
        int fwd22;
        int inv22;
        int x0;
        int y0;
        int res;
    };

    struct qr_finder {
        int                size[2];
        int                eversion[2];
        qr_finder_edge_pt* edge_pts[4];
        int                nedge_pts[4];
        int                ninliers[4];
        qr_point           o;
        qr_finder_center* c;
    };

    typedef int qr_line[3];

    struct qr_hom_cell {
        int fwd[3][3];
        int x0;
        int y0;
        int u0;
        int v0;
    };

    typedef struct qr_sampling_grid {
        qr_hom_cell* cells[6];
        unsigned* fpmask;
        int             cell_limits[6];
        int             ncells;
    } qr_sampling_grid;

    typedef struct qr_pack_buf {
        const unsigned char* buf;
        int                  endbyte;
        int                  endbit;
        int                  storage;
    }  qr_pack_buf;

    typedef void* iconv_t;

    typedef enum qr_eci_encoding {
        QR_ECI_GLI0,
        QR_ECI_GLI1,
        QR_ECI_CP437,
        QR_ECI_ISO8859_1,
        QR_ECI_ISO8859_2,
        QR_ECI_ISO8859_3,
        QR_ECI_ISO8859_4,
        QR_ECI_ISO8859_5,
        QR_ECI_ISO8859_6,
        QR_ECI_ISO8859_7,
        QR_ECI_ISO8859_8,
        QR_ECI_ISO8859_9,
        QR_ECI_ISO8859_10,
        QR_ECI_ISO8859_11,
        QR_ECI_ISO8859_13 = QR_ECI_ISO8859_11 + 2,
        QR_ECI_ISO8859_14,
        QR_ECI_ISO8859_15,
        QR_ECI_ISO8859_16,
        QR_ECI_SJIS = 20,
        QR_ECI_UTF8 = 26
    }qr_eci_encoding;
    
} // 结构体定义结束

// 函数声明
namespace minizbar {

    int _zbar_qr_decode(qr_reader* reader, zbar_image_scanner_t* iscn, zbar_image_t* img);

    int _zbar_get_symbol_hash(zbar_symbol_type_t sym);

    void cache_sym(zbar_image_scanner_t* iscn, zbar_symbol_t* sym);

    void _zbar_image_scanner_add_sym(zbar_image_scanner_t* iscn, zbar_symbol_t* sym);

    zbar_symbol_t* _zbar_image_scanner_alloc_sym(zbar_image_scanner_t* iscn, zbar_symbol_type_t type, int datalen);

    zbar_symbol_t* cache_lookup(zbar_image_scanner_t* iscn, zbar_symbol_t* sym);

    void quiet_border(zbar_image_scanner_t* iscn);

    unsigned calc_thresh(zbar_scanner_t* scn);

    zbar_symbol_type_t zbar_scan_y(zbar_scanner_t* scn, int y);

    char release_lock(zbar_decoder_t* dcode, zbar_symbol_type_t req);
    int decode_e(unsigned e, unsigned s, unsigned n);

    unsigned pair_width(const zbar_decoder_t* dcode, unsigned char offset);

    zbar_symbol_type_t _zbar_find_qr(zbar_decoder_t* dcode);

    zbar_symbol_type_t zbar_decode_width(zbar_decoder_t* dcode, unsigned w);
    zbar_symbol_type_t process_edge(zbar_scanner_t* scn, int y1);

    zbar_symbol_type_t zbar_scanner_flush(zbar_scanner_t* scn);
    zbar_symbol_type_t zbar_scanner_new_scan(zbar_scanner_t* scn);
    zbar_symbol_set_t* _zbar_symbol_set_create();

    int recycle_syms(zbar_image_scanner_t* iscn, zbar_symbol_set_t* syms);

    void _zbar_image_scanner_recycle_syms(zbar_image_scanner_t* iscn, zbar_symbol_t* sym);
    int zbar_parse_config(const char* cfgstr, zbar_symbol_type_t* sym, zbar_config_t* cfg, int* val);

    int decoder_set_config_bool(zbar_decoder_t* dcode, zbar_symbol_type_t sym, zbar_config_t cfg, int val);
    int zbar_decoder_set_config(zbar_decoder_t* dcode, zbar_symbol_type_t sym, zbar_config_t cfg, int val);

    const zbar_format_def_t* _zbar_format_lookup(uint32_t fmt);
    void convert_y_resize(zbar_image_t* dst, const zbar_format_def_t* dstfmt, const zbar_image_t* src, const zbar_format_def_t* srcfmt, size_t n);
    static void cleanup_ref(zbar_image_t* img);
    static void convert_copy(zbar_image_t* dst, const zbar_format_def_t* dstfmt, const zbar_image_t* src, const zbar_format_def_t* srcfmt);
    typedef void (conversion_handler_t)(zbar_image_t*, const zbar_format_def_t*, const zbar_image_t*, const zbar_format_def_t*);

    void zbar_image_free_data(zbar_image_t* img);
    void _zbar_image_refcnt(zbar_image_t* img, int delta);

    void _zbar_symbol_free(zbar_symbol_t* sym);
    int _zbar_refcnt(refcnt_t* cnt, int delta);
    void _zbar_symbol_refcnt(zbar_symbol_t* sym, int delta);
    void _zbar_symbol_set_free(zbar_symbol_set_t* syms);

    void rs_gf256_init(rs_gf256* _gf, unsigned _ppoly);
    static void isaac_update(isaac_ctx* _ctx);
    static void isaac_mix(unsigned _x[8]);
    void isaac_init(isaac_ctx* _ctx, const void* _seed, int _nseed);
    static  void qr_reader_init(qr_reader* reader);
    qr_reader* _zbar_qr_create();

    int _zbar_qr_found_line(qr_reader* reader, int dir, const qr_finder_line* line);
    unsigned zbar_scanner_get_edge(const zbar_scanner_t* scn, unsigned offset, int prec);
    static void qr_handler(zbar_image_scanner_t* iscn);
    static void symbol_handler(zbar_decoder_t* dcode);

    zbar_decoder_handler_t* zbar_decoder_set_handler(zbar_decoder_t* dcode, zbar_decoder_handler_t* handler);
    zbar_symbol_type_t zbar_scanner_reset(zbar_scanner_t* scn);
    zbar_scanner_t* zbar_scanner_create(zbar_decoder_t* dcode);
    void zbar_decoder_reset(zbar_decoder_t* dcode);
    zbar_decoder_t* zbar_decoder_create();

    extern zbar_image_scanner_t* zbar_image_scanner_create();
    extern void zbar_image_scanner_destroy(zbar_image_scanner_t* scanner);
    extern zbar_image_data_handler_t* zbar_image_scanner_set_data_handler(zbar_image_scanner_t* scanner, zbar_image_data_handler_t* handler,const void* userdata);

    extern unsigned zbar_symbol_get_loc_size(const zbar_symbol_t* sym);

    extern int zbar_symbol_get_loc_x(const zbar_symbol_t* sym, unsigned idx);

    extern int zbar_symbol_get_loc_y(const zbar_symbol_t* sym, unsigned idx);

    extern void zbar_symbol_ref(const zbar_symbol_t* sym, int refs);

    extern void zbar_symbol_set_ref(const zbar_symbol_set_t* symbols,int refs);

    extern const zbar_symbol_t* zbar_symbol_set_first_symbol(const zbar_symbol_set_t* symbols);

    extern zbar_image_t* zbar_image_create();

    extern void zbar_image_ref(zbar_image_t* image,int refs);

    static inline unsigned long zbar_fourcc_parse(const char* format)
    {
        unsigned long fourcc = 0;
        if (format) {
            for (int i = 0; i < 4 && format[i]; ++i) {
                fourcc |= ((unsigned long)format[i]) << (i * 8);
            }
        }
        return (fourcc);
    }

    extern void zbar_image_set_crop(zbar_image_t* image,unsigned x,unsigned y,unsigned width,unsigned height);

    extern zbar_image_t* zbar_image_convert_resize(const zbar_image_t* image,unsigned long format,unsigned width,unsigned height);

    extern void zbar_image_set_symbols(zbar_image_t* image,const zbar_symbol_set_t* symbols);

    extern int zbar_image_scanner_set_config(zbar_image_scanner_t* scanner,zbar_symbol_type_t symbology,zbar_config_t config,int value);

    extern int zbar_image_scanner_parse_config(zbar_image_scanner_t* scanner, const char* config_string);

    extern void zbar_image_scanner_enable_cache(zbar_image_scanner_t* scanner,int enable);

    extern void zbar_image_scanner_recycle_image(zbar_image_scanner_t* scanner,zbar_image_t* image);

    extern int zbar_scan_image(zbar_image_scanner_t* scanner,zbar_image_t* image);




} // 函数声明结束


// 类声明与定义
namespace minizbar {


    class SymbolIterator;

    class Symbol;

    class SymbolSet;


    class SymbolSet {

    public:
        SymbolSet(const zbar_symbol_set_t* syms = nullptr)
            : _syms(syms) {
            ref();
        }
        SymbolSet(const SymbolSet& syms)
            : _syms(syms._syms)
        {
            ref();
        }

        ~SymbolSet() {
            ref(-1);
        }

        SymbolSet& operator= (const SymbolSet& syms) {
            syms.ref();
            ref(-1);
            _syms = syms._syms;
            return (*this);
        }

        bool operator! () const {
            return (!_syms || !get_size());
        }

        void ref(int delta =1) const{
            if (_syms) {
                zbar_symbol_set_ref((zbar_symbol_set_t*)_syms, delta);
            }
        }

        operator const zbar_symbol_set_t* () const {
            return (_syms);
        }

        int get_size() const {
            return ((_syms) ? _syms->nsyms : 0);
        }

        SymbolIterator symbol_begin() const;

        const SymbolIterator symbol_end() const;


    private:
        const zbar_symbol_set_t* _syms;
    };



    class Symbol {
    public:
        class Point {
        public:
            int x;
            int y;
            Point(){}
            Point(int x, int y) : x(x), y(y) {}
            Point(const Point& pt)
                :x(pt.x),
                y(pt.y)
            {}
            Point& operator= (const Point& pt) {
                x = pt.x;
                y = pt.y;
                return (*this);
            }
        };

        class PointIterator
            : public std::iterator<std::input_iterator_tag, Point> {
        public:
            PointIterator(const Symbol *sym=nullptr, int index=0)
                : _sym(sym),
                _index(index)
            {
                if (sym) {
                    sym->ref(1);
                }
                if (!sym || (unsigned)_index >= zbar_symbol_get_loc_size(*_sym)) {
                    _index = -1;
                }
            }

            PointIterator(const PointIterator& iter)
                :_sym(iter._sym),
                _index(iter._index) {
                if (_sym) {
					_sym->ref();
				}
            }

            ~PointIterator() {
                if (_sym) {
                    _sym->ref(-1);
                }
            }

            PointIterator& operator= (const PointIterator& iter) {
                if (iter._sym) {
                    iter._sym->ref();
                }
                if (_sym) {
                    _sym->ref(-1);
                }
                _sym = iter._sym;
                _index = iter._index;
                return (*this);
            }

            bool operator! () const {
                return (!_sym || _index < 0);
            }

            PointIterator& operator++ () {
                unsigned int i = ++_index;
                if (!_sym || i >= zbar_symbol_get_loc_size(*_sym)) {
                    _index = -1;
                }
                return (*this);
            }

            const Point operator* () const {
                assert(!!*this);
                if (!*this) {
                    return (Point());
                }
                return (Point(zbar_symbol_get_loc_x(*_sym,_index),
                    					zbar_symbol_get_loc_y(*_sym, _index)));
            }


            bool operator== (const PointIterator& iter) const {
                return(_index == iter._index && 
                    ((_index < 0) || _sym == iter._sym));
            }

            bool operator!= (const PointIterator& iter) const {
                return (!(*this == iter));
            }
        private:
            const Symbol* _sym;
            int _index;
        };


        Symbol (const zbar_symbol_t *sym = nullptr)
            : _xmlbuf(nullptr),
            _xmllen(0){
            init(sym);
            ref();
        }

        Symbol(const Symbol& sym) 
            : _sym(sym._sym),
            _type(sym._type),
            _data(sym._data),
            _xmlbuf(nullptr),
            _xmllen(0)

        {
            ref();
        }

        ~Symbol() {
            if (_xmlbuf) {
                free(_xmlbuf);
            }
            ref(-1);
        }

        Symbol& operator= (const Symbol& sym) {
            sym.ref(1);
            ref(-1);
            _sym = sym._sym;
            _type = sym._type;
            _data = sym._data;
            return (*this);
        }

        Symbol& operator= (const zbar_symbol_t* sym) {
            if (sym) {
                zbar_symbol_ref(sym, 1);
            }
            ref(-1);
            init(sym);
            return (*this);
        }

        bool operator! () const {
            return (!_sym);
        }


        void ref(int delta = 1) const {
            if (_sym) {
                zbar_symbol_ref((zbar_symbol_t*)_sym, delta);
            }

        }
        // 类型转换函数 用于将Symbol转换为zbar_symbol_t
        operator const zbar_symbol_t* () const {
            return (_sym);
        }

        bool operator== (const Symbol& sym) const {
            return (_sym == sym._sym);
        }

        bool operator!= (const Symbol& sym) const {
            return (!(*this == sym));
        }

        zbar_symbol_type_t get_type() const {
            return (_type);
        }

        const std::string get_type_name() const {
            return ("QR-Code");
        }

        const std::string get_addon_name() const {
            return ("");
        }

        const std::string get_data() const {
            return (_data);
        }

        unsigned get_data_length() const {
            return ((_sym) ? _sym->datalen : 0);
        }

        int get_count() const {
            return ((_sym) ? _sym->cache_count : 0);
        }

        int get_quality() const {
            return ((_sym) ? _sym->quality : 0);
        }
        
        SymbolSet get_components() const {
            return (SymbolSet((_sym) ? _sym->syms : nullptr));
        }

        PointIterator point_begin() const {
            return (PointIterator(this));
        }

        const PointIterator point_end() const {
            return (PointIterator());
        }

        int get_location_size() const {
            return((_sym) ? _sym->npts : 0);
        }

        int get_location_x(unsigned index) const {
            return ((_sym) ? _sym->pts[index].x : -1);
        }

        int get_location_y(unsigned index) const {
			return ((_sym) ? _sym->pts[index].y : -1);
		}

        int get_orientation() const {
            return ((_sym) ? _sym->orient : -1);
        }

    protected:
        void init(const zbar_symbol_t* sym = nullptr) {
            _sym = sym;
            if (sym) {
                _type = sym->type;
                _data = std::string(sym->data, sym->datalen);
            }
            else {
                _type = ZBAR_NONE;
                _data = "";
            }
        }

    private:
        const zbar_symbol_t* _sym;
        zbar_symbol_type_t _type;
        std::string _data;
        char* _xmlbuf;
        unsigned _xmllen;

    };

    class SymbolIterator 
        : public std::iterator<std::input_iterator_tag, Symbol>
    {
    public:
        SymbolIterator() {}

        SymbolIterator(const SymbolSet& syms)
            : _syms(syms) {
            const zbar_symbol_set_t* zsyms = _syms;
            if (zsyms) {
                _sym = zbar_symbol_set_first_symbol(zsyms);
            }
        }

        ~SymbolIterator() {}

        SymbolIterator& operator= (const SymbolIterator& iter) {
            _syms = iter._syms;
            _sym = iter._sym;
            return (*this);
        }

        bool operator! () const {
			return (!_syms || !_sym);
		}

        SymbolIterator& operator++ () {
            if (!!_sym) {
                const zbar_symbol_t* tmp_sym = _sym;
                _sym = (tmp_sym)? tmp_sym->next : nullptr;
            }
            else if (!!_syms) {
                _sym = zbar_symbol_set_first_symbol(_syms);
            }
            return(*this);
        }

        const Symbol operator* () const {
            return (_sym);
        }

        const Symbol* operator-> () const {
            return (&_sym);
        }

        bool operator == (const SymbolIterator& iter) const {
            return (_sym == iter._sym);
        }

        bool operator!= (const SymbolIterator& iter) const {
            return (!(*this == iter));
        }

        const SymbolIterator end() const {
            return(SymbolIterator());
        }


    private:
        SymbolSet _syms;
        Symbol _sym;

    };

    inline SymbolIterator SymbolSet::symbol_begin() const {
        return (SymbolIterator(*this));
    }

    inline const SymbolIterator SymbolSet::symbol_end() const {
        return (SymbolIterator());
    }



    class Image {
    public:
        class Handler {
        public:
            virtual ~Handler() { }

            virtual void image_callback(Image& image) = 0;

            operator zbar_image_data_handler_t* () const
            {
                return(_cb);
            }

        private:
            static void _cb(zbar_image_t* zimg,
                const void* userdata)
            {
                if (userdata) {
                    Image* image = (Image*)(zimg->userdata);
                    if (image)
                        ((Handler*)userdata)->image_callback(*image);
                    else {
                        Image tmp(zimg, 1);
                        ((Handler*)userdata)->image_callback(tmp);
                    }
                }
            }
        };

        class SymbolIterator : public minizbar::SymbolIterator {
        public:
            SymbolIterator()
                : minizbar::SymbolIterator()
            { }

            SymbolIterator(const SymbolSet& syms)
                : minizbar::SymbolIterator(syms)
            { }

            SymbolIterator(const SymbolIterator& iter)
                : minizbar::SymbolIterator(iter)
            { }
        };

        Image(unsigned width = 0,
            unsigned height = 0,
            const std::string& format = "",
            const void* data = NULL,
            unsigned long length = 0)
            : _img(zbar_image_create())
        {
            _img->userdata = this;
            if (width && height)
                set_size(width, height);
            if (format.length())
                set_format(format);
            if (data && length)
                set_data(data, length);
        }

        ~Image()
        {
            if (_img->userdata == this) {
                _img->userdata = nullptr;
            }
            zbar_image_ref(_img, -1);
        }

        operator const zbar_image_t* () const
        {
            return(_img);
        }

        operator zbar_image_t* ()
        {
            return(_img);
        }

        unsigned long get_format() const
        {
            return(_img->format);
        }

        void set_format(unsigned long format)
        {
            _img->format = format;
        }

        void set_format(const std::string& format)
        {
            unsigned long fourcc = zbar_fourcc_parse(format.c_str());
            _img->format = fourcc;
        }

        unsigned get_sequence() const
        {
            return(_img->seq);
        }

        void set_sequence(unsigned sequence_num)
        {
            _img->seq = sequence_num;
        }

        unsigned get_width() const
        {
            return(_img->width);
        }

        unsigned get_height() const
        {
            return(_img->height);
        }

        void get_size(unsigned& width,
            unsigned& height) const
        {
            width = _img->width;
            height = _img->height;
        }

        void set_size(unsigned width,
            unsigned height)
        {
            _img->crop_x = _img->crop_y = 0;
            _img->width = _img->crop_w = width;
            _img->height = _img->crop_h = height;
        }

        void get_crop(unsigned& x,
            unsigned& y,
            unsigned& width,
            unsigned& height) const
        {
            x = _img->crop_x;
            y = _img->crop_y;
            width = _img->width;
            height = _img->height;
        }

        void set_crop(unsigned x,
            unsigned y,
            unsigned width,
            unsigned height)
        {
            unsigned img_width = _img->width;
            if (x > img_width) {
                x = img_width;
            }
            if (x + width > img_width) {
                width = img_width - x;
            }
            _img->crop_x = x;
            _img->crop_w = width;

            unsigned img_height = _img->height;
            if (y > img_height) {
                y = img_height;
            }
            if (y + height > img_height) {
                height = img_height - y;
            }
            _img->crop_y = y;
            _img->crop_h = height;

        }

        const void* get_data() const
        {
            return(_img->data);
        }

        unsigned long get_data_length() const
        {
            return(_img->datalen);
        }

        void set_data(const void* data,
            unsigned long length)
        {
            //zbar_image_free_data(_img);
            if (!_img) {
                return;
            }
            if (_img->cleanup && _img->data) {
                free((void*)_img->data);
            }
            _img->data = nullptr;

            _img->data = data;
            _img->datalen = length;
            _img->cleanup = _cleanup;

        }

        Image convert(unsigned long format) const
        {
            zbar_image_t* img = zbar_image_convert_resize(_img, format,_img->width,_img->height);
            if (img)
                return(Image(img));
            throw "ImageFormatError";
        }

        Image convert(std::string format) const
        {
            unsigned long fourcc = zbar_fourcc_parse(format.c_str());
            return(convert(fourcc));
        }

        Image convert(unsigned long format,
            unsigned width,
            unsigned height) const
        {
            zbar_image_t* img = zbar_image_convert_resize(_img, format, width, height);
            if (img)
                return(Image(img));
            throw "ImageForamtError";
        }

        const SymbolSet get_symbols() const {
            return(SymbolSet(_img->syms));
        }

        void set_symbols(const SymbolSet& syms) {
            zbar_image_set_symbols(_img, syms);
        }

        SymbolIterator symbol_begin() const {
            return(SymbolIterator(get_symbols()));
        }

        SymbolIterator symbol_end() const {
            return(SymbolIterator());
        }

    protected:

        friend class Video;

        Image(zbar_image_t* src,
            int refs = 0)
            : _img(src)
        {
            if (refs) {
                zbar_image_ref(_img, refs);
            }
            _img->userdata = this;
        }

        static void _cleanup(zbar_image_t* img)
        {
            assert(img);
        }

    private:
        zbar_image_t* _img;
    };


    class ImageScanner {
    public:
        ImageScanner(zbar_image_scanner_t* scanner = nullptr)
        {
            if (scanner)
                _scanner = scanner;
            else
                _scanner = zbar_image_scanner_create();
        }

        ~ImageScanner()
        {
            zbar_image_scanner_destroy(_scanner);
        }
        operator zbar_image_scanner_t* () const
        {
            return(_scanner);
        }

        void set_handler(Image::Handler& handler)
        {
            zbar_image_scanner_set_data_handler(_scanner, handler, &handler);
        }
        int set_config(zbar_symbol_type_t symbology,
            zbar_config_t config,
            int value)
        {
            return(zbar_image_scanner_set_config(_scanner, symbology,
                config, value));
        }
        int set_config(std::string cfgstr)
        {
            return(zbar_image_scanner_parse_config(_scanner, cfgstr.c_str()));
        }
        void enable_cache(bool enable = true)
        {
            zbar_image_scanner_enable_cache(_scanner, enable);
        }
        void recycle_image(Image& image)
        {
            zbar_image_scanner_recycle_image(_scanner, image);
        }
        const SymbolSet get_results() const {
            return(SymbolSet(_scanner->syms));
        }
        int scan(Image& image)
        {
            return(zbar_scan_image(_scanner, image));
        }
        ImageScanner& operator<< (Image& image)
        {
            scan(image);
            return(*this);
        }

    private:
        zbar_image_scanner_t* _scanner;
    };

} // 类声明与定义结束






