#pragma once
#include "minizbar.h"

// 宏声明定义
namespace minizbar {
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


} // 宏声明定义结束

// 结构体声明
namespace minizbar {
	struct zbar_dcoder_s;
	typedef struct zbar_decoder_s zbar_decoder_t;

	struct zbar_format_def_s;
	typedef struct zbar_format_def_s zbar_format_def_t;

	enum zbar_format_group_e;
	typedef enum zbar_format_group_e zbar_format_group_t;

	enum zbar_color_e;
	typedef enum zbar_color_e zbar_color_t;
 
} // 结构体声明结束

// 结构体与类的定义
namespace minizbar {
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

} // 结构体与类的定义结束


// 函数声明
namespace minizbar {

	int _zbar_qr_decode(qr_reader* reader,zbar_image_scanner_t* iscn,zbar_image_t* img);

	int _zbar_get_symbol_hash(zbar_symbol_type_t sym);

	void cache_sym(zbar_image_scanner_t* iscn,zbar_symbol_t* sym);

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
	zbar_symbol_type_t process_edge(zbar_scanner_t* scn,int y1);

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
	typedef void (conversion_handler_t)(zbar_image_t*,const zbar_format_def_t*,const zbar_image_t*,const zbar_format_def_t*);

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
} // 函数声明结束
