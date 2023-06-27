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

} // 结构体与类的定义结束


// 函数声明
namespace minizbar {

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
