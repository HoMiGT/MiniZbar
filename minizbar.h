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

    typedef struct qr_finder_line {
        qr_point pos;
        int len;
        int boffs;
        int eoffs;
    } qr_finder_line;

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


    
} // 结构体定义结束

// 函数声明
namespace minizbar {

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

    extern unsigned long zbar_fourcc_parse(const char* format);

    typedef void (zbar_image_cleanup_handler_t)(zbar_image_t* image);

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
            zbar_image_free_data(_img);
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


