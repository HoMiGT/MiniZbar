#include "utils.h"
#include <string.h>
using namespace minizbar;


#define TEST_CFG(config, cfg) (((config) >> (cfg)) & 1)

static inline int minizbar::decode_e(unsigned e, unsigned s, unsigned n)
{
    unsigned char E = ((e * n * 2 + 1) / s - 3) / 2;
    return((E >= n - 3) ? -1 : E);
}

static inline unsigned minizbar::pair_width(const zbar_decoder_t* dcode, unsigned char offset)
{
    return (dcode->w[(dcode->idx - offset) & (DECODE_WINDOW - 1)]) + (dcode->w[(dcode->idx - offset+1) & (DECODE_WINDOW - 1)]);
}

zbar_symbol_type_t minizbar::_zbar_find_qr(zbar_decoder_t* dcode)
{
    qr_finder_t* qrf = &dcode->qrf;
    unsigned s, qz, w;
    int ei;
    qrf->s5 -= (dcode->w[(dcode->idx - 6) & (DECODE_WINDOW - 1)]);
    qrf->s5 += (dcode->w[(dcode->idx - 1) & (DECODE_WINDOW - 1)]);
    s = qrf->s5;

    if ((dcode->idx & 1) != ZBAR_SPACE || s < 7) {
        return((zbar_symbol_type_t)0);
    }

    ei = decode_e(pair_width(dcode, 1), s, 7);
    if (ei)
        goto invalid;

    ei = decode_e(pair_width(dcode, 2), s, 7);
    if (ei != 2)
        goto invalid;

    ei = decode_e(pair_width(dcode, 3), s, 7);
    if (ei != 2)
        goto invalid;

    ei = decode_e(pair_width(dcode, 4), s, 7);
    if (ei)
        goto invalid;
    qz = dcode->w[(dcode->idx - 0) & (DECODE_WINDOW - 1)];
    w = dcode->w[(dcode->idx - 1) & (DECODE_WINDOW - 1)]; 
    qrf->line.eoffs = qz + (w + 1) / 2;
    qrf->line.len = qz + w +(dcode->w[(dcode->idx - 2) & (DECODE_WINDOW - 1)]);
    qrf->line.pos[0] = qrf->line.len + (dcode->w[(dcode->idx - 3) & (DECODE_WINDOW - 1)]);
    qrf->line.pos[1] = qrf->line.pos[0];
    w = dcode->w[(dcode->idx - 5) & (DECODE_WINDOW - 1)] ;
    qrf->line.boffs = qrf->line.pos[0] + (dcode->w[(dcode->idx - 4) & (DECODE_WINDOW - 1)]) + (w + 1) / 2;

    dcode->direction = 0;
    dcode->buflen = 0;
    return(ZBAR_QRCODE);

invalid:
    return((zbar_symbol_type_t)0);
}




static inline char release_lock(zbar_decoder_t* dcode, zbar_symbol_type_t req)
{
    zassert(dcode->lock == req, 1, "lock=%d req=%d\n",
        dcode->lock, req);
    dcode->lock = (zbar_symbol_type_t)0;
    return(0);
}

zbar_symbol_type_t minizbar::zbar_decode_width(zbar_decoder_t* dcode, unsigned w)
{
    zbar_symbol_type_t tmp, sym = ZBAR_NONE;

    dcode->w[dcode->idx & (DECODE_WINDOW - 1)] = w;
    
    dcode->s6 -= (dcode->w[(dcode->idx - 7) & (DECODE_WINDOW - 1)]);
    dcode->s6 += (dcode->w[(dcode->idx - 1) & (DECODE_WINDOW - 1)]); 

    if (TEST_CFG(dcode->qrf.config, ZBAR_CFG_ENABLE) && (tmp = _zbar_find_qr(dcode)) > ZBAR_PARTIAL) {
        sym = tmp;
    }

    dcode->idx++;
    dcode->type = sym;
    if (sym) {
        if (dcode->lock && sym > ZBAR_PARTIAL && sym != ZBAR_QRCODE) {
            release_lock(dcode, sym);
        }
        if (dcode->handler) {
            dcode->handler(dcode);
        }
    }
    return(sym);
}

static inline zbar_symbol_type_t minizbar::process_edge(zbar_scanner_t* scn, int y1)
{
    if (!scn->y1_sign)
        scn->last_edge = scn->cur_edge = (1 << ZBAR_FIXED) + ROUND;
    else if (!scn->last_edge)
        scn->last_edge = scn->cur_edge;
    scn->width = scn->cur_edge - scn->last_edge;
    scn->last_edge = scn->cur_edge;

    if (scn->decoder)
        return(zbar_decode_width(scn->decoder, scn->width));
    return(ZBAR_PARTIAL);
}

inline zbar_symbol_type_t minizbar::zbar_scanner_flush(zbar_scanner_t* scn)
{
    unsigned x;
    if (!scn->y1_sign)
        return(ZBAR_NONE);

    x = (scn->x << ZBAR_FIXED) + ROUND;

    if (scn->cur_edge != x || scn->y1_sign > 0) {
        zbar_symbol_type_t edge = process_edge(scn, -scn->y1_sign);
        scn->cur_edge = x;
        scn->y1_sign = -scn->y1_sign;
        return(edge);
    }

    scn->y1_sign = scn->width = 0;
    if (scn->decoder)
        return(zbar_decode_width(scn->decoder, 0));
    return(ZBAR_PARTIAL);
}

zbar_symbol_type_t minizbar::zbar_scanner_new_scan(zbar_scanner_t* scn)
{
    zbar_symbol_type_t edge = ZBAR_NONE;
    while (scn->y1_sign) {
        zbar_symbol_type_t tmp = zbar_scanner_flush(scn);
        if (tmp < 0 || tmp > edge)
            edge = tmp;
    }
    memset(&scn->x, 0, sizeof(zbar_scanner_t) - offsetof(zbar_scanner_t, x));
    scn->y1_thresh = scn->y1_min_thresh;
    if (scn->decoder) {
        zbar_decoder_new_scan(scn->decoder);
    }
    return(edge);
}

zbar_symbol_set_t* minizbar::_zbar_symbol_set_create()
{
    zbar_symbol_set_t* syms = (zbar_symbol_set_t*)calloc(1, sizeof(*syms));
    _zbar_refcnt(&syms->refcnt, 1);
    return (syms);
}

int minizbar::recycle_syms(zbar_image_scanner_t* iscn, zbar_symbol_set_t* syms)
{
    if (_zbar_refcnt(&syms->refcnt, -1))
        return(1);

    _zbar_image_scanner_recycle_syms(iscn, syms->head);
    syms->head = syms->tail = nullptr;
    syms->nsyms = 0;
    return(0);
}

void minizbar::_zbar_image_scanner_recycle_syms(zbar_image_scanner_t* iscn, zbar_symbol_t* sym)
{
    zbar_symbol_t* next = nullptr;
    for (; sym; sym = next) {
        next = sym->next;
        if (sym->refcnt && _zbar_refcnt(&sym->refcnt, -1)) {
            assert(sym->data_alloc);
            sym->next = nullptr;
        }
        else {
            int i;
            recycle_bucket_t* bucket;
            if (!sym->data_alloc) {
                sym->data = nullptr;
                sym->datalen = 0;
            }
            if (sym->syms) {
                if (_zbar_refcnt(&sym->syms->refcnt, -1))
                    assert(0);
                _zbar_image_scanner_recycle_syms(iscn, sym->syms->head);
                sym->syms->head = nullptr;
                _zbar_symbol_set_free(sym->syms);
                sym->syms = nullptr;
            }
            for (i = 0; i < RECYCLE_BUCKETS; i++)
                if (sym->data_alloc < 1 << (i * 2))
                    break;
            if (i == RECYCLE_BUCKETS) {
                assert(sym->data);
                free(sym->data);
                sym->data = nullptr;
                sym->data_alloc = 0;
                i = 0;
            }
            bucket = &iscn->recycle[i];
            bucket->nsyms++;
            sym->next = bucket->head;
            bucket->head = sym;
        }
    }
}

int minizbar::zbar_parse_config(const char* cfgstr, zbar_symbol_type_t* sym, zbar_config_t* cfg, int* val)
{
    const char* dot, * eq;
    int len;
    char negate;
    if (!cfgstr) {
        return (1);
    }
    dot = strchr(cfgstr, '.');
    if (dot) {
        int len = dot - cfgstr;
        if (!len || (len == 1 && !strncmp(cfgstr, "*", len))) {
            *sym = (zbar_symbol_type_t)0;
        }else if (!strncmp(cfgstr, "qrcode", len)) {
            *sym = ZBAR_QRCODE;
        }
        else {
            return (1);
        }
        cfgstr = dot + 1;
    }
    else {
        *sym = (zbar_symbol_type_t)0;
    }
    len = strlen(cfgstr);
    eq = strchr(cfgstr, '=');
    if (eq) {
        len = eq - cfgstr;
    }
    else {
        *val = 1;
    }
    negate = 0;
    if (len > 3 && !strncmp(cfgstr, "no-", 3)) {
        negate = 1;
        cfgstr += 3;
        len -= 3;
    }
    if (len < 1)
        return(1);
    else if (!strncmp(cfgstr, "y-density", len))
        *cfg = ZBAR_CFG_Y_DENSITY;
    else if (!strncmp(cfgstr, "x-density", len))
        *cfg = ZBAR_CFG_X_DENSITY;
    else if (len < 2)
        return(1);
    else if (!strncmp(cfgstr, "enable", len))
        *cfg = ZBAR_CFG_ENABLE;
    else if (len < 3)
        return(1);
    else if (!strncmp(cfgstr, "disable", len)) {
        *cfg = ZBAR_CFG_ENABLE;
        negate = !negate; /* no-disable ?!? */
    }
    else if (!strncmp(cfgstr, "min-length", len))
        *cfg = ZBAR_CFG_MIN_LEN;
    else if (!strncmp(cfgstr, "max-length", len))
        *cfg = ZBAR_CFG_MAX_LEN;
    else if (!strncmp(cfgstr, "ascii", len))
        *cfg = ZBAR_CFG_ASCII;
    else if (!strncmp(cfgstr, "add-check", len))
        *cfg = ZBAR_CFG_ADD_CHECK;
    else if (!strncmp(cfgstr, "emit-check", len))
        *cfg = ZBAR_CFG_EMIT_CHECK;
    else if (!strncmp(cfgstr, "uncertainty", len))
        *cfg = ZBAR_CFG_UNCERTAINTY;
    else if (!strncmp(cfgstr, "position", len))
        *cfg = ZBAR_CFG_POSITION;
    else
        return(1);
    if (eq) {
        *val = strtol(eq + 1, nullptr, 0);
    }
    if (negate) {
        *val = !*val;
    }
    return 0;
}

int minizbar::decoder_set_config_bool(zbar_decoder_t* dcode, zbar_symbol_type_t sym, zbar_config_t cfg, int val)
{
    unsigned* config =  &dcode->qrf.config;
    if (!config || cfg >= ZBAR_CFG_NUM)
        return(1);

    if (!val)
        *config &= ~(1 << cfg);
    else if (val == 1)
        *config |= (1 << cfg);
    else
        return(1);
    return(0);
}

int minizbar::zbar_decoder_set_config(zbar_decoder_t* dcode, zbar_symbol_type_t sym, zbar_config_t cfg, int val)
{
    if (sym == ZBAR_NONE) {
        static const zbar_symbol_type_t all[] = {
         ZBAR_QRCODE,(zbar_symbol_type_t)0
        };
        const zbar_symbol_type_t* symp;
        for (symp = all; *symp; symp++) {
            zbar_decoder_set_config(dcode, *symp, cfg, val);
        }
        return(0);
    }
    if (cfg >= 0 && cfg < ZBAR_CFG_NUM)
        return(decoder_set_config_bool(dcode, sym, cfg, val));
    else if (cfg >= ZBAR_CFG_MIN_LEN && cfg <= ZBAR_CFG_MAX_LEN)
        throw "No Support";
    else
        return(1);
}


static const zbar_format_def_t format_defs[] = {
    { fourcc('Y','8','0','0'), ZBAR_FMT_GRAY, },
    { fourcc('Y','8', 0 , 0), ZBAR_FMT_GRAY, },
    { fourcc('Y','8',' ',' '), ZBAR_FMT_GRAY, },
};


static const int num_format_defs = sizeof(format_defs) / sizeof(zbar_format_def_t);

const zbar_format_def_t* minizbar::_zbar_format_lookup(uint32_t fmt)
{
    const zbar_format_def_t* def = nullptr;
    int i = 0;
    while (i < num_format_defs) {
        def = &format_defs[i];
        if (fmt == def->format)
            return(def);
        i = i * 2 + 1;
        if (fmt > def->format)
            i++;
    }
    return(nullptr);
}

static inline void minizbar::convert_y_resize(zbar_image_t* dst, const zbar_format_def_t* dstfmt, const zbar_image_t* src, const zbar_format_def_t* srcfmt, size_t n)
{
    uint8_t* psrc, * pdst;
    unsigned width, height, xpad, y;

    if (dst->width == src->width && dst->height == src->height) {
        memcpy((void*)dst->data, src->data, n);
        return;
    }
    psrc = static_cast<uint8_t*>(const_cast<void*>(src->data));
    pdst = static_cast<uint8_t*>(const_cast<void*>(dst->data));
    width = (dst->width > src->width) ? src->width : dst->width;
    xpad = (dst->width > src->width) ? dst->width - src->width : 0;
    height = (dst->height > src->height) ? src->height : dst->height;
    for (y = 0; y < height; y++) {
        memcpy(pdst, psrc, width);
        pdst += width;
        psrc += src->width;
        if (xpad) {
            memset(pdst, *(psrc - 1), xpad);
            pdst += xpad;
        }
    }
    psrc -= src->width;
    for (; y < dst->height; y++) {
        memcpy(pdst, psrc, width);
        pdst += width;
        if (xpad) {
            memset(pdst, *(psrc - 1), xpad);
            pdst += xpad;
        }
    }
}

void minizbar::cleanup_ref(zbar_image_t* img)
{
    if (img->next) {
        _zbar_image_refcnt(img->next, -1);
    }
}

void minizbar::convert_copy(zbar_image_t* dst, const zbar_format_def_t* dstfmt, const zbar_image_t* src, const zbar_format_def_t* srcfmt)
{
    if (src->width == dst->width &&
        src->height == dst->height) {
        zbar_image_t* s = (zbar_image_t*)src;
        dst->data = src->data;
        dst->datalen = src->datalen;
        dst->cleanup = cleanup_ref;
        dst->next = s;
        _zbar_image_refcnt(s, 1);
    }
    else {
        convert_y_resize(dst, dstfmt, src, srcfmt, dst->width * dst->height);
    }
        
}

inline void minizbar::zbar_image_free_data(zbar_image_t* img)
{
    if (!img) {
        return;
    }
    if (img->cleanup && img->data) {
        if (img->cleanup != zbar_image_free_data) {
            zbar_image_cleanup_handler_t* cleanup = img->cleanup;
            img->cleanup = zbar_image_free_data;
            cleanup(img);
        }
    }
    else {
        free((void*)img->data);
    }
    img->data = nullptr;
}

static inline void minizbar::_zbar_image_refcnt(zbar_image_t* img, int delta)
{
    if (!_zbar_refcnt(&img->refcnt, delta) && delta <= 0) {
        if (img->cleanup) {
            img->cleanup(img);
        }
    }
}

void minizbar::_zbar_symbol_free(zbar_symbol_t* sym)
{
    if (sym->syms) {
        zbar_symbol_set_ref(sym->syms, -1);
        sym->syms = nullptr;
    }
    if (sym->pts) {
        free(sym->pts);
    }
    if (sym->data_alloc && sym->data) {
        free(sym->data);
    }
    free(sym);
}

static inline int minizbar::_zbar_refcnt(refcnt_t* cnt, int delta)
{
    int rc = (*cnt += delta);
    assert(rc >= 0);
    return (rc);
}

static inline void minizbar::_zbar_symbol_refcnt(zbar_symbol_t* sym, int delta)
{
     if (!_zbar_refcnt(&sym->refcnt, delta) && delta <= 0) {
         _zbar_symbol_free(sym);
     }
}

void minizbar::_zbar_symbol_set_free(zbar_symbol_set_t* syms)
{
    zbar_symbol_t* sym, * next;
    for (sym = syms->head; sym; sym = next) {
        next = sym->next;
        sym->next = nullptr;
        _zbar_symbol_refcnt(sym, -1);
    }
    syms->head = nullptr;
    free(syms);
}

void minizbar::rs_gf256_init(rs_gf256* _gf, unsigned _ppoly)
{
    unsigned p;
    int      i;
    p = 1;
    for (i = 0; i < 256; i++) {
        _gf->exp[i] = _gf->exp[i + 255] = p;
        p = ((p << 1) ^ (-(p >> 7) & _ppoly)) & 0xFF;
    }
    for (i = 0; i < 255; i++) {
        _gf->log[_gf->exp[i]] = i;
    }
    _gf->log[0] = 0;
}

void minizbar::isaac_update(isaac_ctx* _ctx)
{
    unsigned* m;
    unsigned* r;
    unsigned  a;
    unsigned  b;
    unsigned  x;
    unsigned  y;
    int       i;
    m = _ctx->m;
    r = _ctx->r;
    a = _ctx->a;
    b = _ctx->b + (++_ctx->c) & ISAAC_MASK;
    for (i = 0; i < ISAAC_SZ / 2; i++) {
        x = m[i];
        a = (a ^ a << 13) + m[i + ISAAC_SZ / 2] & ISAAC_MASK;
        m[i] = y = m[(x & ISAAC_SZ - 1 << 2) >> 2] + a + b & ISAAC_MASK;
        r[i] = b = m[y >> ISAAC_SZ_LOG + 2 & ISAAC_SZ - 1] + x & ISAAC_MASK;
        x = m[++i];
        a = (a ^ a >> 6) + m[i + ISAAC_SZ / 2] & ISAAC_MASK;
        m[i] = y = m[(x & ISAAC_SZ - 1 << 2) >> 2] + a + b & ISAAC_MASK;
        r[i] = b = m[y >> ISAAC_SZ_LOG + 2 & ISAAC_SZ - 1] + x & ISAAC_MASK;
        x = m[++i];
        a = (a ^ a << 2) + m[i + ISAAC_SZ / 2] & ISAAC_MASK;
        m[i] = y = m[(x & ISAAC_SZ - 1 << 2) >> 2] + a + b & ISAAC_MASK;
        r[i] = b = m[y >> ISAAC_SZ_LOG + 2 & ISAAC_SZ - 1] + x & ISAAC_MASK;
        x = m[++i];
        a = (a ^ a >> 16) + m[i + ISAAC_SZ / 2] & ISAAC_MASK;
        m[i] = y = m[(x & ISAAC_SZ - 1 << 2) >> 2] + a + b & ISAAC_MASK;
        r[i] = b = m[y >> ISAAC_SZ_LOG + 2 & ISAAC_SZ - 1] + x & ISAAC_MASK;
    }
    for (i = ISAAC_SZ / 2; i < ISAAC_SZ; i++) {
        x = m[i];
        a = (a ^ a << 13) + m[i - ISAAC_SZ / 2] & ISAAC_MASK;
        m[i] = y = m[(x & ISAAC_SZ - 1 << 2) >> 2] + a + b & ISAAC_MASK;
        r[i] = b = m[y >> ISAAC_SZ_LOG + 2 & ISAAC_SZ - 1] + x & ISAAC_MASK;
        x = m[++i];
        a = (a ^ a >> 6) + m[i - ISAAC_SZ / 2] & ISAAC_MASK;
        m[i] = y = m[(x & ISAAC_SZ - 1 << 2) >> 2] + a + b & ISAAC_MASK;
        r[i] = b = m[y >> ISAAC_SZ_LOG + 2 & ISAAC_SZ - 1] + x & ISAAC_MASK;
        x = m[++i];
        a = (a ^ a << 2) + m[i - ISAAC_SZ / 2] & ISAAC_MASK;
        m[i] = y = m[(x & ISAAC_SZ - 1 << 2) >> 2] + a + b & ISAAC_MASK;
        r[i] = b = m[y >> ISAAC_SZ_LOG + 2 & ISAAC_SZ - 1] + x & ISAAC_MASK;
        x = m[++i];
        a = (a ^ a >> 16) + m[i - ISAAC_SZ / 2] & ISAAC_MASK;
        m[i] = y = m[(x & ISAAC_SZ - 1 << 2) >> 2] + a + b & ISAAC_MASK;
        r[i] = b = m[y >> ISAAC_SZ_LOG + 2 & ISAAC_SZ - 1] + x & ISAAC_MASK;
    }
    _ctx->b = b;
    _ctx->a = a;
    _ctx->n = ISAAC_SZ;
}

static void minizbar::isaac_mix(unsigned _x[8])
{
    static const unsigned char SHIFT[8] = { 11,2,8,16,10,4,8,9 };
    for (int i = 0; i < 8; ++i) {
        _x[i] ^= _x[i + 1 & 7] << SHIFT[i];
        _x[i + 3 & 7] += _x[i];
        _x[i + 1 & 7] += _x[i + 2 & 7];
        i++;
        _x[i] ^= _x[i + 1 & 7] >> SHIFT[i];
        _x[i + 3 & 7] += _x[i];
        _x[i + 1 & 7] += _x[i + 2 & 7];
    }
}

void minizbar::isaac_init(isaac_ctx* _ctx, const void* _seed, int _nseed)
{
    const unsigned char* seed;
    unsigned* m;
    unsigned* r;
    unsigned             x[8];
    int                  i;
    int                  j;
    _ctx->a = _ctx->b = _ctx->c = 0;
    m = _ctx->m;
    r = _ctx->r;
    x[0] = x[1] = x[2] = x[3] = x[4] = x[5] = x[6] = x[7] = 0x9E3779B9;
    for (i = 0; i < 4; i++) {
        isaac_mix(x);
    }
    if (_nseed > ISAAC_SEED_SZ_MAX)_nseed = ISAAC_SEED_SZ_MAX;
    seed = (const unsigned char*)_seed;
    for (i = 0; i < _nseed >> 2; i++) {
        r[i] = seed[i << 2 | 3] << 24 | seed[i << 2 | 2] << 16 | seed[i << 2 | 1] << 8 | seed[i << 2];
    }
    if (_nseed & 3) {
        r[i] = seed[i << 2];
        for (j = 1; j < (_nseed & 3); j++)r[i] += seed[i << 2 | j] << (j << 3);
        i++;
    }
    memset(r + i, 0, (ISAAC_SZ - i) * sizeof(*r));
    for (i = 0; i < ISAAC_SZ; i += 8) {
        for (j = 0; j < 8; j++)x[j] += r[i + j];
        isaac_mix(x);
        memcpy(m + i, x, sizeof(x));
    }
    for (i = 0; i < ISAAC_SZ; i += 8) {
        for (j = 0; j < 8; j++)x[j] += m[i + j];
        isaac_mix(x);
        memcpy(m + i, x, sizeof(x));
    }
    isaac_update(_ctx);
}

static void minizbar::qr_reader_init(qr_reader* reader)
{
    isaac_init(&reader->isaac, nullptr, 0);
    rs_gf256_init(&reader->gf, QR_PPOLY);
}

qr_reader* minizbar::_zbar_qr_create()
{
    qr_reader* reader = (qr_reader*)calloc(1, sizeof(*reader));
    qr_reader_init(reader);
    return (reader);
}

int minizbar::_zbar_qr_found_line(qr_reader* reader, int dir, const qr_finder_line* line)
{
    qr_finder_lines* lines = &reader->finder_lines[dir];
    if (lines->nlines >= lines->clines) {
        lines->clines *= 2;
        lines->lines = (qr_finder_line*)realloc(lines->lines, ++lines->clines * sizeof(*lines->lines));
    }
    memcpy(lines->lines + lines->nlines++, line, sizeof(*line));
    return (0);
}

unsigned minizbar::zbar_scanner_get_edge(const zbar_scanner_t* scn, unsigned offset, int prec)
{
    unsigned edge = scn->last_edge - offset - (1 << ZBAR_FIXED) - ROUND;
    prec = ZBAR_FIXED - prec;
    if (prec > 0) {
        return (edge >> prec);
    }
    else if (!prec) {
        return (edge);
    }
    else {
        return (edge << -prec);
    }
}

static inline void minizbar::qr_handler(zbar_image_scanner_t* iscn)
{
    unsigned u;
    int vert;
    qr_finder_line* line = &(iscn->dcode->qrf.line);
    assert(line);
    u = zbar_scanner_get_edge(iscn->scn, line->pos[0], QR_FINDER_SUBPREC);
    line->boffs = u - zbar_scanner_get_edge(iscn->scn,line->boffs,QR_FINDER_SUBPREC);
    line->len = zbar_scanner_get_edge(iscn->scn, line->len, QR_FINDER_SUBPREC);
    line->eoffs = zbar_scanner_get_edge(iscn->scn, line->eoffs, QR_FINDER_SUBPREC);
    line->len -= u;
    u = QR_FIXED(iscn->umin, 0) + iscn->du * u;
    if (iscn->du < 0) {
        int tmp = line->boffs;
        line->boffs = line->eoffs;
        line->eoffs = tmp;
        u -= line->len;
    }
    vert = !iscn->dx;
    line->pos[vert] = u;
    line->pos[!vert] = QR_FIXED(iscn->v, 1);
    _zbar_qr_found_line(iscn->qr, vert, line);
}

static void minizbar::symbol_handler(zbar_decoder_t* dcode)
{
    zbar_image_scanner_t* iscn = (zbar_image_scanner_t*)dcode->userdata;
    zbar_symbol_type_t type = dcode->type;
    int x = 0, y = 0, dir;
    const char* data;
    unsigned datalen;
    zbar_symbol_t* sym;
    if (type != ZBAR_QRCODE) {
        throw "NO SUPPORT OTHER TYPE";
    }
    qr_handler(iscn);
}

zbar_decoder_handler_t* minizbar::zbar_decoder_set_handler(zbar_decoder_t* dcode, zbar_decoder_handler_t* handler)
{
    zbar_decoder_handler_t* result = dcode->handler;
    dcode->handler = handler;
    return (result);
}

zbar_symbol_type_t minizbar::zbar_scanner_reset(zbar_scanner_t* scn)
{
    memset(&scn->x, 0, sizeof(zbar_scanner_t) - offsetof(zbar_scanner_t, x));
    scn->y1_thresh = scn->y1_min_thresh;
    if (scn->decoder) {
        zbar_decoder_reset(scn->decoder);
    }
    return(ZBAR_NONE);
}

zbar_scanner_t* minizbar::zbar_scanner_create(zbar_decoder_t* dcode)
{
    zbar_scanner_t* scn = (zbar_scanner_t*)malloc(sizeof(zbar_scanner_t));
    scn->decoder = dcode;
    scn->y1_min_thresh = ZBAR_SCANNER_THRESH_MIN;
    zbar_scanner_reset(scn);
    return nullptr;
}

void minizbar::zbar_decoder_reset(zbar_decoder_t* dcode)
{
    memset(dcode, 0, (long)&dcode->buf_alloc - (long)dcode);
    dcode->qrf.s5 = 0;
}

zbar_decoder_t* minizbar::zbar_decoder_create()
{
    zbar_decoder_t* dcode = (zbar_decoder_t*)calloc(1, sizeof(zbar_decoder_t));
    dcode->buf_alloc = BUFFER_MIN;
    dcode->buf = (unsigned char*)malloc(dcode->buf_alloc);
    dcode->qrf.config = 1 << ZBAR_CFG_ENABLE;
    zbar_decoder_reset(dcode);
    return(dcode);
}
