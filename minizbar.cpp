#include "minizbar.h"
#include <chrono>
#include <iconv.h>
#include <string.h>

using namespace minizbar;

void zbar_decoder_set_userdata(zbar_decoder_t* dcode,
    void* userdata)
{
    dcode->userdata = userdata;
}

zbar_image_scanner_t* minizbar::zbar_image_scanner_create()
{
    zbar_image_scanner_t* iscn = (zbar_image_scanner_t*)calloc(1, sizeof(zbar_image_scanner_t));
    if (!iscn)
        return(nullptr);
    iscn->dcode = zbar_decoder_create();
    iscn->scn = zbar_scanner_create(iscn->dcode);
    if (!iscn->dcode || !iscn->scn) {
        zbar_image_scanner_destroy(iscn);
        return(nullptr);
    }
    zbar_decoder_set_userdata(iscn->dcode, iscn);
    zbar_decoder_set_handler(iscn->dcode, symbol_handler);

    iscn->qr = _zbar_qr_create();

    CFG(iscn, ZBAR_CFG_X_DENSITY) = 1;
    CFG(iscn, ZBAR_CFG_Y_DENSITY) = 1;
    zbar_image_scanner_set_config(iscn, (zbar_symbol_type_t)0, ZBAR_CFG_POSITION, 1);
    zbar_image_scanner_set_config(iscn, (zbar_symbol_type_t)0, ZBAR_CFG_UNCERTAINTY, 2);
    zbar_image_scanner_set_config(iscn, ZBAR_QRCODE, ZBAR_CFG_UNCERTAINTY, 0);
    zbar_image_scanner_set_config(iscn, ZBAR_CODE128, ZBAR_CFG_UNCERTAINTY, 0);
    zbar_image_scanner_set_config(iscn, ZBAR_CODE93, ZBAR_CFG_UNCERTAINTY, 0);
    zbar_image_scanner_set_config(iscn, ZBAR_CODE39, ZBAR_CFG_UNCERTAINTY, 0);
    zbar_image_scanner_set_config(iscn, ZBAR_CODABAR, ZBAR_CFG_UNCERTAINTY, 1);
    zbar_image_scanner_set_config(iscn, ZBAR_COMPOSITE, ZBAR_CFG_UNCERTAINTY, 0);
    return(iscn);
}

void minizbar::zbar_image_scanner_destroy(zbar_image_scanner_t* iscn)
{
    if (iscn->syms) {
        if (iscn->syms->refcnt) {
            zbar_symbol_set_ref(iscn->syms, -1);
        }
        else {
            _zbar_symbol_set_free(iscn->syms);
        }
        iscn->syms = nullptr;
    }
    if (iscn->scn) {
        free(iscn->scn);
    }
    iscn->scn = nullptr;
    if (iscn->dcode) {
        if (iscn->dcode->buf) {
            free(iscn->dcode->buf);
        }
        iscn->dcode->buf = nullptr;
        free(iscn->dcode);
    }
    iscn->dcode = nullptr;
    for (int i = 0; i < RECYCLE_BUCKETS; ++i) {
        zbar_symbol_t* sym, * next;
        for (sym = iscn->recycle[i].head; sym; sym = next) {
            next = sym->next;
            _zbar_symbol_free(sym);
        }
    }
    if (iscn->qr) {
        if (iscn->qr->finder_lines[0].lines) {
            free(iscn->qr->finder_lines[0].lines);
        }
        iscn->qr->finder_lines[0].lines = nullptr;
        if (iscn->qr->finder_lines[1].lines) {
            free(iscn->qr->finder_lines[1].lines);
        }
        iscn->qr->finder_lines[1].lines = nullptr;
        free(iscn->qr);
    }
    iscn->qr = nullptr;
    free(iscn);
    iscn = nullptr;
}

zbar_image_data_handler_t* minizbar::zbar_image_scanner_set_data_handler(zbar_image_scanner_t* iscn, zbar_image_data_handler_t* handler, const void* userdata)
{
    zbar_image_data_handler_t* result = iscn->handler;
    iscn->handler = handler;
    iscn->userdata = userdata;
    return (result);
}

unsigned minizbar::zbar_symbol_get_loc_size(const zbar_symbol_t* sym)
{
	return (sym->npts);
}

int minizbar::zbar_symbol_get_loc_x(const zbar_symbol_t* sym, unsigned idx)
{
	if (idx < sym->npts) {
		return(sym->pts[idx].x);
	}
	else {
		return (-1);
	}
}

int minizbar::zbar_symbol_get_loc_y(const zbar_symbol_t* sym, unsigned idx)
{
	if (idx < sym->npts) {
		return(sym->pts[idx].y);
	}
	else {
		return (-1);
	}
}

void minizbar::zbar_symbol_ref(const zbar_symbol_t* sym, int refs)
{
    zbar_symbol_t* ncsym = (zbar_symbol_t*)sym;
    _zbar_symbol_refcnt(ncsym, refs);
}

void minizbar::zbar_symbol_set_ref(const zbar_symbol_set_t* syms, int delta)
{
    zbar_symbol_set_t* ncsyms = (zbar_symbol_set_t*)syms;
    if (!_zbar_refcnt(&ncsyms->refcnt, delta) && delta <= 0) {
        _zbar_symbol_set_free(ncsyms);
    }
}

const zbar_symbol_t* minizbar::zbar_symbol_set_first_symbol(const zbar_symbol_set_t* syms)
{
    zbar_symbol_t* sym = syms->tail;
    if (sym) {
        return (sym->next);
    }
    return (syms->head);
}

zbar_image_t* minizbar::zbar_image_create()
{
    zbar_image_t* img = (zbar_image_t*)calloc(1, sizeof(zbar_image_t));
    _zbar_image_refcnt(img, 1);
    img->srcidx = -1;
    return (img);
}

void minizbar::zbar_image_ref(zbar_image_t* image, int refs)
{
    _zbar_image_refcnt(image, refs);
}



void minizbar::zbar_image_set_crop(zbar_image_t* img, unsigned x, unsigned y, unsigned w, unsigned h)
{
    unsigned img_w = img->width;
    if (x > img_w) { x = img_w; }
    if (x + w > img_w) { w = img_w - x; }
    img->crop_x = x;
    img->crop_w = w;

    unsigned img_h = img->height;
    if (y > img_h) {y = img_h;}
    if (y + h > img_h) { h = img_h - y; }
    img->crop_y = y;
    img->crop_h = h;
}

zbar_image_t* minizbar::zbar_image_convert_resize(const zbar_image_t* src, unsigned long fmt, unsigned width, unsigned height)
{
    const zbar_format_def_t* srcfmt, * dstfmt;
    conversion_handler_t* func;
    zbar_image_t* dst = zbar_image_create();
    dst->format = fmt;
    dst->width = width;
    dst->height = height;
    zbar_image_set_crop(dst, src->crop_x, src->crop_y,src->crop_w, src->crop_h);

    if (src->format == fmt &&
        src->width == width &&
        src->height == height) {
        convert_copy(dst, nullptr, src, nullptr);
        return(dst);
    }

    srcfmt = _zbar_format_lookup(src->format);
    dstfmt = _zbar_format_lookup(dst->format);
    if (!srcfmt || !dstfmt)
        return(nullptr);

    if (srcfmt->group == dstfmt->group &&
        srcfmt->p.cmp == dstfmt->p.cmp &&
        src->width == width &&
        src->height == height) {
        convert_copy(dst, nullptr, src, nullptr);
        return(dst);
    }
    else {
        throw "Image formats other than gray are not supported";
    }
}

void minizbar::zbar_image_set_symbols(zbar_image_t* img, const zbar_symbol_set_t* syms)
{
    if (syms) {
        zbar_symbol_set_ref(syms, 1);
    }
    if (img->syms) {
        zbar_symbol_set_ref(img->syms, -1);
    }
    img->syms = (zbar_symbol_set_t*)syms;
}

int minizbar::zbar_image_scanner_set_config(zbar_image_scanner_t* iscn, zbar_symbol_type_t sym, zbar_config_t cfg, int val)
{
    if ((sym == 0 || sym == ZBAR_COMPOSITE) && cfg == ZBAR_CFG_ENABLE) {
        iscn->ean_config = !!val;
        if (sym)
            return(0);
    }

    if (cfg < ZBAR_CFG_UNCERTAINTY)
        return(zbar_decoder_set_config(iscn->dcode, sym, cfg, val));

    if (cfg < ZBAR_CFG_POSITION) {
        int c, i;
        if (cfg > ZBAR_CFG_UNCERTAINTY)
            return(1);
        c = cfg - ZBAR_CFG_UNCERTAINTY;
        if (sym > ZBAR_PARTIAL) {
            i = _zbar_get_symbol_hash(sym);
            iscn->sym_configs[c][i] = val;
        }
        else
            for (i = 0; i < NUM_SYMS; i++)
                iscn->sym_configs[c][i] = val;
        return(0);
    }

    if (sym > ZBAR_PARTIAL)
        return(1);

    if (cfg >= ZBAR_CFG_X_DENSITY && cfg <= ZBAR_CFG_Y_DENSITY) {
        CFG(iscn, cfg) = val;
        return(0);
    }

    if (cfg > ZBAR_CFG_POSITION)
        return(1);
    //cfg -= ZBAR_CFG_POSITION;
    cfg = static_cast<zbar_config_t>(cfg & ~ZBAR_CFG_POSITION);


    if (!val)
        iscn->config &= ~(1 << cfg);
    else if (val == 1)
        iscn->config |= (1 << cfg);
    else
        return(1);

    return(0);
}

static inline int minizbar::zbar_image_scanner_parse_config(zbar_image_scanner_t* scanner, const char* config_string)
{
    zbar_symbol_type_t sym;
    zbar_config_t cfg;
    int val;
    return(zbar_parse_config(config_string, &sym, &cfg, &val) ||
        zbar_image_scanner_set_config(scanner, sym, cfg, val));
}

void minizbar::zbar_image_scanner_enable_cache(zbar_image_scanner_t* iscn, int enable)
{
    if (iscn->cache) {
        _zbar_image_scanner_recycle_syms(iscn, iscn->cache);
        iscn->cache = nullptr;
    }
    iscn->enable_cache = (enable) ? 1 : 0;
}

void minizbar::zbar_image_scanner_recycle_image(zbar_image_scanner_t* iscn, zbar_image_t* img)
{
    zbar_symbol_set_t* syms = iscn->syms;
    if (syms && syms->refcnt) {
        if (recycle_syms(iscn, syms)) {
            STAT(iscn_syms_inuse);
            iscn->syms = nullptr;
        }
        else
            STAT(iscn_syms_recycle);
    }

    syms = img->syms;
    img->syms = nullptr;
    if (syms && recycle_syms(iscn, syms))
        STAT(img_syms_inuse);
    else if (syms) {
        STAT(img_syms_recycle);
        if (iscn->syms)
            _zbar_symbol_set_free(syms);
        else
            iscn->syms = syms;
    }
}

static inline int _zbar_timer_now()
{
    auto now = std::chrono::steady_clock::now();
    auto ms = std::chrono::time_point_cast<std::chrono::milliseconds>(now);
    return ms.time_since_epoch().count();
}

int minizbar::zbar_scan_image(zbar_image_scanner_t* iscn, zbar_image_t* img)
{
    zbar_symbol_set_t* syms;
    const uint8_t* data;
    zbar_scanner_t* scn = iscn->scn;
    unsigned w, h, cx1, cy1;
    int density;
    iscn->time = _zbar_timer_now();

    iscn->qr->finder_lines[0].nlines = 0;
    iscn->qr->finder_lines[1].nlines = 1;

    if (img->format != fourcc('Y', '8', '0', '0') && img->format != fourcc('G', 'R', 'E', 'Y'))
        return (-1);
    iscn->img = img;

    zbar_image_scanner_recycle_image(iscn, img);
    syms = iscn->syms;

    if (!syms) {
        syms = iscn->syms = _zbar_symbol_set_create();
        STAT(syms_new);
        zbar_symbol_set_ref(syms, 1);
    }
    else {
        zbar_symbol_set_ref(syms, 2);
    }
    img->syms = syms;

    w = img->width;
    h = img->height;
    cx1 = img->crop_x + img->crop_w;
    assert(cx1 <= w);
    cy1 = img->crop_y + img->crop_h;
    assert(cy1 <= h);
    data = (const uint8_t*)img->data;

    zbar_scanner_new_scan(scn);
    density = CFG(iscn, ZBAR_CFG_Y_DENSITY);
    if (density > 0) {
        const uint8_t* p = data;
        int x = 0, y = 0;

        int border = (((img->crop_h - 1) % density) + 1) / 2;
        if (border > img->crop_h / 2)
            border = img->crop_h / 2;
        border += img->crop_y;
        assert(border <= h);
        iscn->dy = 0;
        movedelta(img->crop_x, border);
        iscn->v = y;

        while (y < cy1) {
            int cx0 = img->crop_x;
            iscn->dx = iscn->du = 1;
            iscn->umin = cx0;
            while (x < cx1) {
                uint8_t d = *p;
                movedelta(1, 0);
                zbar_scan_y(scn, d);
            }
            ASSERT_POS;
            quiet_border(iscn);
            movedelta(-1, density);
            iscn->v = y;
            if (y >= cy1)
                break;
            iscn->dx = iscn->du = -1;
            iscn->umin = cx1;
            while (x >= cx0) {
                uint8_t d = *p;
                movedelta(-1, 0);
                zbar_scan_y(scn, d);
            }
            ASSERT_POS;
            quiet_border(iscn);
            movedelta(1, density);
            iscn->v = y;
        }

    }
    iscn->dx = 0;
    density = CFG(iscn, ZBAR_CFG_X_DENSITY);
    if (density > 0) {
        const uint8_t* p = data;
        int x = 0, y = 0;

        int border = (((img->crop_w - 1) % density) + 1) / 2;
        if (border > img->crop_w / 2)
            border = img->crop_w / 2;
        border += img->crop_x;
        assert(border <= w);
        movedelta(border, img->crop_y);
        iscn->v = x;

        while (x < cx1) {
            int cy0 = img->crop_y;
            iscn->dy = iscn->du = 1;
            iscn->umin = cy0;
            while (y < cy1) {
                uint8_t d = *p;
                movedelta(0, 1);
                zbar_scan_y(scn, d);
            }
            ASSERT_POS;
            quiet_border(iscn);

            movedelta(density, -1);
            iscn->v = x;
            if (x >= cx1)
                break;

            iscn->dy = iscn->du = -1;
            iscn->umin = cy1;
            while (y >= cy0) {
                uint8_t d = *p;
                movedelta(0, -1);
                zbar_scan_y(scn, d);
            }
            ASSERT_POS;
            quiet_border(iscn);

            movedelta(density, 1);
            iscn->v = x;
        }
    }
    iscn->dy = 0;
    iscn->img = nullptr;

    _zbar_qr_decode(iscn->qr, iscn, img);

    char filter = (!iscn->enable_cache &&
        (density == 1 || CFG(iscn, ZBAR_CFG_Y_DENSITY) == 1));
    int nean = 0, naddon = 0;
    if (syms->nsyms) {
        zbar_symbol_t** symp;
        for (symp = &syms->head; *symp; ) {
            zbar_symbol_t* sym = *symp;
            if (sym->cache_count <= 0 &&
                ((sym->type < ZBAR_COMPOSITE && sym->type > ZBAR_PARTIAL) ||
                    sym->type == ZBAR_DATABAR ||
                    sym->type == ZBAR_DATABAR_EXP ||
                    sym->type == ZBAR_CODABAR))
            {
                if ((sym->type == ZBAR_CODABAR || filter) && sym->quality < 4) {
                    if (iscn->enable_cache) {
                        zbar_symbol_t* entry = cache_lookup(iscn, sym);
                        if (entry)
                            entry->cache_count--;
                        else
                            assert(0);
                    }
                    *symp = sym->next;
                    syms->nsyms--;
                    sym->next = nullptr;
                    _zbar_image_scanner_recycle_syms(iscn, sym);
                    continue;
                }
                else if (sym->type < ZBAR_COMPOSITE &&
                    sym->type != ZBAR_ISBN10)
                {
                    if (sym->type > ZBAR_EAN5)
                        nean++;
                    else
                        naddon++;
                }
            }
            symp = &sym->next;
        }

        if (nean == 1 && naddon == 1 && iscn->ean_config) {
            zbar_symbol_t* ean = nullptr, * addon = nullptr;
            for (symp = &syms->head; *symp; ) {
                zbar_symbol_t* sym = *symp;
                if (sym->type < ZBAR_COMPOSITE && sym->type > ZBAR_PARTIAL) {
                    *symp = sym->next;
                    syms->nsyms--;
                    sym->next = nullptr;
                    if (sym->type <= ZBAR_EAN5)
                        addon = sym;
                    else
                        ean = sym;
                }
                else
                    symp = &sym->next;
            }
            assert(ean);
            assert(addon);

            int datalen = ean->datalen + addon->datalen + 1;
            zbar_symbol_t* ean_sym =
                _zbar_image_scanner_alloc_sym(iscn, ZBAR_COMPOSITE, datalen);
            ean_sym->orient = ean->orient;
            ean_sym->syms = _zbar_symbol_set_create();
            memcpy(ean_sym->data, ean->data, ean->datalen);
            memcpy(ean_sym->data + ean->datalen,
                addon->data, addon->datalen + 1);
            ean_sym->syms->head = ean;
            ean->next = addon;
            ean_sym->syms->nsyms = 2;
            _zbar_image_scanner_add_sym(iscn, ean_sym);
        }
    }
    if (syms->nsyms && iscn->handler)
        iscn->handler(img, iscn->userdata);
    return(syms->nsyms);
}

#define TEST_CFG(config, cfg) (((config) >> (cfg)) & 1)

static inline int minizbar::decode_e(unsigned e, unsigned s, unsigned n)
{
    unsigned char E = ((e * n * 2 + 1) / s - 3) / 2;
    return((E >= n - 3) ? -1 : E);
}

static inline unsigned minizbar::pair_width(const zbar_decoder_t* dcode, unsigned char offset)
{
    return (dcode->w[(dcode->idx - offset) & (DECODE_WINDOW - 1)]) + (dcode->w[(dcode->idx - offset + 1) & (DECODE_WINDOW - 1)]);
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
    qrf->line.len = qz + w + (dcode->w[(dcode->idx - 2) & (DECODE_WINDOW - 1)]);
    qrf->line.pos[0] = qrf->line.len + (dcode->w[(dcode->idx - 3) & (DECODE_WINDOW - 1)]);
    qrf->line.pos[1] = qrf->line.pos[0];
    w = dcode->w[(dcode->idx - 5) & (DECODE_WINDOW - 1)];
    qrf->line.boffs = qrf->line.pos[0] + (dcode->w[(dcode->idx - 4) & (DECODE_WINDOW - 1)]) + (w + 1) / 2;

    dcode->direction = 0;
    dcode->buflen = 0;
    return(ZBAR_QRCODE);

invalid:
    return((zbar_symbol_type_t)0);
}


static int qr_finder_cluster_lines(qr_finder_cluster* _clusters,
    qr_finder_line** _neighbors, qr_finder_line* _lines, int _nlines, int _v) {
    unsigned char* mark;
    qr_finder_line** neighbors;
    int              nneighbors;
    int              nclusters;
    int              i;
    mark = (unsigned char*)calloc(_nlines, sizeof(*mark));
    neighbors = _neighbors;
    nclusters = 0;
    for (i = 0; i < _nlines - 1; i++)if (!mark[i]) {
        int len;
        int j;
        nneighbors = 1;
        neighbors[0] = _lines + i;
        len = _lines[i].len;
        for (j = i + 1; j < _nlines; j++)if (!mark[j]) {
            const qr_finder_line* a;
            const qr_finder_line* b;
            int                   thresh;
            a = neighbors[nneighbors - 1];
            b = _lines + j;
            thresh = a->len + 7 >> 2;
            if (abs(a->pos[1 - _v] - b->pos[1 - _v]) > thresh)break;
            if (abs(a->pos[_v] - b->pos[_v]) > thresh)continue;
            if (abs(a->pos[_v] + a->len - b->pos[_v] - b->len) > thresh)continue;
            if (a->boffs > 0 && b->boffs > 0 &&
                abs(a->pos[_v] - a->boffs - b->pos[_v] + b->boffs) > thresh) {
                continue;
            }
            if (a->eoffs > 0 && b->eoffs > 0 &&
                abs(a->pos[_v] + a->len + a->eoffs - b->pos[_v] - b->len - b->eoffs) > thresh) {
                continue;
            }
            neighbors[nneighbors++] = _lines + j;
            len += b->len;
        }
        if (nneighbors < 3)continue;
        len = ((len << 1) + nneighbors) / (nneighbors << 1);
        if (nneighbors * (5 << QR_FINDER_SUBPREC) >= len) {
            _clusters[nclusters].lines = neighbors;
            _clusters[nclusters].nlines = nneighbors;
            for (j = 0; j < nneighbors; j++)mark[neighbors[j] - _lines] = 1;
            neighbors += nneighbors;
            nclusters++;
        }
    }
    free(mark);
    return nclusters;
}

static int qr_finder_vline_cmp(const void* _a, const void* _b) {
    const qr_finder_line* a;
    const qr_finder_line* b;
    a = (const qr_finder_line*)_a;
    b = (const qr_finder_line*)_b;
    return ((a->pos[0] > b->pos[0]) - (a->pos[0] < b->pos[0]) << 1) +
        (a->pos[1] > b->pos[1]) - (a->pos[1] < b->pos[1]);
}

static int qr_finder_lines_are_crossing(const qr_finder_line* _hline,
    const qr_finder_line* _vline) {
    return
        _hline->pos[0] <= _vline->pos[0] && _vline->pos[0] < _hline->pos[0] + _hline->len &&
        _vline->pos[1] <= _hline->pos[1] && _hline->pos[1] < _vline->pos[1] + _vline->len;
}

static int qr_finder_edge_pts_fill(qr_finder_edge_pt* _edge_pts, int _nedge_pts,
    qr_finder_cluster** _neighbors, int _nneighbors, int _v) {
    int i;
    for (i = 0; i < _nneighbors; i++) {
        qr_finder_cluster* c;
        int                j;
        c = _neighbors[i];
        for (j = 0; j < c->nlines; j++) {
            qr_finder_line* l;
            l = c->lines[j];
            if (l->boffs > 0) {
                _edge_pts[_nedge_pts].pos[0] = l->pos[0];
                _edge_pts[_nedge_pts].pos[1] = l->pos[1];
                _edge_pts[_nedge_pts].pos[_v] -= l->boffs;
                _nedge_pts++;
            }
            if (l->eoffs > 0) {
                _edge_pts[_nedge_pts].pos[0] = l->pos[0];
                _edge_pts[_nedge_pts].pos[1] = l->pos[1];
                _edge_pts[_nedge_pts].pos[_v] += l->len + l->eoffs;
                _nedge_pts++;
            }
        }
    }
    return _nedge_pts;
}

static int qr_finder_center_cmp(const void* _a, const void* _b) {
    const qr_finder_center* a;
    const qr_finder_center* b;
    a = (const qr_finder_center*)_a;
    b = (const qr_finder_center*)_b;
    return ((b->nedge_pts > a->nedge_pts) - (b->nedge_pts < a->nedge_pts) << 2) +
        ((a->pos[1] > b->pos[1]) - (a->pos[1] < b->pos[1]) << 1) +
        (a->pos[0] > b->pos[0]) - (a->pos[0] < b->pos[0]);
}

static int qr_finder_find_crossings(qr_finder_center* _centers,
    qr_finder_edge_pt* _edge_pts, qr_finder_cluster* _hclusters, int _nhclusters,
    qr_finder_cluster* _vclusters, int _nvclusters) {
    qr_finder_cluster** hneighbors;
    qr_finder_cluster** vneighbors;
    unsigned char* hmark;
    unsigned char* vmark;
    int                 ncenters;
    int                 i;
    int                 j;
    hneighbors = (qr_finder_cluster**)malloc(_nhclusters * sizeof(*hneighbors));
    vneighbors = (qr_finder_cluster**)malloc(_nvclusters * sizeof(*vneighbors));
    hmark = (unsigned char*)calloc(_nhclusters, sizeof(*hmark));
    vmark = (unsigned char*)calloc(_nvclusters, sizeof(*vmark));
    ncenters = 0;
    for (i = 0; i < _nhclusters; i++)if (!hmark[i]) {
        qr_finder_line* a;
        qr_finder_line* b;
        int             nvneighbors;
        int             nedge_pts;
        int             y;
        a = _hclusters[i].lines[_hclusters[i].nlines >> 1];
        y = nvneighbors = 0;
        for (j = 0; j < _nvclusters; j++)if (!vmark[j]) {
            b = _vclusters[j].lines[_vclusters[j].nlines >> 1];
            if (qr_finder_lines_are_crossing(a, b)) {
                vmark[j] = 1;
                y += (b->pos[1] << 1) + b->len;
                if (b->boffs > 0 && b->eoffs > 0)y += b->eoffs - b->boffs;
                vneighbors[nvneighbors++] = _vclusters + j;
            }
        }
        if (nvneighbors > 0) {
            qr_finder_center* c;
            int               nhneighbors;
            int               x;
            x = (a->pos[0] << 1) + a->len;
            if (a->boffs > 0 && a->eoffs > 0)x += a->eoffs - a->boffs;
            hneighbors[0] = _hclusters + i;
            nhneighbors = 1;
            j = nvneighbors >> 1;
            b = vneighbors[j]->lines[vneighbors[j]->nlines >> 1];
            for (j = i + 1; j < _nhclusters; j++)if (!hmark[j]) {
                a = _hclusters[j].lines[_hclusters[j].nlines >> 1];
                if (qr_finder_lines_are_crossing(a, b)) {
                    hmark[j] = 1;
                    x += (a->pos[0] << 1) + a->len;
                    if (a->boffs > 0 && a->eoffs > 0)x += a->eoffs - a->boffs;
                    hneighbors[nhneighbors++] = _hclusters + j;
                }
            }
            c = _centers + ncenters++;
            c->pos[0] = (x + nhneighbors) / (nhneighbors << 1);
            c->pos[1] = (y + nvneighbors) / (nvneighbors << 1);
            c->edge_pts = _edge_pts;
            nedge_pts = qr_finder_edge_pts_fill(_edge_pts, 0,
                hneighbors, nhneighbors, 0);
            nedge_pts = qr_finder_edge_pts_fill(_edge_pts, nedge_pts,
                vneighbors, nvneighbors, 1);
            c->nedge_pts = nedge_pts;
            _edge_pts += nedge_pts;
        }
    }
    free(vmark);
    free(hmark);
    free(vneighbors);
    free(hneighbors);
    qsort(_centers, ncenters, sizeof(*_centers), qr_finder_center_cmp);
    return ncenters;
}

static int qr_finder_centers_locate(qr_finder_center** _centers,
    qr_finder_edge_pt** _edge_pts, qr_reader* reader,
    int _width, int _height) {
    qr_finder_line* hlines = reader->finder_lines[0].lines;
    int                 nhlines = reader->finder_lines[0].nlines;
    qr_finder_line* vlines = reader->finder_lines[1].lines;
    int                 nvlines = reader->finder_lines[1].nlines;

    qr_finder_line** hneighbors;
    qr_finder_cluster* hclusters;
    int                 nhclusters;
    qr_finder_line** vneighbors;
    qr_finder_cluster* vclusters;
    int                 nvclusters;
    int                 ncenters;

    hneighbors = (qr_finder_line**)malloc(nhlines * sizeof(*hneighbors));
    hclusters = (qr_finder_cluster*)malloc((nhlines >> 1) * sizeof(*hclusters));
    nhclusters = qr_finder_cluster_lines(hclusters, hneighbors, hlines, nhlines, 0);
    qsort(vlines, nvlines, sizeof(*vlines), qr_finder_vline_cmp);
    vneighbors = (qr_finder_line**)malloc(nvlines * sizeof(*vneighbors));
    vclusters = (qr_finder_cluster*)malloc((nvlines >> 1) * sizeof(*vclusters));
    nvclusters = qr_finder_cluster_lines(vclusters, vneighbors, vlines, nvlines, 1);
    if (nhclusters >= 3 && nvclusters >= 3) {
        qr_finder_edge_pt* edge_pts;
        qr_finder_center* centers;
        int                 nedge_pts;
        int                 i;
        nedge_pts = 0;
        for (i = 0; i < nhclusters; i++)nedge_pts += hclusters[i].nlines;
        for (i = 0; i < nvclusters; i++)nedge_pts += vclusters[i].nlines;
        nedge_pts <<= 1;
        edge_pts = (qr_finder_edge_pt*)malloc(nedge_pts * sizeof(*edge_pts));
        centers = (qr_finder_center*)malloc(
            QR_MINI(nhclusters, nvclusters) * sizeof(*centers));
        ncenters = qr_finder_find_crossings(centers, edge_pts,
            hclusters, nhclusters, vclusters, nvclusters);
        *_centers = centers;
        *_edge_pts = edge_pts;
    }
    else ncenters = 0;
    free(vclusters);
    free(vneighbors);
    free(hclusters);
    free(hneighbors);
    return ncenters;
}

unsigned char* qr_binarize(const unsigned char* _img, int _width, int _height) {
    unsigned char* mask = nullptr;
    if (_width > 0 && _height > 0) {
        unsigned* col_sums;
        int            logwindw;
        int            logwindh;
        int            windw;
        int            windh;
        int            y0offs;
        int            y1offs;
        unsigned       g;
        int            x;
        int            y;
        mask = (unsigned char*)malloc(_width * _height * sizeof(*mask));
        for (logwindw = 4; logwindw < 8 && (1 << logwindw) < (_width + 7 >> 3); logwindw++);
        for (logwindh = 4; logwindh < 8 && (1 << logwindh) < (_height + 7 >> 3); logwindh++);
        windw = 1 << logwindw;
        windh = 1 << logwindh;
        col_sums = (unsigned*)malloc(_width * sizeof(*col_sums));
        for (x = 0; x < _width; x++) {
            g = _img[x];
            col_sums[x] = (g << logwindh - 1) + g;
        }
        for (y = 1; y < (windh >> 1); y++) {
            y1offs = QR_MINI(y, _height - 1) * _width;
            for (x = 0; x < _width; x++) {
                g = _img[y1offs + x];
                col_sums[x] += g;
            }
        }
        for (y = 0; y < _height; y++) {
            unsigned m;
            int      x0;
            int      x1;
            m = (col_sums[0] << logwindw - 1) + col_sums[0];
            for (x = 1; x < (windw >> 1); x++) {
                x1 = QR_MINI(x, _width - 1);
                m += col_sums[x1];
            }
            for (x = 0; x < _width; x++) {
                g = _img[y * _width + x];
                mask[y * _width + x] = -(g + 3 << logwindw + logwindh < m) & 0xFF;
                if (x + 1 < _width) {
                    x0 = QR_MAXI(0, x - (windw >> 1));
                    x1 = QR_MINI(x + (windw >> 1), _width - 1);
                    m += col_sums[x1] - col_sums[x0];
                }
            }
            /*Update the column sums.*/
            if (y + 1 < _height) {
                y0offs = QR_MAXI(0, y - (windh >> 1)) * _width;
                y1offs = QR_MINI(y + (windh >> 1), _height - 1) * _width;
                for (x = 0; x < _width; x++) {
                    col_sums[x] -= _img[y0offs + x];
                    col_sums[x] += _img[y1offs + x];
                }
            }
        }
        free(col_sums);
    }
    return(mask);
}

inline void qr_code_data_list_init(qr_code_data_list* _qrlist) {
    _qrlist->qrdata = nullptr;
    _qrlist->nqrdata = _qrlist->cqrdata = 0;
}


static int qr_point_ccw(const qr_point _p0,
    const qr_point _p1, const qr_point _p2) {
    return (_p1[0] - _p0[0]) * (_p2[1] - _p0[1]) - (_p1[1] - _p0[1]) * (_p2[0] - _p0[0]);
}

static unsigned qr_point_distance2(const qr_point _p1, const qr_point _p2) {
    return (_p1[0] - _p2[0]) * (_p1[0] - _p2[0]) + (_p1[1] - _p2[1]) * (_p1[1] - _p2[1]);
}

int qr_ilog(unsigned _v) {
    int ret;
    int m;
    m = !!(_v & 0xFFFF0000) << 4;
    _v >>= m;
    ret = m;
    m = !!(_v & 0xFF00) << 3;
    _v >>= m;
    ret |= m;
    m = !!(_v & 0xF0) << 2;
    _v >>= m;
    ret |= m;
    m = !!(_v & 0xC) << 1;
    _v >>= m;
    ret |= m;
    ret |= !!(_v & 0x2);
    return ret + !!_v;
}

static void qr_aff_init(qr_aff* _aff,
    const qr_point _p0, const qr_point _p1, const qr_point _p2, int _res) {
    int det;
    int ires;
    int dx1;
    int dy1;
    int dx2;
    int dy2;
    /*det is ensured to be positive by our caller.*/
    dx1 = _p1[0] - _p0[0];
    dx2 = _p2[0] - _p0[0];
    dy1 = _p1[1] - _p0[1];
    dy2 = _p2[1] - _p0[1];
    det = dx1 * dy2 - dy1 * dx2;
    ires = QR_MAXI((qr_ilog(abs(det)) >> 1) - 2, 0);
    _aff->fwd[0][0] = dx1;
    _aff->fwd[0][1] = dx2;
    _aff->fwd[1][0] = dy1;
    _aff->fwd[1][1] = dy2;
    _aff->inv[0][0] = QR_DIVROUND(dy2 << _res, det >> ires);
    _aff->inv[0][1] = QR_DIVROUND(-dx2 << _res, det >> ires);
    _aff->inv[1][0] = QR_DIVROUND(-dy1 << _res, det >> ires);
    _aff->inv[1][1] = QR_DIVROUND(dx1 << _res, det >> ires);
    _aff->x0 = _p0[0];
    _aff->y0 = _p0[1];
    _aff->res = _res;
    _aff->ires = ires;
}


static void qr_aff_unproject(qr_point _q, const qr_aff* _aff,
    int _x, int _y) {
    _q[0] = _aff->inv[0][0] * (_x - _aff->x0) + _aff->inv[0][1] * (_y - _aff->y0)
        + (1 << _aff->ires >> 1) >> _aff->ires;
    _q[1] = _aff->inv[1][0] * (_x - _aff->x0) + _aff->inv[1][1] * (_y - _aff->y0)
        + (1 << _aff->ires >> 1) >> _aff->ires;
}

static void qr_point_translate(qr_point _point, int _dx, int _dy) {
    _point[0] += _dx;
    _point[1] += _dy;
}

static int qr_cmp_edge_pt(const void* _a, const void* _b) {
    const qr_finder_edge_pt* a;
    const qr_finder_edge_pt* b;
    a = (const qr_finder_edge_pt*)_a;
    b = (const qr_finder_edge_pt*)_b;
    return ((a->edge > b->edge) - (a->edge < b->edge) << 1) +
        (a->extent > b->extent) - (a->extent < b->extent);
}

static void qr_finder_edge_pts_aff_classify(qr_finder* _f, const qr_aff* _aff) {
    qr_finder_center* c;
    int               i;
    int               e;
    c = _f->c;
    for (e = 0; e < 4; e++)_f->nedge_pts[e] = 0;
    for (i = 0; i < c->nedge_pts; i++) {
        qr_point q;
        int      d;
        qr_aff_unproject(q, _aff, c->edge_pts[i].pos[0], c->edge_pts[i].pos[1]);
        qr_point_translate(q, -_f->o[0], -_f->o[1]);
        d = abs(q[1]) > abs(q[0]);
        e = d << 1 | (q[d] >= 0);
        _f->nedge_pts[e]++;
        c->edge_pts[i].edge = e;
        c->edge_pts[i].extent = q[d];
    }
    qsort(c->edge_pts, c->nedge_pts, sizeof(*c->edge_pts), qr_cmp_edge_pt);
    _f->edge_pts[0] = c->edge_pts;
    for (e = 1; e < 4; e++)_f->edge_pts[e] = _f->edge_pts[e - 1] + _f->nedge_pts[e - 1];
}

static int qr_finder_estimate_module_size_and_version(qr_finder* _f,
    int _width, int _height) {
    qr_point offs;
    int      sums[4];
    int      nsums[4];
    int      usize;
    int      nusize;
    int      vsize;
    int      nvsize;
    int      uversion;
    int      vversion;
    int      e;
    offs[0] = offs[1] = 0;
    for (e = 0; e < 4; e++)if (_f->nedge_pts[e] > 0) {
        qr_finder_edge_pt* edge_pts;
        int                sum;
        int                mean;
        int                n;
        int                i;
        edge_pts = _f->edge_pts[e];
        n = _f->nedge_pts[e];
        sum = 0;
        for (i = (n >> 2); i < n - (n >> 2); i++)sum += edge_pts[i].extent;
        n = n - ((n >> 2) << 1);
        mean = QR_DIVROUND(sum, n);
        offs[e >> 1] += mean;
        sums[e] = sum;
        nsums[e] = n;
    }
    else nsums[e] = sums[e] = 0;
    if (_f->nedge_pts[0] > 0 && _f->nedge_pts[1] > 0) {
        _f->o[0] -= offs[0] >> 1;
        sums[0] -= offs[0] * nsums[0] >> 1;
        sums[1] -= offs[0] * nsums[1] >> 1;
    }
    if (_f->nedge_pts[2] > 0 && _f->nedge_pts[3] > 0) {
        _f->o[1] -= offs[1] >> 1;
        sums[2] -= offs[1] * nsums[2] >> 1;
        sums[3] -= offs[1] * nsums[3] >> 1;
    }
    nusize = nsums[0] + nsums[1];
    if (nusize <= 0)return -1;
    nusize *= 3;
    usize = sums[1] - sums[0];
    usize = ((usize << 1) + nusize) / (nusize << 1);
    if (usize <= 0)return -1;
    uversion = (_width - 8 * usize) / (usize << 2);
    if (uversion < 1 || uversion>40 + QR_LARGE_VERSION_SLACK)return -1;
    /*Now do the same for the other axis.*/
    nvsize = nsums[2] + nsums[3];
    if (nvsize <= 0)return -1;
    nvsize *= 3;
    vsize = sums[3] - sums[2];
    vsize = ((vsize << 1) + nvsize) / (nvsize << 1);
    if (vsize <= 0)return -1;
    vversion = (_height - 8 * vsize) / (vsize << 2);
    if (vversion < 1 || vversion>40 + QR_LARGE_VERSION_SLACK)return -1;
    if (abs(uversion - vversion) > QR_LARGE_VERSION_SLACK)return -1;
    _f->size[0] = usize;
    _f->size[1] = vsize;
    _f->eversion[0] = uversion;
    _f->eversion[1] = vversion;
    return 0;
}


unsigned isaac_next_uint32(isaac_ctx* _ctx) {
    if (!_ctx->n)isaac_update(_ctx);
    return _ctx->r[--_ctx->n];
}

unsigned isaac_next_uint(isaac_ctx* _ctx, unsigned _n) {
    unsigned r;
    unsigned v;
    unsigned d;
    do {
        r = isaac_next_uint32(_ctx);
        v = r % _n;
        d = r - v;
    } while ((d + _n - 1 & ISAAC_MASK) < d);
    return v;
}

unsigned qr_isqrt(unsigned _val) {
    unsigned g;
    unsigned b;
    int      bshift;
    g = 0;
    b = 0x8000;
    for (bshift = 16; bshift-- > 0;) {
        unsigned t;
        t = (g << 1) + b << bshift;
        if (t <= _val) {
            g += b;
            _val -= t;
        }
        b >>= 1;
    }
    return g;
}

static void qr_finder_ransac(qr_finder* _f, const qr_aff* _hom,
    isaac_ctx* _isaac, int _e) {
    qr_finder_edge_pt* edge_pts;
    int                best_ninliers;
    int                n;
    edge_pts = _f->edge_pts[_e];
    n = _f->nedge_pts[_e];
    best_ninliers = 0;
    if (n > 1) {
        int max_iters;
        int i;
        int j;
        max_iters = 17;
        for (i = 0; i < max_iters; i++) {
            qr_point  q0;
            qr_point  q1;
            int       ninliers;
            int       thresh;
            int       p0i;
            int       p1i;
            int* p0;
            int* p1;
            int       j;
            p0i = isaac_next_uint(_isaac, n);
            p1i = isaac_next_uint(_isaac, n - 1);
            if (p1i >= p0i)p1i++;
            p0 = edge_pts[p0i].pos;
            p1 = edge_pts[p1i].pos;
            qr_aff_unproject(q0, _hom, p0[0], p0[1]);
            qr_aff_unproject(q1, _hom, p1[0], p1[1]);
            qr_point_translate(q0, -_f->o[0], -_f->o[1]);
            qr_point_translate(q1, -_f->o[0], -_f->o[1]);
            if (abs(q0[_e >> 1] - q1[_e >> 1]) > abs(q0[1 - (_e >> 1)] - q1[1 - (_e >> 1)]))continue;
            thresh = qr_isqrt(qr_point_distance2(p0, p1) << 2 * QR_FINDER_SUBPREC + 1);
            ninliers = 0;
            for (j = 0; j < n; j++) {
                if (abs(qr_point_ccw(p0, p1, edge_pts[j].pos)) <= thresh) {
                    edge_pts[j].extent |= 1;
                    ninliers++;
                }
                else edge_pts[j].extent &= ~1;
            }
            if (ninliers > best_ninliers) {
                for (j = 0; j < n; j++)edge_pts[j].extent <<= 1;
                best_ninliers = ninliers;
                /*The actual number of iterations required is
                    log(1-\alpha)/log(1-r*r),
                   where \alpha is the required probability of taking a sample with
                    no outliers (e.g., 0.99) and r is the estimated ratio of inliers
                    (e.g. ninliers/n).
                  This is just a rough (but conservative) approximation, but it
                   should be good enough to stop the iteration early when we find
                   a good set of inliers.*/
                if (ninliers > n >> 1)max_iters = (67 * n - 63 * ninliers - 1) / (n << 1);
            }
        }
        /*Now collect all the inliers at the beginning of the list.*/
        for (i = j = 0; j < best_ninliers; i++)if (edge_pts[i].extent & 2) {
            if (j < i) {
                qr_finder_edge_pt tmp;
                *&tmp = *(edge_pts + i);
                *(edge_pts + j) = *(edge_pts + i);
                *(edge_pts + i) = *&tmp;
            }
            j++;
        }
    }
    _f->ninliers[_e] = best_ninliers;
}

static void qr_aff_project(qr_point _p, const qr_aff* _aff,
    int _u, int _v) {
    _p[0] = (_aff->fwd[0][0] * _u + _aff->fwd[0][1] * _v + (1 << _aff->res - 1) >> _aff->res)
        + _aff->x0;
    _p[1] = (_aff->fwd[1][0] * _u + _aff->fwd[1][1] * _v + (1 << _aff->res - 1) >> _aff->res)
        + _aff->y0;
}

unsigned qr_ihypot(int _x, int _y) {
    unsigned x;
    unsigned y;
    int      mask;
    int      shift;
    int      u;
    int      v;
    int      i;
    x = _x = abs(_x);
    y = _y = abs(_y);
    mask = -(x > y) & (_x ^ _y);
    x ^= mask;
    y ^= mask;
    _y ^= mask;
    shift = 31 - qr_ilog(y);
    shift = QR_MAXI(shift, 0);
    x = (unsigned)((x << shift) * 0x9B74EDAAULL >> 32);
    _y = (int)((_y << shift) * 0x9B74EDA9LL >> 32);
    u = x;
    mask = -(_y < 0);
    x += _y + mask ^ mask;
    _y -= u + mask ^ mask;
    u = x + 1 >> 1;
    v = _y + 1 >> 1;
    mask = -(_y < 0);
    x += v + mask ^ mask;
    _y -= u + mask ^ mask;
    for (i = 1; i < 16; i++) {
        int r;
        u = x + 1 >> 2;
        r = (1 << 2 * i) >> 1;
        v = _y + r >> 2 * i;
        mask = -(_y < 0);
        x += v + mask ^ mask;
        _y = _y - (u + mask ^ mask) << 1;
    }
    return x + ((1U << shift) >> 1) >> shift;
}

static void qr_line_fit(qr_line _l, int _x0, int _y0,
    int _sxx, int _sxy, int _syy, int _res) {
    int dshift;
    int dround;
    int u;
    int v;
    int w;
    u = abs(_sxx - _syy);
    v = -_sxy << 1;
    w = qr_ihypot(u, v);
    dshift = QR_MAXI(0, QR_MAXI(qr_ilog(u), qr_ilog(abs(v))) + 1 - (_res + 1 >> 1));
    dround = (1 << dshift) >> 1;
    if (_sxx > _syy) {
        _l[0] = v + dround >> dshift;
        _l[1] = u + w + dround >> dshift;
    }
    else {
        _l[0] = u + w + dround >> dshift;
        _l[1] = v + dround >> dshift;
    }
    _l[2] = -(_x0 * _l[0] + _y0 * _l[1]);
}

static void qr_line_fit_points(qr_line _l, qr_point* _p, int _np, int _res) {
    int sx;
    int sy;
    int xmin;
    int xmax;
    int ymin;
    int ymax;
    int xbar;
    int ybar;
    int dx;
    int dy;
    int sxx;
    int sxy;
    int syy;
    int sshift;
    int sround;
    int i;
    sx = sy = 0;
    ymax = xmax = INT_MIN;
    ymin = xmin = INT_MAX;
    for (i = 0; i < _np; i++) {
        sx += _p[i][0];
        xmin = QR_MINI(xmin, _p[i][0]);
        xmax = QR_MAXI(xmax, _p[i][0]);
        sy += _p[i][1];
        ymin = QR_MINI(ymin, _p[i][1]);
        ymax = QR_MAXI(ymax, _p[i][1]);
    }
    xbar = (sx + (_np >> 1)) / _np;
    ybar = (sy + (_np >> 1)) / _np;
    sshift = QR_MAXI(0, qr_ilog(_np * QR_MAXI(QR_MAXI(xmax - xbar, xbar - xmin),
        QR_MAXI(ymax - ybar, ybar - ymin))) - (QR_INT_BITS - 1 >> 1));
    sround = (1 << sshift) >> 1;
    sxx = sxy = syy = 0;
    for (i = 0; i < _np; i++) {
        dx = _p[i][0] - xbar + sround >> sshift;
        dy = _p[i][1] - ybar + sround >> sshift;
        sxx += dx * dx;
        sxy += dx * dy;
        syy += dy * dy;
    }
    qr_line_fit(_l, xbar, ybar, sxx, sxy, syy, _res);
}

static int qr_line_eval(qr_line _line, int _x, int _y) {
    return _line[0] * _x + _line[1] * _y + _line[2];
}

static void qr_line_orient(qr_line _l, int _x, int _y) {
    if (qr_line_eval(_l, _x, _y) < 0) {
        _l[0] = -_l[0];
        _l[1] = -_l[1];
        _l[2] = -_l[2];
    }
}

static void qr_line_fit_finder_pair(qr_line _l, const qr_aff* _aff,
    const qr_finder* _f0, const qr_finder* _f1, int _e) {
    qr_point* pts;
    int                npts;
    qr_finder_edge_pt* edge_pts;
    qr_point           q;
    int                n0;
    int                n1;
    int                i;
    n0 = _f0->ninliers[_e];
    n1 = _f1->ninliers[_e];
    npts = QR_MAXI(n0, 1) + QR_MAXI(n1, 1);
    pts = (qr_point*)malloc(npts * sizeof(*pts));
    if (n0 > 0) {
        edge_pts = _f0->edge_pts[_e];
        for (i = 0; i < n0; i++) {
            pts[i][0] = edge_pts[i].pos[0];
            pts[i][1] = edge_pts[i].pos[1];
        }
    }
    else {
        q[0] = _f0->o[0];
        q[1] = _f0->o[1];
        q[_e >> 1] += _f0->size[_e >> 1] * (2 * (_e & 1) - 1);
        qr_aff_project(pts[0], _aff, q[0], q[1]);
        n0++;
    }
    if (n1 > 0) {
        edge_pts = _f1->edge_pts[_e];
        for (i = 0; i < n1; i++) {
            pts[n0 + i][0] = edge_pts[i].pos[0];
            pts[n0 + i][1] = edge_pts[i].pos[1];
        }
    }
    else {
        q[0] = _f1->o[0];
        q[1] = _f1->o[1];
        q[_e >> 1] += _f1->size[_e >> 1] * (2 * (_e & 1) - 1);
        qr_aff_project(pts[n0], _aff, q[0], q[1]);
        n1++;
    }
    qr_line_fit_points(_l, pts, npts, _aff->res);
    qr_line_orient(_l, _f0->c->pos[0], _f0->c->pos[1]);
    free(pts);
}

static int qr_line_fit_finder_edge(qr_line _l,
    const qr_finder* _f, int _e, int _res) {
    qr_finder_edge_pt* edge_pts;
    qr_point* pts;
    int                npts;
    int                i;
    npts = _f->ninliers[_e];
    if (npts < 2)return -1;
    pts = (qr_point*)malloc(npts * sizeof(*pts));
    edge_pts = _f->edge_pts[_e];
    for (i = 0; i < npts; i++) {
        pts[i][0] = edge_pts[i].pos[0];
        pts[i][1] = edge_pts[i].pos[1];
    }
    qr_line_fit_points(_l, pts, npts, _res);
    qr_line_orient(_l, _f->c->pos[0], _f->c->pos[1]);
    free(pts);
    return 0;
}

static int qr_aff_line_step(const qr_aff* _aff, qr_line _l,
    int _v, int _du, int* _dv) {
    int shift;
    int round;
    int dv;
    int n;
    int d;
    n = _aff->fwd[0][_v] * _l[0] + _aff->fwd[1][_v] * _l[1];
    d = _aff->fwd[0][1 - _v] * _l[0] + _aff->fwd[1][1 - _v] * _l[1];
    if (d < 0) {
        n = -n;
        d = -d;
    }
    shift = QR_MAXI(0, qr_ilog(_du) + qr_ilog(abs(n)) + 3 - QR_INT_BITS);
    round = (1 << shift) >> 1;
    n = n + round >> shift;
    d = d + round >> shift;
    if (abs(n) >= d)return -1;
    n = -_du * n;
    dv = QR_DIVROUND(n, d);
    if (abs(dv) >= _du)return -1;
    *_dv = dv;
    return 0;
}

static int qr_finder_quick_crossing_check(const unsigned char* _img,
    int _width, int _height, int _x0, int _y0, int _x1, int _y1, int _v) {
    if (_x0 < 0 || _x0 >= _width || _y0 < 0 || _y0 >= _height ||
        _x1 < 0 || _x1 >= _width || _y1 < 0 || _y1 >= _height) {
        return -1;
    }
    if (!_img[_y0 * _width + _x0] != _v || !_img[_y1 * _width + _x1] != _v)return 1;
    if (!_img[(_y0 + _y1 >> 1) * _width + (_x0 + _x1 >> 1)] == _v)return -1;
    return 0;
}

static int qr_finder_locate_crossing(const unsigned char* _img,
    int _width, int _height, int _x0, int _y0, int _x1, int _y1, int _v, qr_point _p) {
    qr_point x0;
    qr_point x1;
    qr_point dx;
    int      step[2];
    int      steep;
    int      err;
    int      derr;
    x0[0] = _x0;
    x0[1] = _y0;
    x1[0] = _x1;
    x1[1] = _y1;
    dx[0] = abs(_x1 - _x0);
    dx[1] = abs(_y1 - _y0);
    steep = dx[1] > dx[0];
    err = 0;
    derr = dx[1 - steep];
    step[0] = ((_x0 < _x1) << 1) - 1;
    step[1] = ((_y0 < _y1) << 1) - 1;
    for (;;) {
        if (x0[steep] == x1[steep])return -1;
        x0[steep] += step[steep];
        err += derr;
        if (err << 1 > dx[steep]) {
            x0[1 - steep] += step[1 - steep];
            err -= dx[steep];
        }
        if (!_img[x0[1] * _width + x0[0]] != _v)break;
    }
    err = 0;
    for (;;) {
        if (x0[steep] == x1[steep])break;
        x1[steep] -= step[steep];
        err += derr;
        if (err << 1 > dx[steep]) {
            x1[1 - steep] -= step[1 - steep];
            err -= dx[steep];
        }
        if (!_img[x1[1] * _width + x1[0]] != _v)break;
    }
    _p[0] = (x0[0] + x1[0] + 1 << QR_FINDER_SUBPREC) >> 1;
    _p[1] = (x0[1] + x1[1] + 1 << QR_FINDER_SUBPREC) >> 1;
    return 0;
}

static int qr_line_isect(qr_point _p, const qr_line _l0, const qr_line _l1) {
    int d;
    int x;
    int y;
    d = _l0[0] * _l1[1] - _l0[1] * _l1[0];
    if (d == 0)return -1;
    x = _l0[1] * _l1[2] - _l1[1] * _l0[2];
    y = _l1[0] * _l0[2] - _l0[0] * _l1[2];
    if (d < 0) {
        x = -x;
        y = -y;
        d = -d;
    }
    _p[0] = QR_DIVROUND(x, d);
    _p[1] = QR_DIVROUND(y, d);
    return 0;
}

static void qr_hom_cell_init(qr_hom_cell* _cell, int _u0, int _v0,
    int _u1, int _v1, int _u2, int _v2, int _u3, int _v3, int _x0, int _y0,
    int _x1, int _y1, int _x2, int _y2, int _x3, int _y3) {
    int du10;
    int du20;
    int du30;
    int du31;
    int du32;
    int dv10;
    int dv20;
    int dv30;
    int dv31;
    int dv32;
    int dx10;
    int dx20;
    int dx30;
    int dx31;
    int dx32;
    int dy10;
    int dy20;
    int dy30;
    int dy31;
    int dy32;
    int a00;
    int a01;
    int a02;
    int a10;
    int a11;
    int a12;
    int a20;
    int a21;
    int a22;
    int i00;
    int i01;
    int i10;
    int i11;
    int i20;
    int i21;
    int i22;
    int b0;
    int b1;
    int b2;
    int shift;
    int round;
    int x;
    int y;
    int w;
    du10 = _u1 - _u0;
    du20 = _u2 - _u0;
    du30 = _u3 - _u0;
    du31 = _u3 - _u1;
    du32 = _u3 - _u2;
    dv10 = _v1 - _v0;
    dv20 = _v2 - _v0;
    dv30 = _v3 - _v0;
    dv31 = _v3 - _v1;
    dv32 = _v3 - _v2;
    a20 = du32 * dv10 - du10 * dv32;
    a21 = du20 * dv31 - du31 * dv20;
    if (a20 || a21)a22 = du32 * dv31 - du31 * dv32;
    else a22 = 1;
    a00 = du10 * (a20 + a22);
    a01 = du20 * (a21 + a22);
    a10 = dv10 * (a20 + a22);
    a11 = dv20 * (a21 + a22);
    i00 = a11 * a22;
    i01 = -a01 * a22;
    i10 = -a10 * a22;
    i11 = a00 * a22;
    i20 = a10 * a21 - a11 * a20;
    i21 = a01 * a20 - a00 * a21;
    i22 = a00 * a11 - a01 * a10;
    if (i00)i00 = QR_FLIPSIGNI(QR_DIVROUND(i22, abs(i00)), i00);
    if (i01)i01 = QR_FLIPSIGNI(QR_DIVROUND(i22, abs(i01)), i01);
    if (i10)i10 = QR_FLIPSIGNI(QR_DIVROUND(i22, abs(i10)), i10);
    if (i11)i11 = QR_FLIPSIGNI(QR_DIVROUND(i22, abs(i11)), i11);
    if (i20)i20 = QR_FLIPSIGNI(QR_DIVROUND(i22, abs(i20)), i20);
    if (i21)i21 = QR_FLIPSIGNI(QR_DIVROUND(i22, abs(i21)), i21);
    dx10 = _x1 - _x0;
    dx20 = _x2 - _x0;
    dx30 = _x3 - _x0;
    dx31 = _x3 - _x1;
    dx32 = _x3 - _x2;
    dy10 = _y1 - _y0;
    dy20 = _y2 - _y0;
    dy30 = _y3 - _y0;
    dy31 = _y3 - _y1;
    dy32 = _y3 - _y2;
    a20 = dx32 * dy10 - dx10 * dy32;
    a21 = dx20 * dy31 - dx31 * dy20;
    a22 = dx32 * dy31 - dx31 * dy32;
    b0 = qr_ilog(QR_MAXI(abs(dx10), abs(dy10))) + qr_ilog(abs(a20 + a22));
    b1 = qr_ilog(QR_MAXI(abs(dx20), abs(dy20))) + qr_ilog(abs(a21 + a22));
    b2 = qr_ilog(QR_MAXI(QR_MAXI(abs(a20), abs(a21)), abs(a22)));
    shift = QR_MAXI(0, QR_MAXI(QR_MAXI(b0, b1), b2) - (QR_INT_BITS - 3 - QR_ALIGN_SUBPREC));
    round = (1 << shift) >> 1;
    a00 = QR_FIXMUL(dx10, a20 + a22, round, shift);
    a01 = QR_FIXMUL(dx20, a21 + a22, round, shift);
    a10 = QR_FIXMUL(dy10, a20 + a22, round, shift);
    a11 = QR_FIXMUL(dy20, a21 + a22, round, shift);
    _cell->fwd[0][0] = (i00 ? QR_DIVROUND(a00, i00) : 0) + (i10 ? QR_DIVROUND(a01, i10) : 0);
    _cell->fwd[0][1] = (i01 ? QR_DIVROUND(a00, i01) : 0) + (i11 ? QR_DIVROUND(a01, i11) : 0);
    _cell->fwd[1][0] = (i00 ? QR_DIVROUND(a10, i00) : 0) + (i10 ? QR_DIVROUND(a11, i10) : 0);
    _cell->fwd[1][1] = (i01 ? QR_DIVROUND(a10, i01) : 0) + (i11 ? QR_DIVROUND(a11, i11) : 0);
    _cell->fwd[2][0] = (i00 ? QR_DIVROUND(a20, i00) : 0) + (i10 ? QR_DIVROUND(a21, i10) : 0)
        + (i20 ? QR_DIVROUND(a22, i20) : 0) + round >> shift;
    _cell->fwd[2][1] = (i01 ? QR_DIVROUND(a20, i01) : 0) + (i11 ? QR_DIVROUND(a21, i11) : 0)
        + (i21 ? QR_DIVROUND(a22, i21) : 0) + round >> shift;
    _cell->fwd[2][2] = a22 + round >> shift;
    x = _cell->fwd[0][0] * du10 + _cell->fwd[0][1] * dv10;
    y = _cell->fwd[1][0] * du10 + _cell->fwd[1][1] * dv10;
    w = _cell->fwd[2][0] * du10 + _cell->fwd[2][1] * dv10 + _cell->fwd[2][2];
    a02 = dx10 * w - x;
    a12 = dy10 * w - y;
    x = _cell->fwd[0][0] * du20 + _cell->fwd[0][1] * dv20;
    y = _cell->fwd[1][0] * du20 + _cell->fwd[1][1] * dv20;
    w = _cell->fwd[2][0] * du20 + _cell->fwd[2][1] * dv20 + _cell->fwd[2][2];
    a02 += dx20 * w - x;
    a12 += dy20 * w - y;
    x = _cell->fwd[0][0] * du30 + _cell->fwd[0][1] * dv30;
    y = _cell->fwd[1][0] * du30 + _cell->fwd[1][1] * dv30;
    w = _cell->fwd[2][0] * du30 + _cell->fwd[2][1] * dv30 + _cell->fwd[2][2];
    a02 += dx30 * w - x;
    a12 += dy30 * w - y;
    _cell->fwd[0][2] = a02 + 2 >> 2;
    _cell->fwd[1][2] = a12 + 2 >> 2;
    _cell->x0 = _x0;
    _cell->y0 = _y0;
    _cell->u0 = _u0;
    _cell->v0 = _v0;
}

static void qr_hom_cell_fproject(qr_point _p, const qr_hom_cell* _cell,
    int _x, int _y, int _w) {
    if (_w == 0) {
        _p[0] = _x < 0 ? INT_MIN : INT_MAX;
        _p[1] = _y < 0 ? INT_MIN : INT_MAX;
    }
    else {
        if (_w < 0) {
            _x = -_x;
            _y = -_y;
            _w = -_w;
        }
        _p[0] = QR_DIVROUND(_x, _w) + _cell->x0;
        _p[1] = QR_DIVROUND(_y, _w) + _cell->y0;
    }
}

static int qr_img_get_bit(const unsigned char* _img, int _width, int _height,
    int _x, int _y) {
    _x >>= QR_FINDER_SUBPREC;
    _y >>= QR_FINDER_SUBPREC;
    return _img[QR_CLAMPI(0, _y, _height - 1) * _width + QR_CLAMPI(0, _x, _width - 1)] != 0;
}

static unsigned qr_alignment_pattern_fetch(qr_point _p[5][5], int _x0, int _y0,
    const unsigned char* _img, int _width, int _height) {
    unsigned v;
    int      i;
    int      j;
    int      k;
    int      dx;
    int      dy;
    dx = _x0 - _p[2][2][0];
    dy = _y0 - _p[2][2][1];
    v = 0;
    for (k = i = 0; i < 5; i++)for (j = 0; j < 5; j++, k++) {
        v |= qr_img_get_bit(_img, _width, _height, _p[i][j][0] + dx, _p[i][j][1] + dy) << k;
    }
    return v;
}

static int qr_hamming_dist(unsigned _y1, unsigned _y2, int _maxdiff) {
    unsigned y;
    int      ret;
    y = _y1 ^ _y2;
    for (ret = 0; ret < _maxdiff && y; ret++)y &= y - 1;
    return ret;
}

static int qr_alignment_pattern_search(qr_point _p, const qr_hom_cell* _cell,
    int _u, int _v, int _r, const unsigned char* _img, int _width, int _height) {
    qr_point c[4];
    int      nc[4];
    qr_point p[5][5];
    qr_point pc;
    unsigned best_match;
    int      best_dist;
    int      bestx;
    int      besty;
    unsigned match;
    int      dist;
    int      u;
    int      v;
    int      x0;
    int      y0;
    int      w0;
    int      x;
    int      y;
    int      w;
    int      dxdu;
    int      dydu;
    int      dwdu;
    int      dxdv;
    int      dydv;
    int      dwdv;
    int      dx;
    int      dy;
    int      i;
    int      j;
    u = (_u - 2) - _cell->u0;
    v = (_v - 2) - _cell->v0;
    x0 = _cell->fwd[0][0] * u + _cell->fwd[0][1] * v + _cell->fwd[0][2];
    y0 = _cell->fwd[1][0] * u + _cell->fwd[1][1] * v + _cell->fwd[1][2];
    w0 = _cell->fwd[2][0] * u + _cell->fwd[2][1] * v + _cell->fwd[2][2];
    dxdu = _cell->fwd[0][0];
    dydu = _cell->fwd[1][0];
    dwdu = _cell->fwd[2][0];
    dxdv = _cell->fwd[0][1];
    dydv = _cell->fwd[1][1];
    dwdv = _cell->fwd[2][1];
    for (i = 0; i < 5; i++) {
        x = x0;
        y = y0;
        w = w0;
        for (j = 0; j < 5; j++) {
            qr_hom_cell_fproject(p[i][j], _cell, x, y, w);
            x += dxdu;
            y += dydu;
            w += dwdu;
        }
        x0 += dxdv;
        y0 += dydv;
        w0 += dwdv;
    }
    bestx = p[2][2][0];
    besty = p[2][2][1];
    best_match = qr_alignment_pattern_fetch(p, bestx, besty, _img, _width, _height);
    best_dist = qr_hamming_dist(best_match, 0x1F8D63F, 25);
    if (best_dist > 0) {
        u = _u - _cell->u0;
        v = _v - _cell->v0;
        x = _cell->fwd[0][0] * u + _cell->fwd[0][1] * v + _cell->fwd[0][2] << QR_ALIGN_SUBPREC;
        y = _cell->fwd[1][0] * u + _cell->fwd[1][1] * v + _cell->fwd[1][2] << QR_ALIGN_SUBPREC;
        w = _cell->fwd[2][0] * u + _cell->fwd[2][1] * v + _cell->fwd[2][2] << QR_ALIGN_SUBPREC;
        /*Search an area at most _r modules around the target location, in
           concentric squares..*/
        for (i = 1; i < _r << QR_ALIGN_SUBPREC; i++) {
            int side_len;
            side_len = (i << 1) - 1;
            x -= dxdu + dxdv;
            y -= dydu + dydv;
            w -= dwdu + dwdv;
            for (j = 0; j < 4 * side_len; j++) {
                int      dir;
                qr_hom_cell_fproject(pc, _cell, x, y, w);
                match = qr_alignment_pattern_fetch(p, pc[0], pc[1], _img, _width, _height);
                dist = qr_hamming_dist(match, 0x1F8D63F, best_dist + 1);
                if (dist < best_dist) {
                    best_match = match;
                    best_dist = dist;
                    bestx = pc[0];
                    besty = pc[1];
                }
                if (j < 2 * side_len) {
                    dir = j >= side_len;
                    x += _cell->fwd[0][dir];
                    y += _cell->fwd[1][dir];
                    w += _cell->fwd[2][dir];
                }
                else {
                    dir = j >= 3 * side_len;
                    x -= _cell->fwd[0][dir];
                    y -= _cell->fwd[1][dir];
                    w -= _cell->fwd[2][dir];
                }
                if (!best_dist)break;
            }
            if (!best_dist)break;
        }
    }
    if (best_dist > 6) {
        _p[0] = p[2][2][0];
        _p[1] = p[2][2][1];
        return -1;
    }
    dx = bestx - p[2][2][0];
    dy = besty - p[2][2][1];
    memset(nc, 0, sizeof(nc));
    memset(c, 0, sizeof(c));
    for (i = 0; i < 8; i++) {
        static const unsigned MASK_TESTS[8][2] = {
          {0x1040041,0x1000001},{0x0041040,0x0001000},
          {0x0110110,0x0100010},{0x0011100,0x0001000},
          {0x0420084,0x0400004},{0x0021080,0x0001000},
          {0x0006C00,0x0004400},{0x0003800,0x0001000},
        };
        static const unsigned char MASK_COORDS[8][2] = {
          {0,0},{1,1},{4,0},{3,1},{2,0},{2,1},{0,2},{1,2}
        };
        if ((best_match & MASK_TESTS[i][0]) == MASK_TESTS[i][1]) {
            int x0;
            int y0;
            int x1;
            int y1;
            x0 = p[MASK_COORDS[i][1]][MASK_COORDS[i][0]][0] + dx >> QR_FINDER_SUBPREC;
            if (x0 < 0 || x0 >= _width)continue;
            y0 = p[MASK_COORDS[i][1]][MASK_COORDS[i][0]][1] + dy >> QR_FINDER_SUBPREC;
            if (y0 < 0 || y0 >= _height)continue;
            x1 = p[4 - MASK_COORDS[i][1]][4 - MASK_COORDS[i][0]][0] + dx >> QR_FINDER_SUBPREC;
            if (x1 < 0 || x1 >= _width)continue;
            y1 = p[4 - MASK_COORDS[i][1]][4 - MASK_COORDS[i][0]][1] + dy >> QR_FINDER_SUBPREC;
            if (y1 < 0 || y1 >= _height)continue;
            if (!qr_finder_locate_crossing(_img, _width, _height, x0, y0, x1, y1, i & 1, pc)) {
                int w;
                int cx;
                int cy;
                cx = pc[0] - bestx;
                cy = pc[1] - besty;
                if (i & 1) {
                    /*Weight crossings around the center dot more highly, as they are
                       generally more reliable.*/
                    w = 3;
                    cx += cx << 1;
                    cy += cy << 1;
                }
                else w = 1;
                nc[i >> 1] += w;
                c[i >> 1][0] += cx;
                c[i >> 1][1] += cy;
            }
        }
    }
    /*Sum offsets from lines in orthogonal directions.*/
    for (i = 0; i < 2; i++) {
        int a;
        int b;
        a = nc[i << 1];
        b = nc[i << 1 | 1];
        if (a && b) {
            int w;
            w = QR_MAXI(a, b);
            c[i << 1][0] = QR_DIVROUND(w * (b * c[i << 1][0] + a * c[i << 1 | 1][0]), a * b);
            c[i << 1][1] = QR_DIVROUND(w * (b * c[i << 1][1] + a * c[i << 1 | 1][1]), a * b);
            nc[i << 1] = w << 1;
        }
        else {
            c[i << 1][0] += c[i << 1 | 1][0];
            c[i << 1][1] += c[i << 1 | 1][1];
            nc[i << 1] += b;
        }
    }
    /*Average offsets from pairs of orthogonal lines.*/
    c[0][0] += c[2][0];
    c[0][1] += c[2][1];
    nc[0] += nc[2];
    /*If we actually found any such lines, apply the adjustment.*/
    if (nc[0]) {
        dx = QR_DIVROUND(c[0][0], nc[0]);
        dy = QR_DIVROUND(c[0][1], nc[0]);
        /*But only if it doesn't make things too much worse.*/
        match = qr_alignment_pattern_fetch(p, bestx + dx, besty + dy, _img, _width, _height);
        dist = qr_hamming_dist(match, 0x1F8D63F, best_dist + 1);
        if (dist <= best_dist + 1) {
            bestx += dx;
            besty += dy;
        }
    }
    _p[0] = bestx;
    _p[1] = besty;
    return 0;
}

static void qr_hom_init(qr_hom* _hom, int _x0, int _y0,
    int _x1, int _y1, int _x2, int _y2, int _x3, int _y3, int _res) {
    int dx10;
    int dx20;
    int dx30;
    int dx31;
    int dx32;
    int dy10;
    int dy20;
    int dy30;
    int dy31;
    int dy32;
    int a20;
    int a21;
    int a22;
    int b0;
    int b1;
    int b2;
    int s1;
    int s2;
    int r1;
    int r2;
    dx10 = _x1 - _x0;
    dx20 = _x2 - _x0;
    dx30 = _x3 - _x0;
    dx31 = _x3 - _x1;
    dx32 = _x3 - _x2;
    dy10 = _y1 - _y0;
    dy20 = _y2 - _y0;
    dy30 = _y3 - _y0;
    dy31 = _y3 - _y1;
    dy32 = _y3 - _y2;
    a20 = dx32 * dy10 - dx10 * dy32;
    a21 = dx20 * dy31 - dx31 * dy20;
    a22 = dx32 * dy31 - dx31 * dy32;
    b0 = qr_ilog(QR_MAXI(abs(dx10), abs(dy10))) + qr_ilog(abs(a20 + a22));
    b1 = qr_ilog(QR_MAXI(abs(dx20), abs(dy20))) + qr_ilog(abs(a21 + a22));
    b2 = qr_ilog(QR_MAXI(QR_MAXI(abs(a20), abs(a21)), abs(a22)));
    s1 = QR_MAXI(0, _res + QR_MAXI(QR_MAXI(b0, b1), b2) - (QR_INT_BITS - 2));
    r1 = (1 << s1) >> 1;
    _hom->fwd[0][0] = QR_FIXMUL(dx10, a20 + a22, r1, s1);
    _hom->fwd[0][1] = QR_FIXMUL(dx20, a21 + a22, r1, s1);
    _hom->x0 = _x0;
    _hom->fwd[1][0] = QR_FIXMUL(dy10, a20 + a22, r1, s1);
    _hom->fwd[1][1] = QR_FIXMUL(dy20, a21 + a22, r1, s1);
    _hom->y0 = _y0;
    _hom->fwd[2][0] = a20 + r1 >> s1;
    _hom->fwd[2][1] = a21 + r1 >> s1;
    _hom->fwd22 = s1 > _res ? a22 + (r1 >> _res) >> s1 - _res : a22 << _res - s1;
    b0 = qr_ilog(QR_MAXI(QR_MAXI(abs(dx10), abs(dx20)), abs(dx30))) +
        qr_ilog(QR_MAXI(abs(_hom->fwd[0][0]), abs(_hom->fwd[1][0])));
    b1 = qr_ilog(QR_MAXI(QR_MAXI(abs(dy10), abs(dy20)), abs(dy30))) +
        qr_ilog(QR_MAXI(abs(_hom->fwd[0][1]), abs(_hom->fwd[1][1])));
    b2 = qr_ilog(abs(a22)) - s1;
    s2 = QR_MAXI(0, QR_MAXI(b0, b1) + b2 - (QR_INT_BITS - 3));
    r2 = (1 << s2) >> 1;
    s1 += s2;
    r1 <<= s2;
    _hom->inv[0][0] = QR_FIXMUL(_hom->fwd[1][1], a22, r1, s1);
    _hom->inv[0][1] = QR_FIXMUL(-_hom->fwd[0][1], a22, r1, s1);
    _hom->inv[1][0] = QR_FIXMUL(-_hom->fwd[1][0], a22, r1, s1);
    _hom->inv[1][1] = QR_FIXMUL(_hom->fwd[0][0], a22, r1, s1);
    _hom->inv[2][0] = QR_FIXMUL(_hom->fwd[1][0], _hom->fwd[2][1],
        -QR_EXTMUL(_hom->fwd[1][1], _hom->fwd[2][0], r2), s2);
    _hom->inv[2][1] = QR_FIXMUL(_hom->fwd[0][1], _hom->fwd[2][0],
        -QR_EXTMUL(_hom->fwd[0][0], _hom->fwd[2][1], r2), s2);
    _hom->inv22 = QR_FIXMUL(_hom->fwd[0][0], _hom->fwd[1][1],
        -QR_EXTMUL(_hom->fwd[0][1], _hom->fwd[1][0], r2), s2);
    _hom->res = _res;
}


static int qr_hom_fit(qr_hom* _hom, qr_finder* _ul, qr_finder* _ur,
    qr_finder* _dl, qr_point _p[4], const qr_aff* _aff, isaac_ctx* _isaac,
    const unsigned char* _img, int _width, int _height) {
    qr_point* b;
    int       nb;
    int       cb;
    qr_point* r;
    int       nr;
    int       cr;
    qr_line   l[4];
    qr_point  q;
    qr_point  p;
    int       ox;
    int       oy;
    int       ru;
    int       rv;
    int       dru;
    int       drv;
    int       bu;
    int       bv;
    int       dbu;
    int       dbv;
    int       rx;
    int       ry;
    int       drxi;
    int       dryi;
    int       drxj;
    int       dryj;
    int       rdone;
    int       nrempty;
    int       rlastfit;
    int       bx;
    int       by;
    int       dbxi;
    int       dbyi;
    int       dbxj;
    int       dbyj;
    int       bdone;
    int       nbempty;
    int       blastfit;
    int       shift;
    int       round;
    int       version4;
    int       brx;
    int       bry;
    int       i;
    qr_finder_ransac(_ul, _aff, _isaac, 0);
    qr_finder_ransac(_dl, _aff, _isaac, 0);
    qr_line_fit_finder_pair(l[0], _aff, _ul, _dl, 0);
    if (qr_line_eval(l[0], _dl->c->pos[0], _dl->c->pos[1]) < 0 ||
        qr_line_eval(l[0], _ur->c->pos[0], _ur->c->pos[1]) < 0) {
        return -1;
    }
    qr_finder_ransac(_ul, _aff, _isaac, 2);
    qr_finder_ransac(_ur, _aff, _isaac, 2);
    qr_line_fit_finder_pair(l[2], _aff, _ul, _ur, 2);
    if (qr_line_eval(l[2], _dl->c->pos[0], _dl->c->pos[1]) < 0 ||
        qr_line_eval(l[2], _ur->c->pos[0], _ur->c->pos[1]) < 0) {
        return -1;
    }
    drv = _ur->size[1] >> 1;
    qr_finder_ransac(_ur, _aff, _isaac, 1);
    if (qr_line_fit_finder_edge(l[1], _ur, 1, _aff->res) >= 0) {
        if (qr_line_eval(l[1], _ul->c->pos[0], _ul->c->pos[1]) < 0 ||
            qr_line_eval(l[1], _dl->c->pos[0], _dl->c->pos[1]) < 0) {
            return -1;
        }
        if (qr_aff_line_step(_aff, l[1], 1, drv, &dru) < 0)return -1;
    }
    else dru = 0;
    ru = _ur->o[0] + 3 * _ur->size[0] - 2 * dru;
    rv = _ur->o[1] - 2 * drv;
    dbu = _dl->size[0] >> 1;
    qr_finder_ransac(_dl, _aff, _isaac, 3);
    if (qr_line_fit_finder_edge(l[3], _dl, 3, _aff->res) >= 0) {
        if (qr_line_eval(l[3], _ul->c->pos[0], _ul->c->pos[1]) < 0 ||
            qr_line_eval(l[3], _ur->c->pos[0], _ur->c->pos[1]) < 0) {
            return -1;
        }
        if (qr_aff_line_step(_aff, l[3], 0, dbu, &dbv) < 0)return -1;
    }
    else dbv = 0;
    bu = _dl->o[0] - 2 * dbu;
    bv = _dl->o[1] + 3 * _dl->size[1] - 2 * dbv;
    nr = rlastfit = _ur->ninliers[1];
    cr = nr + (_dl->o[1] - rv + drv - 1) / drv;
    r = (qr_point*)malloc(cr * sizeof(*r));
    for (i = 0; i < _ur->ninliers[1]; i++) {
        memcpy(r[i], _ur->edge_pts[1][i].pos, sizeof(r[i]));
    }
    nb = blastfit = _dl->ninliers[3];
    cb = nb + (_ur->o[0] - bu + dbu - 1) / dbu;
    b = (qr_point*)malloc(cb * sizeof(*b));
    for (i = 0; i < _dl->ninliers[3]; i++) {
        memcpy(b[i], _dl->edge_pts[3][i].pos, sizeof(b[i]));
    }
    ox = (_aff->x0 << _aff->res) + (1 << _aff->res - 1);
    oy = (_aff->y0 << _aff->res) + (1 << _aff->res - 1);
    rx = _aff->fwd[0][0] * ru + _aff->fwd[0][1] * rv + ox;
    ry = _aff->fwd[1][0] * ru + _aff->fwd[1][1] * rv + oy;
    drxi = _aff->fwd[0][0] * dru + _aff->fwd[0][1] * drv;
    dryi = _aff->fwd[1][0] * dru + _aff->fwd[1][1] * drv;
    drxj = _aff->fwd[0][0] * _ur->size[0];
    dryj = _aff->fwd[1][0] * _ur->size[0];
    bx = _aff->fwd[0][0] * bu + _aff->fwd[0][1] * bv + ox;
    by = _aff->fwd[1][0] * bu + _aff->fwd[1][1] * bv + oy;
    dbxi = _aff->fwd[0][0] * dbu + _aff->fwd[0][1] * dbv;
    dbyi = _aff->fwd[1][0] * dbu + _aff->fwd[1][1] * dbv;
    dbxj = _aff->fwd[0][1] * _dl->size[1];
    dbyj = _aff->fwd[1][1] * _dl->size[1];
    nrempty = nbempty = 0;
    for (;;) {
        int ret;
        int x0;
        int y0;
        int x1;
        int y1;
        rdone = rv >= QR_MINI(bv, _dl->o[1] + bv >> 1) || nrempty > 14;
        bdone = bu >= QR_MINI(ru, _ur->o[0] + ru >> 1) || nbempty > 14;
        if (!rdone && (bdone || rv < bu)) {
            x0 = rx + drxj >> _aff->res + QR_FINDER_SUBPREC;
            y0 = ry + dryj >> _aff->res + QR_FINDER_SUBPREC;
            x1 = rx - drxj >> _aff->res + QR_FINDER_SUBPREC;
            y1 = ry - dryj >> _aff->res + QR_FINDER_SUBPREC;
            if (nr >= cr) {
                cr = cr << 1 | 1;
                r = (qr_point*)realloc(r, cr * sizeof(*r));
            }
            ret = qr_finder_quick_crossing_check(_img, _width, _height, x0, y0, x1, y1, 1);
            if (!ret) {
                ret = qr_finder_locate_crossing(_img, _width, _height, x0, y0, x1, y1, 1, r[nr]);
            }
            if (ret >= 0) {
                if (!ret) {
                    qr_aff_unproject(q, _aff, r[nr][0], r[nr][1]);
                    ru = ru + q[0] >> 1;
                    if (q[1] + drv > rv)rv = rv + q[1] >> 1;
                    rx = _aff->fwd[0][0] * ru + _aff->fwd[0][1] * rv + ox;
                    ry = _aff->fwd[1][0] * ru + _aff->fwd[1][1] * rv + oy;
                    nr++;
                    if (nr > QR_MAXI(1, rlastfit + (rlastfit >> 2))) {
                        qr_line_fit_points(l[1], r, nr, _aff->res);
                        if (qr_aff_line_step(_aff, l[1], 1, drv, &dru) >= 0) {
                            drxi = _aff->fwd[0][0] * dru + _aff->fwd[0][1] * drv;
                            dryi = _aff->fwd[1][0] * dru + _aff->fwd[1][1] * drv;
                        }
                        rlastfit = nr;
                    }
                }
                nrempty = 0;
            }
            else nrempty++;
            ru += dru;
            if (rv + drv > rv)rv += drv;
            else nrempty = INT_MAX;
            rx += drxi;
            ry += dryi;
        }
        else if (!bdone) {
            x0 = bx + dbxj >> _aff->res + QR_FINDER_SUBPREC;
            y0 = by + dbyj >> _aff->res + QR_FINDER_SUBPREC;
            x1 = bx - dbxj >> _aff->res + QR_FINDER_SUBPREC;
            y1 = by - dbyj >> _aff->res + QR_FINDER_SUBPREC;
            if (nb >= cb) {
                cb = cb << 1 | 1;
                b = (qr_point*)realloc(b, cb * sizeof(*b));
            }
            ret = qr_finder_quick_crossing_check(_img, _width, _height, x0, y0, x1, y1, 1);
            if (!ret) {
                ret = qr_finder_locate_crossing(_img, _width, _height, x0, y0, x1, y1, 1, b[nb]);
            }
            if (ret >= 0) {
                if (!ret) {
                    qr_aff_unproject(q, _aff, b[nb][0], b[nb][1]);
                    /*Move the current point halfway towards the crossing.
                      We don't move the whole way to give us some robustness to noise.*/
                      /*But ensure that bu monotonically increases.*/
                    if (q[0] + dbu > bu)bu = bu + q[0] >> 1;
                    bv = bv + q[1] >> 1;
                    bx = _aff->fwd[0][0] * bu + _aff->fwd[0][1] * bv + ox;
                    by = _aff->fwd[1][0] * bu + _aff->fwd[1][1] * bv + oy;
                    nb++;
                    /*Re-fit the line to update the step direction periodically.*/
                    if (nb > QR_MAXI(1, blastfit + (blastfit >> 2))) {
                        qr_line_fit_points(l[3], b, nb, _aff->res);
                        if (qr_aff_line_step(_aff, l[3], 0, dbu, &dbv) >= 0) {
                            dbxi = _aff->fwd[0][0] * dbu + _aff->fwd[0][1] * dbv;
                            dbyi = _aff->fwd[1][0] * dbu + _aff->fwd[1][1] * dbv;
                        }
                        blastfit = nb;
                    }
                }
                nbempty = 0;
            }
            else nbempty++;
            if (bu + dbu > bu)bu += dbu;
            else nbempty = INT_MAX;
            bv += dbv;
            bx += dbxi;
            by += dbyi;
        }
        else break;
    }
    if (nr > 1)qr_line_fit_points(l[1], r, nr, _aff->res);
    else {
        qr_aff_project(p, _aff, _ur->o[0] + 3 * _ur->size[0], _ur->o[1]);
        shift = QR_MAXI(0,
            qr_ilog(QR_MAXI(abs(_aff->fwd[0][1]), abs(_aff->fwd[1][1])))
            - (_aff->res + 1 >> 1));
        round = (1 << shift) >> 1;
        l[1][0] = _aff->fwd[1][1] + round >> shift;
        l[1][1] = -_aff->fwd[0][1] + round >> shift;
        l[1][2] = -(l[1][0] * p[0] + l[1][1] * p[1]);
    }
    free(r);
    if (nb > 1)qr_line_fit_points(l[3], b, nb, _aff->res);
    else {
        qr_aff_project(p, _aff, _dl->o[0], _dl->o[1] + 3 * _dl->size[1]);
        shift = QR_MAXI(0,
            qr_ilog(QR_MAXI(abs(_aff->fwd[0][1]), abs(_aff->fwd[1][1])))
            - (_aff->res + 1 >> 1));
        round = (1 << shift) >> 1;
        l[3][0] = _aff->fwd[1][0] + round >> shift;
        l[3][1] = -_aff->fwd[0][0] + round >> shift;
        l[3][2] = -(l[1][0] * p[0] + l[1][1] * p[1]);
    }
    free(b);
    for (i = 0; i < 4; i++) {
        if (qr_line_isect(_p[i], l[i & 1], l[2 + (i >> 1)]) < 0)return -1;
        if (_p[i][0] < -_width << QR_FINDER_SUBPREC ||
            _p[i][0] >= _width << QR_FINDER_SUBPREC + 1 ||
            _p[i][1] < -_height << QR_FINDER_SUBPREC ||
            _p[i][1] >= _height << QR_FINDER_SUBPREC + 1) {
            return -1;
        }
    }
    brx = _p[3][0];
    bry = _p[3][1];
    version4 = _ul->eversion[0] + _ul->eversion[1] + _ur->eversion[0] + _dl->eversion[1];
    if (version4 > 4) {
        qr_hom_cell cell;
        qr_point    p3;
        int         dim;
        dim = 17 + version4;
        qr_hom_cell_init(&cell, 0, 0, dim - 1, 0, 0, dim - 1, dim - 1, dim - 1,
            _p[0][0], _p[0][1], _p[1][0], _p[1][1],
            _p[2][0], _p[2][1], _p[3][0], _p[3][1]);
        if (qr_alignment_pattern_search(p3, &cell, dim - 7, dim - 7, 4,
            _img, _width, _height) >= 0) {
            long long w;
            long long mask;
            int       c21;
            int       dx21;
            int       dy21;
            c21 = _p[2][0] * _p[1][1] - _p[2][1] * _p[1][0];
            dx21 = _p[2][0] - _p[1][0];
            dy21 = _p[2][1] - _p[1][1];
            w = QR_EXTMUL(dim - 7, c21,
                QR_EXTMUL(dim - 13, _p[0][0] * dy21 - _p[0][1] * dx21,
                    QR_EXTMUL(6, p3[0] * dy21 - p3[1] * dx21, 0)));
            if (w == 0)return -1;
            mask = QR_SIGNMASK(w);
            w = w + mask ^ mask;
            brx = (int)QR_DIVROUND(QR_EXTMUL((dim - 7) * _p[0][0], p3[0] * dy21,
                QR_EXTMUL((dim - 13) * p3[0], c21 - _p[0][1] * dx21,
                    QR_EXTMUL(6 * _p[0][0], c21 - p3[1] * dx21, 0))) + mask ^ mask, w);
            bry = (int)QR_DIVROUND(QR_EXTMUL((dim - 7) * _p[0][1], -p3[1] * dx21,
                QR_EXTMUL((dim - 13) * p3[1], c21 + _p[0][0] * dy21,
                    QR_EXTMUL(6 * _p[0][1], c21 + p3[0] * dy21, 0))) + mask ^ mask, w);
        }
    }

    qr_hom_init(_hom, _p[0][0], _p[0][1], _p[1][0], _p[1][1],
        _p[2][0], _p[2][1], brx, bry, QR_HOM_BITS);
    return 0;
}


static int qr_hom_unproject(qr_point _q, const qr_hom* _hom, int _x, int _y) {
    int x;
    int y;
    int w;
    _x -= _hom->x0;
    _y -= _hom->y0;
    x = _hom->inv[0][0] * _x + _hom->inv[0][1] * _y;
    y = _hom->inv[1][0] * _x + _hom->inv[1][1] * _y;
    w = _hom->inv[2][0] * _x + _hom->inv[2][1] * _y
        + _hom->inv22 + (1 << _hom->res - 1) >> _hom->res;
    if (w == 0) {
        _q[0] = x < 0 ? INT_MIN : INT_MAX;
        _q[1] = y < 0 ? INT_MIN : INT_MAX;
        return -1;
    }
    else {
        if (w < 0) {
            x = -x;
            y = -y;
            w = -w;
        }
        _q[0] = QR_DIVROUND(x, w);
        _q[1] = QR_DIVROUND(y, w);
    }
    return 0;
}

static void qr_finder_edge_pts_hom_classify(qr_finder* _f, const qr_hom* _hom) {
    qr_finder_center* c;
    int               i;
    int               e;
    c = _f->c;
    for (e = 0; e < 4; e++)_f->nedge_pts[e] = 0;
    for (i = 0; i < c->nedge_pts; i++) {
        qr_point q;
        int      d;
        if (qr_hom_unproject(q, _hom,
            c->edge_pts[i].pos[0], c->edge_pts[i].pos[1]) >= 0) {
            qr_point_translate(q, -_f->o[0], -_f->o[1]);
            d = abs(q[1]) > abs(q[0]);
            e = d << 1 | (q[d] >= 0);
            _f->nedge_pts[e]++;
            c->edge_pts[i].edge = e;
            c->edge_pts[i].extent = q[d];
        }
        else {
            c->edge_pts[i].edge = 4;
            c->edge_pts[i].extent = q[0];
        }
    }
    qsort(c->edge_pts, c->nedge_pts, sizeof(*c->edge_pts), qr_cmp_edge_pt);
    _f->edge_pts[0] = c->edge_pts;
    for (e = 1; e < 4; e++)_f->edge_pts[e] = _f->edge_pts[e - 1] + _f->nedge_pts[e - 1];
}

static void qr_hom_fproject(qr_point _p, const qr_hom* _hom,
    int _x, int _y, int _w) {
    if (_w == 0) {
        _p[0] = _x < 0 ? INT_MIN : INT_MAX;
        _p[1] = _y < 0 ? INT_MIN : INT_MAX;
    }
    else {
        if (_w < 0) {
            _x = -_x;
            _y = -_y;
            _w = -_w;
        }
        _p[0] = QR_DIVROUND(_x, _w) + _hom->x0;
        _p[1] = QR_DIVROUND(_y, _w) + _hom->y0;
    }
}

static const unsigned BCH18_6_CODES[34] = {
                                                          0x07C94,
  0x085BC,0x09A99,0x0A4D3,0x0BBF6,0x0C762,0x0D847,0x0E60D,0x0F928,
  0x10B78,0x1145D,0x12A17,0x13532,0x149A6,0x15683,0x168C9,0x177EC,
  0x18EC4,0x191E1,0x1AFAB,0x1B08E,0x1CC1A,0x1D33F,0x1ED75,0x1F250,
  0x209D5,0x216F0,0x228BA,0x2379F,0x24B0B,0x2542E,0x26A64,0x27541,
  0x28C69
};

static int bch18_6_correct(unsigned* _y) {
    unsigned x;
    unsigned y;
    int      nerrs;
    y = *_y;
    x = y >> 12;
    if (x >= 7 && x <= 40) {
        nerrs = qr_hamming_dist(y, BCH18_6_CODES[x - 7], 4);
        if (nerrs < 4) {
            *_y = BCH18_6_CODES[x - 7];
            return nerrs;
        }
    }
    for (x = 0; x < 34; x++)if (x + 7 != y >> 12) {
        nerrs = qr_hamming_dist(y, BCH18_6_CODES[x], 4);
        if (nerrs < 4) {
            *_y = BCH18_6_CODES[x];
            return nerrs;
        }
    }
    return -1;
}

static int qr_finder_version_decode(qr_finder* _f, const qr_hom* _hom,
    const unsigned char* _img, int _width, int _height, int _dir) {
    qr_point q;
    unsigned v;
    int      x0;
    int      y0;
    int      w0;
    int      dxi;
    int      dyi;
    int      dwi;
    int      dxj;
    int      dyj;
    int      dwj;
    int      ret;
    int      i;
    int      j;
    int      k;
    v = 0;
    q[_dir] = _f->o[_dir] - 7 * _f->size[_dir];
    q[1 - _dir] = _f->o[1 - _dir] - 3 * _f->size[1 - _dir];
    x0 = _hom->fwd[0][0] * q[0] + _hom->fwd[0][1] * q[1];
    y0 = _hom->fwd[1][0] * q[0] + _hom->fwd[1][1] * q[1];
    w0 = _hom->fwd[2][0] * q[0] + _hom->fwd[2][1] * q[1] + _hom->fwd22;
    dxi = _hom->fwd[0][1 - _dir] * _f->size[1 - _dir];
    dyi = _hom->fwd[1][1 - _dir] * _f->size[1 - _dir];
    dwi = _hom->fwd[2][1 - _dir] * _f->size[1 - _dir];
    dxj = _hom->fwd[0][_dir] * _f->size[_dir];
    dyj = _hom->fwd[1][_dir] * _f->size[_dir];
    dwj = _hom->fwd[2][_dir] * _f->size[_dir];
    for (k = i = 0; i < 6; i++) {
        int x;
        int y;
        int w;
        x = x0;
        y = y0;
        w = w0;
        for (j = 0; j < 3; j++, k++) {
            qr_point p;
            qr_hom_fproject(p, _hom, x, y, w);
            v |= qr_img_get_bit(_img, _width, _height, p[0], p[1]) << k;
            x += dxj;
            y += dyj;
            w += dwj;
        }
        x0 += dxi;
        y0 += dyi;
        w0 += dwi;
    }
    ret = bch18_6_correct(&v);
    return ret >= 0 ? (int)(v >> 12) : ret;
}

static const unsigned char gf16_exp[31] = {
  1,2,4,8,3,6,12,11,5,10,7,14,15,13,9,1,2,4,8,3,6,12,11,5,10,7,14,15,13,9,1
};

static int bch15_5_calc_syndrome(unsigned _s[3], unsigned _y) {
    unsigned p;
    int      i;
    int      j;
    p = 0;
    for (i = 0; i < 15; i++)if (_y & 1 << i)p ^= gf16_exp[i];
    _s[0] = p;
    p = 0;
    for (i = 0; i < 3; i++)for (j = 0; j < 5; j++)if (_y & 1 << 5 * i + j)p ^= gf16_exp[j * 3];
    _s[1] = p;
    p = 0;
    for (i = 0; i < 5; i++)for (j = 0; j < 3; j++)if (_y & 1 << 3 * i + j)p ^= gf16_exp[j * 5];
    _s[2] = p;
    return _s[0] != 0 || _s[1] != 0 || _s[2] != 0;
}

static const signed char gf16_log[16] = {
  -1,0,1,4,2,8,5,10,3,14,9,7,6,13,11,12
};

static unsigned gf16_mul(unsigned _a, unsigned _b) {
    return _a == 0 || _b == 0 ? 0 : gf16_exp[gf16_log[_a] + gf16_log[_b]];
}

static unsigned gf16_div(unsigned _a, unsigned _b) {
    return _a == 0 ? 0 : gf16_exp[gf16_log[_a] + 15 - gf16_log[_b]];
}

static int bch15_5_calc_omega(unsigned _o[3], unsigned _s[3]) {
    unsigned s02;
    unsigned tt;
    unsigned dd;
    int      d;
    _o[0] = _s[0];
    s02 = gf16_mul(_s[0], _s[0]);
    dd = _s[1] ^ gf16_mul(_s[0], s02);
    tt = _s[2] ^ gf16_mul(s02, _s[1]);
    _o[1] = dd ? gf16_div(tt, dd) : 0;
    _o[2] = dd ^ gf16_mul(_s[0], _o[1]);
    for (d = 3; d > 0 && !_o[d - 1]; d--);
    return d;
}

static unsigned gf16_hmul(unsigned _a, unsigned _logb) {
    return _a == 0 ? 0 : gf16_exp[gf16_log[_a] + _logb];
}

static int bch15_5_calc_epos(unsigned _epos[3], unsigned _s[3]) {
    unsigned o[3];
    int      nerrors;
    int      d;
    int      i;
    d = bch15_5_calc_omega(o, _s);
    nerrors = 0;
    if (d == 1)_epos[nerrors++] = gf16_log[o[0]];
    else if (d > 0) {
        for (i = 0; i < 15; i++) {
            int i2;
            i2 = gf16_log[gf16_exp[i << 1]];
            if (!(gf16_exp[i + i2] ^ gf16_hmul(o[0], i2) ^ gf16_hmul(o[1], i) ^ o[2])) {
                _epos[nerrors++] = i;
            }
        }
        if (nerrors < d)return -1;
    }
    return nerrors;
}

unsigned bch15_5_encode(unsigned _x) {
    return (-(_x & 1) & 0x0537) ^ (-(_x >> 1 & 1) & 0x0A6E) ^ (-(_x >> 2 & 1) & 0x11EB) ^
        (-(_x >> 3 & 1) & 0x23D6) ^ (-(_x >> 4 & 1) & 0x429B);
}

int bch15_5_correct(unsigned* _y) {
    unsigned s[3];
    unsigned epos[3];
    unsigned y;
    int      nerrors;
    int      i;
    y = *_y;
    if (!bch15_5_calc_syndrome(s, y))return 0;
    nerrors = bch15_5_calc_epos(epos, s);
    if (nerrors > 0) {
        for (i = 0; i < nerrors; i++)y ^= 1 << epos[i];
        if (bch15_5_encode(y >> 10) == y) {
            *_y = y;
            return nerrors;
        }
    }
    return -1;
}

static int qr_finder_fmt_info_decode(qr_finder* _ul, qr_finder* _ur,
    qr_finder* _dl, const qr_hom* _hom,
    const unsigned char* _img, int _width, int _height) {
    qr_point p;
    unsigned lo[2];
    unsigned hi[2];
    int      u;
    int      v;
    int      x;
    int      y;
    int      w;
    int      dx;
    int      dy;
    int      dw;
    int      fmt_info[4];
    int      count[4];
    int      nerrs[4];
    int      nfmt_info;
    int      besti;
    int      imax;
    int      di;
    int      i;
    int      k;
    lo[0] = 0;
    u = _ul->o[0] + 5 * _ul->size[0];
    v = _ul->o[1] - 3 * _ul->size[1];
    x = _hom->fwd[0][0] * u + _hom->fwd[0][1] * v;
    y = _hom->fwd[1][0] * u + _hom->fwd[1][1] * v;
    w = _hom->fwd[2][0] * u + _hom->fwd[2][1] * v + _hom->fwd22;
    dx = _hom->fwd[0][1] * _ul->size[1];
    dy = _hom->fwd[1][1] * _ul->size[1];
    dw = _hom->fwd[2][1] * _ul->size[1];
    for (k = i = 0;; i++) {
        if (i != 6) {
            qr_hom_fproject(p, _hom, x, y, w);
            lo[0] |= qr_img_get_bit(_img, _width, _height, p[0], p[1]) << k++;
            if (i >= 8)break;
        }
        x += dx;
        y += dy;
        w += dw;
    }
    hi[0] = 0;
    dx = -_hom->fwd[0][0] * _ul->size[0];
    dy = -_hom->fwd[1][0] * _ul->size[0];
    dw = -_hom->fwd[2][0] * _ul->size[0];
    while (i-- > 0) {
        x += dx;
        y += dy;
        w += dw;
        if (i != 6) {
            qr_hom_fproject(p, _hom, x, y, w);
            hi[0] |= qr_img_get_bit(_img, _width, _height, p[0], p[1]) << k++;
        }
    }
    lo[1] = 0;
    u = _ur->o[0] + 3 * _ur->size[0];
    v = _ur->o[1] + 5 * _ur->size[1];
    x = _hom->fwd[0][0] * u + _hom->fwd[0][1] * v;
    y = _hom->fwd[1][0] * u + _hom->fwd[1][1] * v;
    w = _hom->fwd[2][0] * u + _hom->fwd[2][1] * v + _hom->fwd22;
    dx = -_hom->fwd[0][0] * _ur->size[0];
    dy = -_hom->fwd[1][0] * _ur->size[0];
    dw = -_hom->fwd[2][0] * _ur->size[0];
    for (k = 0; k < 8; k++) {
        qr_hom_fproject(p, _hom, x, y, w);
        lo[1] |= qr_img_get_bit(_img, _width, _height, p[0], p[1]) << k;
        x += dx;
        y += dy;
        w += dw;
    }
    hi[1] = 0;
    u = _dl->o[0] + 5 * _dl->size[0];
    v = _dl->o[1] - 3 * _dl->size[1];
    x = _hom->fwd[0][0] * u + _hom->fwd[0][1] * v;
    y = _hom->fwd[1][0] * u + _hom->fwd[1][1] * v;
    w = _hom->fwd[2][0] * u + _hom->fwd[2][1] * v + _hom->fwd22;
    dx = _hom->fwd[0][1] * _dl->size[1];
    dy = _hom->fwd[1][1] * _dl->size[1];
    dw = _hom->fwd[2][1] * _dl->size[1];
    for (k = 8; k < 15; k++) {
        qr_hom_fproject(p, _hom, x, y, w);
        hi[1] |= qr_img_get_bit(_img, _width, _height, p[0], p[1]) << k;
        x += dx;
        y += dy;
        w += dw;
    }
    imax = 2 << (hi[0] != hi[1]);
    di = 1 + (lo[0] == lo[1]);
    nfmt_info = 0;
    for (i = 0; i < imax; i += di) {
        unsigned v;
        int      ret;
        int      j;
        v = (lo[i & 1] | hi[i >> 1]) ^ 0x5412;
        ret = bch15_5_correct(&v);
        v >>= 10;
        if (ret < 0)ret = 4;
        for (j = 0;; j++) {
            if (j >= nfmt_info) {
                fmt_info[j] = v;
                count[j] = 1;
                nerrs[j] = ret;
                nfmt_info++;
                break;
            }
            if (fmt_info[j] == (int)v) {
                count[j]++;
                if (ret < nerrs[j])nerrs[j] = ret;
                break;
            }
        }
    }
    besti = 0;
    for (i = 1; i < nfmt_info; i++) {
        if (nerrs[besti] > 3 && nerrs[i] <= 3 ||
            count[i] > count[besti] || count[i] == count[besti] && nerrs[i] < nerrs[besti]) {
            besti = i;
        }
    }
    return nerrs[besti] < 4 ? fmt_info[besti] : -1;
}

static void qr_sampling_grid_fp_mask_rect(qr_sampling_grid* _grid, int _dim,
    int _u, int _v, int _w, int _h) {
    int i;
    int j;
    int stride;
    stride = _dim + QR_INT_BITS - 1 >> QR_INT_LOGBITS;
    for (j = _u; j < _u + _w; j++)for (i = _v; i < _v + _h; i++) {
        _grid->fpmask[j * stride + (i >> QR_INT_LOGBITS)] |= 1 << (i & QR_INT_BITS - 1);
    }
}

static const unsigned char QR_ALIGNMENT_SPACING[34] = {
  16,18,20,22,24,26,28,
  20,22,24,24,26,28,28,
  22,24,24,26,26,28,28,
  24,24,26,26,26,28,28,
  24,26,26,26,28,28
};

static void qr_hom_cell_project(qr_point _p, const qr_hom_cell* _cell,
    int _u, int _v, int _res) {
    _u -= _cell->u0 << _res;
    _v -= _cell->v0 << _res;
    qr_hom_cell_fproject(_p, _cell,
        _cell->fwd[0][0] * _u + _cell->fwd[0][1] * _v + (_cell->fwd[0][2] << _res),
        _cell->fwd[1][0] * _u + _cell->fwd[1][1] * _v + (_cell->fwd[1][2] << _res),
        _cell->fwd[2][0] * _u + _cell->fwd[2][1] * _v + (_cell->fwd[2][2] << _res));
}

static void qr_sampling_grid_init(qr_sampling_grid* _grid, int _version,
    const qr_point _ul_pos, const qr_point _ur_pos, const qr_point _dl_pos,
    qr_point _p[4], const unsigned char* _img, int _width, int _height) {
    qr_hom_cell          base_cell;
    int                  align_pos[7];
    int                  dim;
    int                  nalign;
    int                  i;
    dim = 17 + (_version << 2);
    nalign = (_version / 7) + 2;
    qr_hom_cell_init(&base_cell, 0, 0, dim - 1, 0, 0, dim - 1, dim - 1, dim - 1,
        _p[0][0], _p[0][1], _p[1][0], _p[1][1], _p[2][0], _p[2][1], _p[3][0], _p[3][1]);
    _grid->ncells = nalign - 1;
    _grid->cells[0] = (qr_hom_cell*)malloc(
        (nalign - 1) * (nalign - 1) * sizeof(*_grid->cells[0]));
    for (i = 1; i < _grid->ncells; i++)_grid->cells[i] = _grid->cells[i - 1] + _grid->ncells;
    _grid->fpmask = (unsigned*)calloc(dim,
        (dim + QR_INT_BITS - 1 >> QR_INT_LOGBITS) * sizeof(*_grid->fpmask));
    qr_sampling_grid_fp_mask_rect(_grid, dim, 0, 0, 9, 9);
    qr_sampling_grid_fp_mask_rect(_grid, dim, 0, dim - 8, 9, 8);
    qr_sampling_grid_fp_mask_rect(_grid, dim, dim - 8, 0, 8, 9);
    if (_version > 6) {
        qr_sampling_grid_fp_mask_rect(_grid, dim, 0, dim - 11, 6, 3);
        qr_sampling_grid_fp_mask_rect(_grid, dim, dim - 11, 0, 3, 6);
    }
    qr_sampling_grid_fp_mask_rect(_grid, dim, 9, 6, dim - 17, 1);
    qr_sampling_grid_fp_mask_rect(_grid, dim, 6, 9, 1, dim - 17);
    if (_version < 2)memcpy(_grid->cells[0], &base_cell, sizeof(base_cell));
    else {
        qr_point* q;
        qr_point* p;
        int       j;
        int       k;
        q = (qr_point*)malloc(nalign * nalign * sizeof(*q));
        p = (qr_point*)malloc(nalign * nalign * sizeof(*p));
        align_pos[0] = 6;
        align_pos[nalign - 1] = dim - 7;
        if (_version > 6) {
            int d;
            d = QR_ALIGNMENT_SPACING[_version - 7];
            for (i = nalign - 1; i-- > 1;)align_pos[i] = align_pos[i + 1] - d;
        }
        q[0][0] = 3;
        q[0][1] = 3;
        p[0][0] = _ul_pos[0];
        p[0][1] = _ul_pos[1];
        q[nalign - 1][0] = dim - 4;
        q[nalign - 1][1] = 3;
        p[nalign - 1][0] = _ur_pos[0];
        p[nalign - 1][1] = _ur_pos[1];
        q[(nalign - 1) * nalign][0] = 3;
        q[(nalign - 1) * nalign][1] = dim - 4;
        p[(nalign - 1) * nalign][0] = _dl_pos[0];
        p[(nalign - 1) * nalign][1] = _dl_pos[1];
        for (k = 1; k < 2 * nalign - 1; k++) {
            int jmin;
            int jmax;
            jmax = QR_MINI(k, nalign - 1) - (k == nalign - 1);
            jmin = QR_MAXI(0, k - (nalign - 1)) + (k == nalign - 1);
            for (j = jmin; j <= jmax; j++) {
                qr_hom_cell* cell;
                int          u;
                int          v;
                int          k;
                i = jmax - (j - jmin);
                k = i * nalign + j;
                u = align_pos[j];
                v = align_pos[i];
                q[k][0] = u;
                q[k][1] = v;
                qr_sampling_grid_fp_mask_rect(_grid, dim, u - 2, v - 2, 5, 5);
                if (i > 1 && j > 1) {
                    qr_point p0;
                    qr_point p1;
                    qr_point p2;
                    qr_hom_cell_project(p0, _grid->cells[i - 2] + j - 1, u, v, 0);
                    qr_hom_cell_project(p1, _grid->cells[i - 2] + j - 2, u, v, 0);
                    qr_hom_cell_project(p2, _grid->cells[i - 1] + j - 2, u, v, 0);
                    QR_SORT2I(p0[0], p1[0]);
                    QR_SORT2I(p0[1], p1[1]);
                    QR_SORT2I(p1[0], p2[0]);
                    QR_SORT2I(p1[1], p2[1]);
                    QR_SORT2I(p0[0], p1[0]);
                    QR_SORT2I(p0[1], p1[1]);
                    cell = _grid->cells[i - 1] + j - 1;
                    qr_hom_cell_init(cell,
                        q[k - nalign - 1][0], q[k - nalign - 1][1], q[k - nalign][0], q[k - nalign][1],
                        q[k - 1][0], q[k - 1][1], q[k][0], q[k][1],
                        p[k - nalign - 1][0], p[k - nalign - 1][1], p[k - nalign][0], p[k - nalign][1],
                        p[k - 1][0], p[k - 1][1], p1[0], p1[1]);
                }
                else if (i > 1 && j > 0)cell = _grid->cells[i - 2] + j - 1;
                else if (i > 0 && j > 1)cell = _grid->cells[i - 1] + j - 2;
                else cell = &base_cell;
                qr_alignment_pattern_search(p[k], cell, u, v, 2, _img, _width, _height);
                if (i > 0 && j > 0) {
                    qr_hom_cell_init(_grid->cells[i - 1] + j - 1,
                        q[k - nalign - 1][0], q[k - nalign - 1][1], q[k - nalign][0], q[k - nalign][1],
                        q[k - 1][0], q[k - 1][1], q[k][0], q[k][1],
                        p[k - nalign - 1][0], p[k - nalign - 1][1], p[k - nalign][0], p[k - nalign][1],
                        p[k - 1][0], p[k - 1][1], p[k][0], p[k][1]);
                }
            }
        }
        free(q);
        free(p);
    }
    memcpy(_grid->cell_limits, align_pos + 1,
        (_grid->ncells - 1) * sizeof(*_grid->cell_limits));
    _grid->cell_limits[_grid->ncells - 1] = dim;
    qr_hom_cell_project(_p[0], _grid->cells[0] + 0, -1, -1, 1);
    qr_hom_cell_project(_p[1], _grid->cells[0] + _grid->ncells - 1, (dim << 1) - 1, -1, 1);
    qr_hom_cell_project(_p[2], _grid->cells[_grid->ncells - 1] + 0, -1, (dim << 1) - 1, 1);
    qr_hom_cell_project(_p[3], _grid->cells[_grid->ncells - 1] + _grid->ncells - 1,
        (dim << 1) - 1, (dim << 1) - 1, 1);
    for (i = 0; i < 4; i++) {
        _p[i][0] = QR_CLAMPI(-_width << QR_FINDER_SUBPREC, _p[i][0],
            _width << QR_FINDER_SUBPREC + 1);
        _p[i][1] = QR_CLAMPI(-_height << QR_FINDER_SUBPREC, _p[i][1],
            _height << QR_FINDER_SUBPREC + 1);
    }
}

static void qr_data_mask_fill(unsigned* _mask, int _dim, int _pattern) {
    int stride;
    int i;
    int j;
    stride = _dim + QR_INT_BITS - 1 >> QR_INT_LOGBITS;
    switch (_pattern) {
    case 0: {
        int m;
        m = 0x55;
        for (j = 0; j < _dim; j++) {
            memset(_mask + j * stride, m, stride * sizeof(*_mask));
            m ^= 0xFF;
        }
    }break;
    case 1:memset(_mask, 0x55, _dim * stride * sizeof(*_mask)); break;
    case 2: {
        unsigned m;
        m = 0xFF;
        for (j = 0; j < _dim; j++) {
            memset(_mask + j * stride, m & 0xFF, stride * sizeof(*_mask));
            m = m << 8 | m >> 16;
        }
    }break;
    case 3: {
        unsigned mi;
        unsigned mj;
        mj = 0;
        for (i = 0; i < (QR_INT_BITS + 2) / 3; i++)mj |= 1 << 3 * i;
        for (j = 0; j < _dim; j++) {
            mi = mj;
            for (i = 0; i < stride; i++) {
                _mask[j * stride + i] = mi;
                mi = mi >> QR_INT_BITS % 3 | mi << 3 - QR_INT_BITS % 3;
            }
            mj = mj >> 1 | mj << 2;
        }
    }break;
    case 4: {
        unsigned m;
        m = 7;
        for (j = 0; j < _dim; j++) {
            memset(_mask + j * stride, (0xCC ^ -(m & 1)) & 0xFF, stride * sizeof(*_mask));
            m = m >> 1 | m << 5;
        }
    }break;
    case 5: {
        for (j = 0; j < _dim; j++) {
            unsigned m;
            m = 0;
            for (i = 0; i < 6; i++)m |= !((i * j) % 6) << i;
            for (i = 6; i < QR_INT_BITS; i <<= 1)m |= m << i;
            for (i = 0; i < stride; i++) {
                _mask[j * stride + i] = m;
                m = m >> QR_INT_BITS % 6 | m << 6 - QR_INT_BITS % 6;
            }
        }
    }break;
    case 6: {
        for (j = 0; j < _dim; j++) {
            unsigned m;
            m = 0;
            for (i = 0; i < 6; i++)m |= ((i * j) % 3 + i * j + 1 & 1) << i;
            for (i = 6; i < QR_INT_BITS; i <<= 1)m |= m << i;
            for (i = 0; i < stride; i++) {
                _mask[j * stride + i] = m;
                m = m >> QR_INT_BITS % 6 | m << 6 - QR_INT_BITS % 6;
            }
        }
    }break;
    default: {
        for (j = 0; j < _dim; j++) {
            unsigned m;
            m = 0;
            for (i = 0; i < 6; i++)m |= ((i * j) % 3 + i + j + 1 & 1) << i;
            for (i = 6; i < QR_INT_BITS; i <<= 1)m |= m << i;
            for (i = 0; i < stride; i++) {
                _mask[j * stride + i] = m;
                m = m >> QR_INT_BITS % 6 | m << 6 - QR_INT_BITS % 6;
            }
        }
    }break;
    }
}

static int qr_sampling_grid_is_in_fp(const qr_sampling_grid* _grid, int _dim,
    int _u, int _v) {
    return _grid->fpmask[_u * (_dim + QR_INT_BITS - 1 >> QR_INT_LOGBITS)
        + (_v >> QR_INT_LOGBITS)] >> (_v & QR_INT_BITS - 1) & 1;
}

static void qr_sampling_grid_sample(const qr_sampling_grid* _grid,
    unsigned* _data_bits, int _dim, int _fmt_info,
    const unsigned char* _img, int _width, int _height) {
    int stride;
    int u0;
    int u1;
    int j;
    qr_data_mask_fill(_data_bits, _dim, _fmt_info & 7);
    stride = _dim + QR_INT_BITS - 1 >> QR_INT_LOGBITS;
    u0 = 0;
    for (j = 0; j < _grid->ncells; j++) {
        int i;
        int v0;
        int v1;
        u1 = _grid->cell_limits[j];
        v0 = 0;
        for (i = 0; i < _grid->ncells; i++) {
            qr_hom_cell* cell;
            int          x0;
            int          y0;
            int          w0;
            int          u;
            int          du;
            int          dv;
            v1 = _grid->cell_limits[i];
            cell = _grid->cells[i] + j;
            du = u0 - cell->u0;
            dv = v0 - cell->v0;
            x0 = cell->fwd[0][0] * du + cell->fwd[0][1] * dv + cell->fwd[0][2];
            y0 = cell->fwd[1][0] * du + cell->fwd[1][1] * dv + cell->fwd[1][2];
            w0 = cell->fwd[2][0] * du + cell->fwd[2][1] * dv + cell->fwd[2][2];
            for (u = u0; u < u1; u++) {
                int x;
                int y;
                int w;
                int v;
                x = x0;
                y = y0;
                w = w0;
                for (v = v0; v < v1; v++) {
                    if (!qr_sampling_grid_is_in_fp(_grid, _dim, u, v)) {
                        qr_point p;
                        qr_hom_cell_fproject(p, cell, x, y, w);
                        _data_bits[u * stride + (v >> QR_INT_LOGBITS)] ^=
                            qr_img_get_bit(_img, _width, _height, p[0], p[1]) << (v & QR_INT_BITS - 1);
                    }
                    x += cell->fwd[0][1];
                    y += cell->fwd[1][1];
                    w += cell->fwd[2][1];
                }
                x0 += cell->fwd[0][0];
                y0 += cell->fwd[1][0];
                w0 += cell->fwd[2][0];
            }
            v0 = v1;
        }
        u0 = u1;
    }
}

static const unsigned char QR_RS_NBLOCKS[40][4] = {
  { 1, 1, 1, 1},{ 1, 1, 1, 1},{ 1, 1, 2, 2},{ 1, 2, 2, 4},
  { 1, 2, 4, 4},{ 2, 4, 4, 4},{ 2, 4, 6, 5},{ 2, 4, 6, 6},
  { 2, 5, 8, 8},{ 4, 5, 8, 8},{ 4, 5, 8,11},{ 4, 8,10,11},
  { 4, 9,12,16},{ 4, 9,16,16},{ 6,10,12,18},{ 6,10,17,16},
  { 6,11,16,19},{ 6,13,18,21},{ 7,14,21,25},{ 8,16,20,25},
  { 8,17,23,25},{ 9,17,23,34},{ 9,18,25,30},{10,20,27,32},
  {12,21,29,35},{12,23,34,37},{12,25,34,40},{13,26,35,42},
  {14,28,38,45},{15,29,40,48},{16,31,43,51},{17,33,45,54},
  {18,35,48,57},{19,37,51,60},{19,38,53,63},{20,40,56,66},
  {21,43,59,70},{22,45,62,74},{24,47,65,77},{25,49,68,81}
};

static const unsigned char QR_RS_NPAR_VALS[71] = {
    /*[ 0]*/ 7,10,13,17,
    /*[ 4]*/10,16,22, 28,26,26, 26,22, 24,22,22, 26,24,18,22,
    /*[19]*/15,26,18, 22,24, 30,24,20,24,
    /*[28]*/18,16,24, 28, 28, 28,28,30,24,
    /*[37]*/20,18, 18,26, 24,28,24, 30,26,28, 28, 26,28,30, 30,22,20,24,
    /*[55]*/20,18,26,16,
    /*[59]*/20,30,28, 24,22,26, 28,26, 30,28,30,30
};

static const unsigned char QR_RS_NPAR_OFFS[40] = {
   0, 4,19,55,15,28,37,12,51,39,
  59,62,10,24,22,41,31,44, 7,65,
  47,33,67,67,48,32,67,67,67,67,
  67,67,67,67,67,67,67,67,67,67
};

static int qr_code_ncodewords(unsigned _version) {
    unsigned nalign;
    if (_version == 1)return 26;
    nalign = (_version / 7) + 2;
    return (_version << 4) * (_version + 8)
        - (5 * nalign) * (5 * nalign - 2) + 36 * (_version < 7) + 83 >> 3;
}

static void qr_samples_unpack(unsigned char** _blocks, int _nblocks,
    int _nshort_data, int _nshort_blocks, const unsigned* _data_bits,
    const unsigned* _fp_mask, int _dim) {
    unsigned bits;
    int      biti;
    int      stride;
    int      blocki;
    int      blockj;
    int      i;
    int      j;
    stride = _dim + QR_INT_BITS - 1 >> QR_INT_LOGBITS;
    if (_nshort_blocks >= _nblocks)_nshort_blocks = 0;
    bits = 0;
    for (blocki = blockj = biti = 0, j = _dim - 1; j > 0; j -= 2) {
        unsigned data1;
        unsigned data2;
        unsigned fp_mask1;
        unsigned fp_mask2;
        int      nbits;
        int      l;
        nbits = (_dim - 1 & QR_INT_BITS - 1) + 1;
        l = j * stride;
        for (i = stride; i-- > 0;) {
            data1 = _data_bits[l + i];
            fp_mask1 = _fp_mask[l + i];
            data2 = _data_bits[l + i - stride];
            fp_mask2 = _fp_mask[l + i - stride];
            while (nbits-- > 0) {
                if (!(fp_mask1 >> nbits & 1)) {
                    bits = bits << 1 | data1 >> nbits & 1;
                    biti++;
                }
                if (!(fp_mask2 >> nbits & 1)) {
                    bits = bits << 1 | data2 >> nbits & 1;
                    biti++;
                }
                if (biti >= 8) {
                    biti -= 8;
                    *_blocks[blocki++]++ = (unsigned char)(bits >> biti);
                    if (blocki >= _nblocks)blocki = ++blockj == _nshort_data ? _nshort_blocks : 0;
                }
            }
            nbits = QR_INT_BITS;
        }
        j -= 2;
        if (j == 6)j--;
        l = j * stride;
        for (i = 0; i < stride; i++) {
            data1 = _data_bits[l + i];
            fp_mask1 = _fp_mask[l + i];
            data2 = _data_bits[l + i - stride];
            fp_mask2 = _fp_mask[l + i - stride];
            nbits = QR_MINI(_dim - (i << QR_INT_LOGBITS), QR_INT_BITS);
            while (nbits-- > 0) {
                if (!(fp_mask1 & 1)) {
                    bits = bits << 1 | data1 & 1;
                    biti++;
                }
                data1 >>= 1;
                fp_mask1 >>= 1;
                if (!(fp_mask2 & 1)) {
                    bits = bits << 1 | data2 & 1;
                    biti++;
                }
                data2 >>= 1;
                fp_mask2 >>= 1;
                if (biti >= 8) {
                    biti -= 8;
                    *_blocks[blocki++]++ = (unsigned char)(bits >> biti);
                    if (blocki >= _nblocks)blocki = ++blockj == _nshort_data ? _nshort_blocks : 0;
                }
            }
        }
    }
}

static void qr_sampling_grid_clear(qr_sampling_grid* _grid) {
    free(_grid->fpmask);
    free(_grid->cells[0]);
}

static unsigned rs_hgmul(const rs_gf256* _gf, unsigned _a, unsigned _logb) {
    return _a == 0 ? 0 : _gf->exp[_gf->log[_a] + _logb];
}

static void rs_calc_syndrome(const rs_gf256* _gf, int _m0,
    unsigned char* _s, int _npar, const unsigned char* _data, int _ndata) {
    int i;
    int j;
    for (j = 0; j < _npar; j++) {
        unsigned alphaj;
        unsigned sj;
        sj = 0;
        alphaj = _gf->log[_gf->exp[j + _m0]];
        for (i = 0; i < _ndata; i++)sj = _data[i] ^ rs_hgmul(_gf, sj, alphaj);
        _s[j] = sj;
    }
}

static void rs_poly_zero(unsigned char* _p, int _dp1) {
    memset(_p, 0, _dp1 * sizeof(*_p));
}

static void rs_init_lambda(const rs_gf256* _gf, unsigned char* _lambda, int _npar,
    const unsigned char* _erasures, int _nerasures, int _ndata) {
    int i;
    int j;
    rs_poly_zero(_lambda, (_npar < 4 ? 4 : _npar) + 1);
    _lambda[0] = 1;
    for (i = 0; i < _nerasures; i++)for (j = i + 1; j > 0; j--) {
        _lambda[j] ^= rs_hgmul(_gf, _lambda[j - 1], _ndata - 1 - _erasures[i]);
    }
}

static void rs_poly_copy(unsigned char* _p, const unsigned char* _q, int _dp1) {
    memcpy(_p, _q, _dp1 * sizeof(*_p));
}

static void rs_poly_mul_x(unsigned char* _p, const unsigned char* _q, int _dp1) {
    memmove(_p + 1, _q, (_dp1 - 1) * sizeof(*_p));
    _p[0] = 0;
}

static unsigned rs_gmul(const rs_gf256* _gf, unsigned _a, unsigned _b) {
    return _a == 0 || _b == 0 ? 0 : _gf->exp[_gf->log[_a] + _gf->log[_b]];
}

static void rs_poly_mult(const rs_gf256* _gf, unsigned char* _p, int _dp1,
    const unsigned char* _q, int _ep1, const unsigned char* _r, int _fp1) {
    int m;
    int i;
    rs_poly_zero(_p, _dp1);
    m = _ep1 < _dp1 ? _ep1 : _dp1;
    for (i = 0; i < m; i++)if (_q[i] != 0) {
        unsigned logqi;
        int      n;
        int      j;
        n = _dp1 - i < _fp1 ? _dp1 - i : _fp1;
        logqi = _gf->log[_q[i]];
        for (j = 0; j < n; j++)_p[i + j] ^= rs_hgmul(_gf, _r[j], logqi);
    }
}

static int rs_modified_berlekamp_massey(const rs_gf256* _gf,
    unsigned char* _lambda, const unsigned char* _s, unsigned char* _omega, int _npar,
    const unsigned char* _erasures, int _nerasures, int _ndata) {
    unsigned char tt[256];
    int           n;
    int           l;
    int           k;
    int           i;
    rs_init_lambda(_gf, _lambda, _npar, _erasures, _nerasures, _ndata);
    rs_poly_copy(tt, _lambda, _npar + 1);
    l = _nerasures;
    k = 0;
    for (n = _nerasures + 1; n <= _npar; n++) {
        unsigned d;
        rs_poly_mul_x(tt, tt, n - k + 1);
        d = 0;
        for (i = 0; i <= l; i++)d ^= rs_gmul(_gf, _lambda[i], _s[n - 1 - i]);
        if (d != 0) {
            unsigned logd;
            logd = _gf->log[d];
            if (l < n - k) {
                int t;
                for (i = 0; i <= n - k; i++) {
                    unsigned tti;
                    tti = tt[i];
                    tt[i] = rs_hgmul(_gf, _lambda[i], 255 - logd);
                    _lambda[i] = _lambda[i] ^ rs_hgmul(_gf, tti, logd);
                }
                t = n - k;
                k = n - l;
                l = t;
            }
            else for (i = 0; i <= l; i++)_lambda[i] = _lambda[i] ^ rs_hgmul(_gf, tt[i], logd);
        }
    }
    rs_poly_mult(_gf, _omega, _npar, _lambda, l + 1, _s, _npar);
    return l;
}

static unsigned rs_gsqrt(const rs_gf256* _gf, unsigned _a) {
    unsigned loga;
    if (!_a)return 0;
    loga = _gf->log[_a];
    return _gf->exp[loga + (255 & -(loga & 1)) >> 1];
}

static unsigned rs_gdiv(const rs_gf256* _gf, unsigned _a, unsigned _b) {
    return _a == 0 ? 0 : _gf->exp[_gf->log[_a] + 255 - _gf->log[_b]];
}

static int rs_quadratic_solve(const rs_gf256* _gf, unsigned _b, unsigned _c,
    unsigned char _x[2]) {
    unsigned b;
    unsigned logb;
    unsigned logb2;
    unsigned logb4;
    unsigned logb8;
    unsigned logb12;
    unsigned logb14;
    unsigned logc;
    unsigned logc2;
    unsigned logc4;
    unsigned c8;
    unsigned g3;
    unsigned z3;
    unsigned l3;
    unsigned c0;
    unsigned g2;
    unsigned l2;
    unsigned z2;
    int      inc;
    if (!_b) {
        _x[0] = rs_gsqrt(_gf, _c);
        return 1;
    }
    if (!_c) {
        _x[0] = 0;
        _x[1] = _b;
        return 2;
    }
    logb = _gf->log[_b];
    logc = _gf->log[_c];
    inc = logb % (255 / 15) == 0;
    if (inc) {
        b = _gf->exp[logb + 254];
        logb = _gf->log[b];
        _c = _gf->exp[logc + 253];
        logc = _gf->log[_c];
    }
    else b = _b;
    logb2 = _gf->log[_gf->exp[logb << 1]];
    logb4 = _gf->log[_gf->exp[logb2 << 1]];
    logb8 = _gf->log[_gf->exp[logb4 << 1]];
    logb12 = _gf->log[_gf->exp[logb4 + logb8]];
    logb14 = _gf->log[_gf->exp[logb2 + logb12]];
    logc2 = _gf->log[_gf->exp[logc << 1]];
    logc4 = _gf->log[_gf->exp[logc2 << 1]];
    c8 = _gf->exp[logc4 << 1];
    g3 = rs_hgmul(_gf,
        _gf->exp[logb14 + logc] ^ _gf->exp[logb12 + logc2] ^ _gf->exp[logb8 + logc4] ^ c8, logb);
    if (_gf->log[g3] % (255 / 15) != 0)return 0;
    z3 = rs_gdiv(_gf, g3, _gf->exp[logb8 << 1] ^ b);
    l3 = rs_hgmul(_gf, rs_gmul(_gf, z3, z3) ^ rs_hgmul(_gf, z3, logb) ^ _c, 255 - logb2);
    c0 = rs_hgmul(_gf, l3, 255 - 2 * (255 / 15));
    g2 = rs_hgmul(_gf,
        rs_hgmul(_gf, c0, 255 - 2 * (255 / 15)) ^ rs_gmul(_gf, c0, c0), 255 - 255 / 15);
    z2 = rs_gdiv(_gf, g2, _gf->exp[255 - (255 / 15) * 4] ^ _gf->exp[255 - (255 / 15)]);
    l2 = rs_hgmul(_gf,
        rs_gmul(_gf, z2, z2) ^ rs_hgmul(_gf, z2, 255 - (255 / 15)) ^ c0, 2 * (255 / 15));
    _x[0] = _gf->exp[_gf->log[z3 ^ rs_hgmul(_gf,
        rs_hgmul(_gf, l2, 255 / 3) ^ rs_hgmul(_gf, z2, 255 / 15), logb)] + inc];
    _x[1] = _x[0] ^ _b;
    return 2;
}


static int rs_cubic_solve(const rs_gf256* _gf,
    unsigned _a, unsigned _b, unsigned _c, unsigned char _x[3]) {
    unsigned k;
    unsigned logd;
    unsigned d2;
    unsigned logd2;
    unsigned logw;
    int      nroots;
    if (!_c) {
        nroots = rs_quadratic_solve(_gf, _a, _b, _x);
        if (_b)_x[nroots++] = 0;
        return nroots;
    }
    k = rs_gmul(_gf, _a, _b) ^ _c;
    d2 = rs_gmul(_gf, _a, _a) ^ _b;
    if (!d2) {
        int logx;
        if (!k) {
            _x[0] = _a;
            return 1;
        }
        logx = _gf->log[k];
        if (logx % 3 != 0)return 0;
        logx /= 3;
        _x[0] = _a ^ _gf->exp[logx];
        _x[1] = _a ^ _gf->exp[logx + 255 / 3];
        _x[2] = _a ^ _x[0] ^ _x[1];
        return 3;
    }
    logd2 = _gf->log[d2];
    logd = logd2 + (255 & -(logd2 & 1)) >> 1;
    k = rs_gdiv(_gf, k, _gf->exp[logd + logd2]);
    nroots = rs_quadratic_solve(_gf, k, 1, _x);
    if (nroots < 1) {
        return 0;
    }
    logw = _gf->log[_x[0]];
    if (logw) {
        if (logw % 3 != 0)return 0;
        logw /= 3;
        /*Recover x from w.*/
        _x[0] = _gf->exp[_gf->log[_gf->exp[logw] ^ _gf->exp[255 - logw]] + logd] ^ _a;
        logw += 255 / 3;
        _x[1] = _gf->exp[_gf->log[_gf->exp[logw] ^ _gf->exp[255 - logw]] + logd] ^ _a;
        _x[2] = _x[0] ^ _x[1] ^ _a;
        return 3;
    }
    else {
        _x[0] = _a;
        /*In this case _x[1] is a double root, so we know the Reed-Solomon code is
           invalid.
          Note that we still have to return at least one root, because if we're
           being called by the quartic solver, the quartic might still have 4
           distinct roots.
          But we don't need more than one root, so we can avoid computing the
           expensive one.*/
           /*_x[1]=_gf->exp[_gf->log[_gf->exp[255/3]^_gf->exp[2*(255/3)]]+logd]^_a;*/
        return 1;
    }
}

static int rs_quartic_solve(const rs_gf256* _gf,
    unsigned _a, unsigned _b, unsigned _c, unsigned _d, unsigned char _x[3]) {
    unsigned r;
    unsigned s;
    unsigned t;
    unsigned b;
    int      nroots;
    int      i;
    if (!_d) {
        nroots = rs_cubic_solve(_gf, _a, _b, _c, _x);
        if (_c)_x[nroots++] = 0;
        return nroots;
    }
    if (_a) {
        unsigned loga;
        loga = _gf->log[_a];
        r = rs_hgmul(_gf, _c, 255 - loga);
        s = rs_gsqrt(_gf, r);
        t = _d ^ rs_gmul(_gf, _b, r) ^ rs_gmul(_gf, r, r);
        if (t) {
            unsigned logti;
            logti = 255 - _gf->log[t];
            nroots = rs_quartic_solve(_gf, 0, rs_hgmul(_gf, _b ^ rs_hgmul(_gf, s, loga), logti),
                _gf->exp[loga + logti], _gf->exp[logti], _x);
            for (i = 0; i < nroots; i++)_x[i] = _gf->exp[255 - _gf->log[_x[i]]] ^ s;
        }
        else {
            nroots = rs_quadratic_solve(_gf, _a, _b ^ r, _x);
            if (nroots != 2 || _x[0] != s && _x[1] != s)_x[nroots++] = s;
        }
        return nroots;
    }
    /*If there are no odd powers, it's really just a quadratic in disguise.*/
    if (!_c)return rs_quadratic_solve(_gf, rs_gsqrt(_gf, _b), rs_gsqrt(_gf, _d), _x);
    /*Factor into (x**2 + r*x + s)*(x**2 + r*x + t) by solving for r, which can
       be shown to satisfy r**3 + _b*r + _c == 0.*/
    nroots = rs_cubic_solve(_gf, 0, _b, _c, _x);
    if (nroots < 1) {
        /*The Reed-Solomon code is only valid if we can find 4 distinct roots in
           GF(2**8).
          If the cubic does not factor into 3 (possibly duplicate) roots, then we
           know that the quartic must have a non-trivial irreducible factor.*/
        return 0;
    }
    r = _x[0];
    /*Now solve for s and t.*/
    b = rs_gdiv(_gf, _c, r);
    nroots = rs_quadratic_solve(_gf, b, _d, _x);
    if (nroots < 2)return 0;
    s = _x[0];
    t = _x[1];
    /*_c=r*(s^t) was non-zero, so s and t must be distinct.
      But if z is a root of z**2 ^ r*z ^ s, then so is (z^r), and s = z*(z^r).
      Hence if z is also a root of z**2 ^ r*z ^ t, then t = s, a contradiction.
      Thus all four roots are distinct, if they exist.*/
    nroots = rs_quadratic_solve(_gf, r, s, _x);
    return nroots + rs_quadratic_solve(_gf, r, t, _x + nroots);
}


static int rs_find_roots(const rs_gf256* _gf, unsigned char* _epos,
    const unsigned char* _lambda, int _nerrors, int _ndata) {
    unsigned alpha;
    int      nroots;
    int      i;
    nroots = 0;
    if (_nerrors <= 4) {
        _nerrors = rs_quartic_solve(_gf, _lambda[1], _lambda[2], _lambda[3], _lambda[4],
            _epos);
        for (i = 0; i < _nerrors; i++)if (_epos[i]) {
            alpha = _gf->log[_epos[i]];
            if ((int)alpha < _ndata)_epos[nroots++] = alpha;
        }
        return nroots;
    }
    else for (alpha = 0; (int)alpha < _ndata; alpha++) {
        unsigned alphai;
        unsigned sum;
        sum = 0;
        alphai = 0;
        for (i = 0; i <= _nerrors; i++) {
            sum ^= rs_hgmul(_gf, _lambda[_nerrors - i], alphai);
            alphai = _gf->log[_gf->exp[alphai + alpha]];
        }
        if (!sum)_epos[nroots++] = alpha;
    }
    return nroots;
}


int rs_correct(const rs_gf256* _gf, int _m0, unsigned char* _data, int _ndata,
    int _npar, const unsigned char* _erasures, int _nerasures) {
    unsigned char lambda[256];
    unsigned char omega[256];
    unsigned char epos[256];
    unsigned char s[256];
    int           i;
    if (_nerasures > _npar)return -1;
    rs_calc_syndrome(_gf, _m0, s, _npar, _data, _ndata);
    for (i = 0; i < _npar; i++)if (s[i]) {
        int nerrors;
        int j;
        nerrors = rs_modified_berlekamp_massey(_gf, lambda, s, omega, _npar,
            _erasures, _nerasures, _ndata);
        if (nerrors <= 0 || nerrors - _nerasures > _npar - _nerasures >> 1)return -1;
        if (rs_find_roots(_gf, epos, lambda, nerrors, _ndata) < nerrors)return -1;
        /*Now compute the error magnitudes.*/
        for (i = 0; i < nerrors; i++) {
            unsigned a;
            unsigned b;
            unsigned alpha;
            unsigned alphan1;
            unsigned alphan2;
            unsigned alphanj;
            alpha = epos[i];
            /*Evaluate omega at alpha**-1.*/
            a = 0;
            alphan1 = 255 - alpha;
            alphanj = 0;
            for (j = 0; j < _npar; j++) {
                a ^= rs_hgmul(_gf, omega[j], alphanj);
                alphanj = _gf->log[_gf->exp[alphanj + alphan1]];
            }
            b = 0;
            alphan2 = _gf->log[_gf->exp[alphan1 << 1]];
            alphanj = alphan1 + _m0 * alpha % 255;
            for (j = 1; j <= _npar; j += 2) {
                b ^= rs_hgmul(_gf, lambda[j], alphanj);
                alphanj = _gf->log[_gf->exp[alphanj + alphan2]];
            }
            _data[_ndata - 1 - alpha] ^= rs_gdiv(_gf, a, b);
        }
        return nerrors;
    }
    return 0;
}

static void qr_pack_buf_init(qr_pack_buf* _b,
    const unsigned char* _data, int _ndata) {
    _b->buf = _data;
    _b->storage = _ndata;
    _b->endbyte = _b->endbit = 0;
}

static int qr_pack_buf_avail(const qr_pack_buf* _b) {
    return (_b->storage - _b->endbyte << 3) - _b->endbit;
}

static int qr_pack_buf_read(qr_pack_buf* _b, int _bits) {
    const unsigned char* p;
    unsigned             ret;
    int                  m;
    int                  d;
    m = 16 - _bits;
    _bits += _b->endbit;
    d = _b->storage - _b->endbyte;
    if (d <= 2) {
        if (d * 8 < _bits) {
            _b->endbyte += _bits >> 3;
            _b->endbit = _bits & 7;
            return -1;
        }
        else if (!_bits)return 0;
    }
    p = _b->buf + _b->endbyte;
    ret = p[0] << 8 + _b->endbit;
    if (_bits > 8) {
        ret |= p[1] << _b->endbit;
        if (_bits > 16)ret |= p[2] >> 8 - _b->endbit;
    }
    _b->endbyte += _bits >> 3;
    _b->endbit = _bits & 7;
    return (ret & 0xFFFF) >> m;
}

static const unsigned char QR_ALNUM_TABLE[45] = {
  '0','1','2','3','4','5','6','7','8','9',
  'A','B','C','D','E','F','G','H','I','J',
  'K','L','M','N','O','P','Q','R','S','T',
  'U','V','W','X','Y','Z',' ','$','%','*',
  '+','-','.','/',':'
};

static const unsigned char LEN_BITS[3][4] = {
              {10, 9, 8, 8},
              {12,11,16,10},
              {14,13,16,12}
};

static int qr_code_data_parse(qr_code_data* _qrdata, int _version,
    const unsigned char* _data, int _ndata) {
    qr_pack_buf qpb;
    unsigned    self_parity;
    int         centries;
    int         len_bits_idx;
    _qrdata->entries = nullptr;
    _qrdata->nentries = 0;
    _qrdata->sa_size = 0;
    self_parity = 0;
    centries = 0;
    len_bits_idx = (_version > 9) + (_version > 26);
    qr_pack_buf_init(&qpb, _data, _ndata);
    while (qr_pack_buf_avail(&qpb) >= 4) {
        qr_code_data_entry* entry;
        int                 mode;
        mode = qr_pack_buf_read(&qpb, 4);
        if (!mode)break;
        if (_qrdata->nentries >= centries) {
            centries = centries << 1 | 1;
            _qrdata->entries = (qr_code_data_entry*)realloc(_qrdata->entries,
                centries * sizeof(*_qrdata->entries));
        }
        entry = _qrdata->entries + _qrdata->nentries++;
        entry->mode = (qr_mode)mode;
        entry->payload.data.buf = nullptr;
        switch (mode) {
            
        case QR_MODE_NUM: {
            unsigned char* buf;
            unsigned       bits;
            unsigned       c;
            int            len;
            int            count;
            int            rem;
            len = qr_pack_buf_read(&qpb, LEN_BITS[len_bits_idx][0]);
            if (len < 0)return -1;
            count = len / 3;
            rem = len % 3;
            if (qr_pack_buf_avail(&qpb) < 10 * count + 7 * (rem >> 1 & 1) + 4 * (rem & 1))return -1;
            entry->payload.data.buf = buf = (unsigned char*)malloc(len * sizeof(*buf));
            entry->payload.data.len = len;
            while (count-- > 0) {
                bits = qr_pack_buf_read(&qpb, 10);
                if (bits >= 1000)return -1;
                c = '0' + bits / 100;
                self_parity ^= c;
                *buf++ = (unsigned char)c;
                bits %= 100;
                c = '0' + bits / 10;
                self_parity ^= c;
                *buf++ = (unsigned char)c;
                c = '0' + bits % 10;
                self_parity ^= c;
                *buf++ = (unsigned char)c;
            }
            if (rem > 1) {
                bits = qr_pack_buf_read(&qpb, 7);
                if (bits >= 100)return -1;
                c = '0' + bits / 10;
                self_parity ^= c;
                *buf++ = (unsigned char)c;
                c = '0' + bits % 10;
                self_parity ^= c;
                *buf++ = (unsigned char)c;
            }
            else if (rem) {
                bits = qr_pack_buf_read(&qpb, 4);
                if (bits >= 10)return -1;
                c = '0' + bits;
                self_parity ^= c;
                *buf++ = (unsigned char)c;
            }
        }break;
        case QR_MODE_ALNUM: {
            unsigned char* buf;
            unsigned       bits;
            unsigned       c;
            int            len;
            int            count;
            int            rem;
            len = qr_pack_buf_read(&qpb, LEN_BITS[len_bits_idx][1]);
            if (len < 0)return -1;
            count = len >> 1;
            rem = len & 1;
            if (qr_pack_buf_avail(&qpb) < 11 * count + 6 * rem)return -1;
            entry->payload.data.buf = buf = (unsigned char*)malloc(len * sizeof(*buf));
            entry->payload.data.len = len;
            while (count-- > 0) {
                bits = qr_pack_buf_read(&qpb, 11);
                if (bits >= 2025)return -1;
                c = QR_ALNUM_TABLE[bits / 45];
                self_parity ^= c;
                *buf++ = (unsigned char)c;
                c = QR_ALNUM_TABLE[bits % 45];
                self_parity ^= c;
                *buf++ = (unsigned char)c;
                len -= 2;
            }
            if (rem) {
                bits = qr_pack_buf_read(&qpb, 6);
                if (bits >= 45)return -1;
                c = QR_ALNUM_TABLE[bits];
                self_parity ^= c;
                *buf++ = (unsigned char)c;
            }
        }break;
        case QR_MODE_STRUCT: {
            int bits;
            bits = qr_pack_buf_read(&qpb, 16);
            if (bits < 0)return -1;
            if (_qrdata->sa_size == 0) {
                _qrdata->sa_index = entry->payload.sa.sa_index =
                    (unsigned char)(bits >> 12 & 0xF);
                _qrdata->sa_size = entry->payload.sa.sa_size =
                    (unsigned char)((bits >> 8 & 0xF) + 1);
                _qrdata->sa_parity = entry->payload.sa.sa_parity =
                    (unsigned char)(bits & 0xFF);
            }
        }break;
        case QR_MODE_BYTE: {
            unsigned char* buf;
            unsigned       c;
            int            len;
            len = qr_pack_buf_read(&qpb, LEN_BITS[len_bits_idx][2]);
            if (len < 0)return -1;
            if (qr_pack_buf_avail(&qpb) < len << 3)return -1;
            entry->payload.data.buf = buf = (unsigned char*)malloc(len * sizeof(*buf));
            entry->payload.data.len = len;
            while (len-- > 0) {
                c = qr_pack_buf_read(&qpb, 8);
                self_parity ^= c;
                *buf++ = (unsigned char)c;
            }
        }break;
        case QR_MODE_FNC1_1ST:break;
        case QR_MODE_ECI: {
            unsigned val;
            int      bits;
            bits = qr_pack_buf_read(&qpb, 8);
            if (bits < 0)return -1;
            if (!(bits & 0x80))val = bits;
            else if (!(bits & 0x40)) {
                val = bits & 0x3F << 8;
                bits = qr_pack_buf_read(&qpb, 8);
                if (bits < 0)return -1;
                val |= bits;
            }
            else if (!(bits & 0x20)) {
                val = bits & 0x1F << 16;
                bits = qr_pack_buf_read(&qpb, 16);
                if (bits < 0)return -1;
                val |= bits;
                if (val >= 1000000)return -1;
            }
            else return -1;
            entry->payload.eci = val;
        }break;
        case QR_MODE_KANJI: {
            unsigned char* buf;
            unsigned       bits;
            int            len;
            len = qr_pack_buf_read(&qpb, LEN_BITS[len_bits_idx][3]);
            if (len < 0)return -1;
            if (qr_pack_buf_avail(&qpb) < 13 * len)return -1;
            entry->payload.data.buf = buf = (unsigned char*)malloc(2 * len * sizeof(*buf));
            entry->payload.data.len = 2 * len;
            while (len-- > 0) {
                bits = qr_pack_buf_read(&qpb, 13);
                bits = (bits / 0xC0 << 8 | bits % 0xC0) + 0x8140;
                if (bits >= 0xA000)bits += 0x4000;
                self_parity ^= bits;
                *buf++ = (unsigned char)(bits >> 8);
                *buf++ = (unsigned char)(bits & 0xFF);
            }
        }break;
        case QR_MODE_FNC1_2ND: {
            int bits;
            bits = qr_pack_buf_read(&qpb, 8);
            if (!(bits >= 0 && bits < 100 || bits >= 165 && bits < 191 || bits >= 197 && bits < 223)) {
                return -1;
            }
            entry->payload.ai = bits;
        }break;
        default: {
            return -1;
        }break;
        }
    }
    _qrdata->self_parity = ((self_parity >> 8) ^ self_parity) & 0xFF;
    _qrdata->entries = (qr_code_data_entry*)realloc(_qrdata->entries,
        _qrdata->nentries * sizeof(*_qrdata->entries));
    return 0;
}

static void qr_code_data_clear(qr_code_data* _qrdata) {
    int i;
    for (i = 0; i < _qrdata->nentries; i++) {
        if (QR_MODE_HAS_DATA(_qrdata->entries[i].mode)) {
            free(_qrdata->entries[i].payload.data.buf);
        }
    }
    free(_qrdata->entries);
}

static int qr_code_decode(qr_code_data* _qrdata, const rs_gf256* _gf,
    const qr_point _ul_pos, const qr_point _ur_pos, const qr_point _dl_pos,
    int _version, int _fmt_info,
    const unsigned char* _img, int _width, int _height) {
    qr_sampling_grid   grid;
    unsigned* data_bits;
    unsigned char** blocks;
    unsigned char* block_data;
    int                nblocks;
    int                nshort_blocks;
    int                ncodewords;
    int                block_sz;
    int                ecc_level;
    int                ndata;
    int                npar;
    int                dim;
    int                ret;
    int                i;
    qr_sampling_grid_init(&grid, _version, _ul_pos, _ur_pos, _dl_pos, _qrdata->bbox,
        _img, _width, _height);
    dim = 17 + (_version << 2);
    data_bits = (unsigned*)malloc(
        dim * (dim + QR_INT_BITS - 1 >> QR_INT_LOGBITS) * sizeof(*data_bits));
    qr_sampling_grid_sample(&grid, data_bits, dim, _fmt_info, _img, _width, _height);
    ecc_level = (_fmt_info >> 3) ^ 1;
    nblocks = QR_RS_NBLOCKS[_version - 1][ecc_level];
    npar = *(QR_RS_NPAR_VALS + QR_RS_NPAR_OFFS[_version - 1] + ecc_level);
    ncodewords = qr_code_ncodewords(_version);
    block_sz = ncodewords / nblocks;
    nshort_blocks = nblocks - (ncodewords % nblocks);
    blocks = (unsigned char**)malloc(nblocks * sizeof(*blocks));
    block_data = (unsigned char*)malloc(ncodewords * sizeof(*block_data));
    blocks[0] = block_data;
    for (i = 1; i < nblocks; i++)blocks[i] = blocks[i - 1] + block_sz + (i > nshort_blocks);
    qr_samples_unpack(blocks, nblocks, block_sz - npar, nshort_blocks,
        data_bits, grid.fpmask, dim);
    qr_sampling_grid_clear(&grid);
    free(blocks);
    free(data_bits);
    ndata = 0;
    ncodewords = 0;
    ret = 0;
    for (i = 0; i < nblocks; i++) {
        int block_szi;
        int ndatai;
        block_szi = block_sz + (i >= nshort_blocks);
        ret = rs_correct(_gf, QR_M0, block_data + ncodewords, block_szi, npar, nullptr, 0);
        if (ret<0 || _version == 1 && ret>ecc_level + 1 << 1 ||
            _version == 2 && ecc_level == 0 && ret > 4) {
            ret = -1;
            break;
        }
        ndatai = block_szi - npar;
        memmove(block_data + ndata, block_data + ncodewords, ndatai * sizeof(*block_data));
        ncodewords += block_szi;
        ndata += ndatai;
    }
    if (ret >= 0) {
        ret = qr_code_data_parse(_qrdata, _version, block_data, ndata);
        if (ret < 0)qr_code_data_clear(_qrdata);
        _qrdata->version = _version;
        _qrdata->ecc_level = ecc_level;
    }
    free(block_data);
    return ret;
}


static int qr_reader_try_configuration(qr_reader* _reader,
    qr_code_data* _qrdata, const unsigned char* _img, int _width, int _height,
    qr_finder_center* _c[3]) {
    int      ci[7];
    unsigned maxd;
    int      ccw;
    int      i0;
    int      i;
    ccw = qr_point_ccw(_c[0]->pos, _c[1]->pos, _c[2]->pos);
    if (!ccw)return -1;
    ci[6] = ci[3] = ci[0] = 0;
    ci[4] = ci[1] = 1 + (ccw < 0);
    ci[5] = ci[2] = 2 - (ccw < 0);
    maxd = qr_point_distance2(_c[1]->pos, _c[2]->pos);
    i0 = 0;
    for (i = 1; i < 3; i++) {
        unsigned d;
        d = qr_point_distance2(_c[ci[i + 1]]->pos, _c[ci[i + 2]]->pos);
        if (d > maxd) {
            i0 = i;
            maxd = d;
        }
    }

    for (i = i0; i < i0 + 3; i++) {
        qr_aff    aff;
        qr_hom    hom;
        qr_finder ul;
        qr_finder ur;
        qr_finder dl;
        qr_point  bbox[4];
        int       res;
        int       ur_version;
        int       dl_version;
        int       fmt_info;
        ul.c = _c[ci[i]];
        ur.c = _c[ci[i + 1]];
        dl.c = _c[ci[i + 2]];
        res = QR_INT_BITS - 2 - QR_FINDER_SUBPREC - qr_ilog(QR_MAXI(_width, _height) - 1);
        qr_aff_init(&aff, ul.c->pos, ur.c->pos, dl.c->pos, res);
        qr_aff_unproject(ur.o, &aff, ur.c->pos[0], ur.c->pos[1]);
        qr_finder_edge_pts_aff_classify(&ur, &aff);
        if (qr_finder_estimate_module_size_and_version(&ur, 1 << res, 1 << res) < 0)continue;
        qr_aff_unproject(dl.o, &aff, dl.c->pos[0], dl.c->pos[1]);
        qr_finder_edge_pts_aff_classify(&dl, &aff);
        if (qr_finder_estimate_module_size_and_version(&dl, 1 << res, 1 << res) < 0)continue;
        if (abs(ur.eversion[1] - dl.eversion[0]) > QR_LARGE_VERSION_SLACK)continue;
        qr_aff_unproject(ul.o, &aff, ul.c->pos[0], ul.c->pos[1]);
        qr_finder_edge_pts_aff_classify(&ul, &aff);
        if (qr_finder_estimate_module_size_and_version(&ul, 1 << res, 1 << res) < 0 ||
            abs(ul.eversion[1] - ur.eversion[1]) > QR_LARGE_VERSION_SLACK ||
            abs(ul.eversion[0] - dl.eversion[0]) > QR_LARGE_VERSION_SLACK) {
            continue;
        }
        if (qr_hom_fit(&hom, &ul, &ur, &dl, bbox, &aff,
            &_reader->isaac, _img, _width, _height) < 0) {
            continue;
        }
        memcpy(_qrdata->bbox, bbox, sizeof(bbox));
        qr_hom_unproject(ul.o, &hom, ul.c->pos[0], ul.c->pos[1]);
        qr_hom_unproject(ur.o, &hom, ur.c->pos[0], ur.c->pos[1]);
        qr_hom_unproject(dl.o, &hom, dl.c->pos[0], dl.c->pos[1]);
        qr_finder_edge_pts_hom_classify(&ur, &hom);
        if (qr_finder_estimate_module_size_and_version(&ur,
            ur.o[0] - ul.o[0], ur.o[0] - ul.o[0]) < 0) {
            continue;
        }
        qr_finder_edge_pts_hom_classify(&dl, &hom);
        if (qr_finder_estimate_module_size_and_version(&dl,
            dl.o[1] - ul.o[1], dl.o[1] - ul.o[1]) < 0) {
            continue;
        }
        if (ur.eversion[1] == dl.eversion[0] && ur.eversion[1] < 7) {
            ur_version = ur.eversion[1];
        }
        else {
            if (abs(ur.eversion[1] - dl.eversion[0]) > QR_LARGE_VERSION_SLACK)continue;

            if (ur.eversion[1] >= 7 - QR_LARGE_VERSION_SLACK) {
                ur_version = qr_finder_version_decode(&ur, &hom, _img, _width, _height, 0);
                if (abs(ur_version - ur.eversion[1]) > QR_LARGE_VERSION_SLACK)ur_version = -1;
            }
            else ur_version = -1;
            if (dl.eversion[0] >= 7 - QR_LARGE_VERSION_SLACK) {
                dl_version = qr_finder_version_decode(&dl, &hom, _img, _width, _height, 1);
                if (abs(dl_version - dl.eversion[0]) > QR_LARGE_VERSION_SLACK)dl_version = -1;
            }
            else dl_version = -1;
            if (ur_version >= 0) {
                if (dl_version >= 0 && dl_version != ur_version)continue;
            }
            else if (dl_version < 0)continue;
            else ur_version = dl_version;
        }
        qr_finder_edge_pts_hom_classify(&ul, &hom);
        if (qr_finder_estimate_module_size_and_version(&ul,
            ur.o[0] - dl.o[0], dl.o[1] - ul.o[1]) < 0 ||
            abs(ul.eversion[1] - ur.eversion[1]) > QR_SMALL_VERSION_SLACK ||
            abs(ul.eversion[0] - dl.eversion[0]) > QR_SMALL_VERSION_SLACK) {
            continue;
        }
        fmt_info = qr_finder_fmt_info_decode(&ul, &ur, &dl, &hom, _img, _width, _height);
        if (fmt_info < 0 ||
            qr_code_decode(_qrdata, &_reader->gf, ul.c->pos, ur.c->pos, dl.c->pos,
                ur_version, fmt_info, _img, _width, _height) < 0) {
            QR_SWAP2I(hom.inv[0][0], hom.inv[1][0]);
            QR_SWAP2I(hom.inv[0][1], hom.inv[1][1]);
            QR_SWAP2I(hom.fwd[0][0], hom.fwd[0][1]);
            QR_SWAP2I(hom.fwd[1][0], hom.fwd[1][1]);
            QR_SWAP2I(hom.fwd[2][0], hom.fwd[2][1]);
            QR_SWAP2I(ul.o[0], ul.o[1]);
            QR_SWAP2I(ul.size[0], ul.size[1]);
            QR_SWAP2I(ur.o[0], ur.o[1]);
            QR_SWAP2I(ur.size[0], ur.size[1]);
            QR_SWAP2I(dl.o[0], dl.o[1]);
            QR_SWAP2I(dl.size[0], dl.size[1]);
            fmt_info = qr_finder_fmt_info_decode(&ul, &dl, &ur, &hom, _img, _width, _height);
            if (fmt_info < 0)continue;
            QR_SWAP2I(bbox[1][0], bbox[2][0]);
            QR_SWAP2I(bbox[1][1], bbox[2][1]);
            memcpy(_qrdata->bbox, bbox, sizeof(bbox));
            if (qr_code_decode(_qrdata, &_reader->gf, ul.c->pos, dl.c->pos, ur.c->pos,
                ur_version, fmt_info, _img, _width, _height) < 0) {
                continue;
            }
        }
        return ur_version;
    }
    return -1;
}

static void qr_code_data_list_add(qr_code_data_list* _qrlist,
    qr_code_data* _qrdata) {
    if (_qrlist->nqrdata >= _qrlist->cqrdata) {
        _qrlist->cqrdata = _qrlist->cqrdata << 1 | 1;
        _qrlist->qrdata = (qr_code_data*)realloc(_qrlist->qrdata,
            _qrlist->cqrdata * sizeof(*_qrlist->qrdata));
    }
    memcpy(_qrlist->qrdata + _qrlist->nqrdata++, _qrdata, sizeof(*_qrdata));
}

void qr_reader_match_centers(qr_reader* _reader, qr_code_data_list* _qrlist,
    qr_finder_center* _centers, int _ncenters,
    const unsigned char* _img, int _width, int _height) {
    unsigned char* mark;
    int            nfailures_max;
    int            nfailures;
    int            i;
    int            j;
    int            k;
    mark = (unsigned char*)calloc(_ncenters, sizeof(*mark));
    nfailures_max = QR_MAXI(8192, _width * _height >> 9);
    nfailures = 0;
    for (i = 0; i < _ncenters; i++) {
        for (j = i + 1; !mark[i] && j < _ncenters; j++) {
            for (k = j + 1; !mark[j] && k < _ncenters; k++)if (!mark[k]) {
                qr_finder_center* c[3];
                qr_code_data      qrdata;
                int               version;
                c[0] = _centers + i;
                c[1] = _centers + j;
                c[2] = _centers + k;
                version = qr_reader_try_configuration(_reader, &qrdata,
                    _img, _width, _height, c);
                if (version >= 0) {
                    int ninside;
                    int l;
                    qr_code_data_list_add(_qrlist, &qrdata);
                    for (l = 0; l < 4; l++) {
                        _qrlist->qrdata[_qrlist->nqrdata - 1].bbox[l][0] >>= QR_FINDER_SUBPREC;
                        _qrlist->qrdata[_qrlist->nqrdata - 1].bbox[l][1] >>= QR_FINDER_SUBPREC;
                    }
                    mark[i] = mark[j] = mark[k] = 1;
                    for (l = ninside = 0; l < _ncenters; l++)if (!mark[l]) {
                        if (qr_point_ccw(qrdata.bbox[0], qrdata.bbox[1], _centers[l].pos) >= 0 &&
                            qr_point_ccw(qrdata.bbox[1], qrdata.bbox[3], _centers[l].pos) >= 0 &&
                            qr_point_ccw(qrdata.bbox[3], qrdata.bbox[2], _centers[l].pos) >= 0 &&
                            qr_point_ccw(qrdata.bbox[2], qrdata.bbox[0], _centers[l].pos) >= 0) {
                            mark[l] = 2;
                            ninside++;
                        }
                    }
                    if (ninside >= 3) {
                        qr_finder_center* inside;
                        inside = (qr_finder_center*)malloc(ninside * sizeof(*inside));
                        for (l = ninside = 0; l < _ncenters; l++) {
                            if (mark[l] == 2)*&inside[ninside++] = *&_centers[l];
                        }
                        qr_reader_match_centers(_reader, _qrlist, inside, ninside,
                            _img, _width, _height);
                        free(inside);
                    }
                    for (l = 0; l < _ncenters; l++)if (mark[l] == 2)mark[l] = 1;
                    nfailures = 0;
                }
                else if (++nfailures > nfailures_max) {
                    i = j = k = _ncenters;
                }
            }
        }
    }
    free(mark);
}

static __inline void sym_add_point(zbar_symbol_t* sym,
    int x,
    int y)
{
    int i = sym->npts;
    if (++sym->npts >= sym->pts_alloc)
        sym->pts = (point_t*)realloc(sym->pts, ++sym->pts_alloc * sizeof(point_t));
    sym->pts[i].x = x;
    sym->pts[i].y = y;
}

static void enc_list_mtf(iconv_t _enc_list[3], iconv_t _enc) {
    int i;
    for (i = 0; i < 3; i++)if (_enc_list[i] == _enc) {
        int j;
        for (j = i; j-- > 0;)_enc_list[j + 1] = _enc_list[j];
        _enc_list[0] = _enc;
        break;
    }
}

static int text_is_ascii(const unsigned char* _text, int _len) {
    int i;
    for (i = 0; i < _len; i++)if (_text[i] >= 0x80)return 0;
    return 1;
}

static int text_is_latin1(const unsigned char* _text, int _len) {
    int i;
    for (i = 0; i < _len; i++) {
        if (_text[i] >= 0x80 && _text[i] < 0xA0)return 0;
    }
    return 1;
}

int qr_code_data_list_extract_text(const qr_code_data_list* _qrlist,
    zbar_image_scanner_t* iscn,
    zbar_image_t* img)
{
    iconv_t              sjis_cd;
    iconv_t              utf8_cd;
    iconv_t              latin1_cd;
    const qr_code_data* qrdata;
    int                  nqrdata;
    unsigned char* mark;
    int                  ntext;
    int                  i;
    qrdata = _qrlist->qrdata;
    nqrdata = _qrlist->nqrdata;
    mark = (unsigned char*)calloc(nqrdata, sizeof(*mark));
    ntext = 0;
    latin1_cd = iconv_open("UTF-8", "ISO8859-1");
    sjis_cd = iconv_open("UTF-8", "SJIS");
    utf8_cd = iconv_open("UTF-8", "UTF-8");
    for (i = 0; i < nqrdata; i++)if (!mark[i]) {
        const qr_code_data* qrdataj;
        const qr_code_data_entry* entry;
        iconv_t                   enc_list[3];
        iconv_t                   eci_cd;
        int                       sa[16];
        int                       sa_size;
        char* sa_text;
        size_t                    sa_ntext;
        size_t                    sa_ctext;
        int                       fnc1;
        int                       fnc1_2ai;
        int                       has_kanji;
        int                       eci;
        int                       err;
        int                       j;
        int                       k;
        zbar_symbol_t* syms = nullptr, ** sym = &syms;
        qr_point dir;
        int horiz;

        if (qrdata[i].sa_size) {
            unsigned sa_parity;
            sa_size = qrdata[i].sa_size;
            sa_parity = qrdata[i].sa_parity;
            for (j = 0; j < sa_size; j++)sa[j] = -1;
            for (j = i; j < nqrdata; j++)if (!mark[j]) {
                if (qrdata[j].sa_size == sa_size && qrdata[j].sa_parity == sa_parity &&
                    sa[qrdata[j].sa_index] < 0) {
                    sa[qrdata[j].sa_index] = j;
                    mark[j] = 1;
                }
            }
        }
        else {
            sa[0] = i;
            sa_size = 1;
        }

        sa_ctext = 0;
        fnc1 = 0;
        fnc1_2ai = 0;
        has_kanji = 0;
        for (j = 0; j < sa_size; j++)if (sa[j] >= 0) {
            qrdataj = qrdata + sa[j];
            for (k = 0; k < qrdataj->nentries; k++) {
                int shift;
                entry = qrdataj->entries + k;
                shift = 0;
                switch (entry->mode) {
                case QR_MODE_FNC1_1ST: {
                    if (!fnc1)fnc1 = MOD(ZBAR_MOD_GS1);
                }break;
                case QR_MODE_FNC1_2ND: {
                    if (!fnc1) {
                        fnc1 = MOD(ZBAR_MOD_AIM);
                        fnc1_2ai = entry->payload.ai;
                        sa_ctext += 2;
                    }
                }break;
                case QR_MODE_KANJI:has_kanji = 1;
                case QR_MODE_BYTE:shift = 2;
                default: {
                    if (QR_MODE_HAS_DATA(entry->mode)) {
                        sa_ctext += entry->payload.data.len << shift;
                    }
                }break;
                }
            }
        }

        sa_text = (char*)malloc((sa_ctext + 1) * sizeof(*sa_text));
        sa_ntext = 0;
        if (fnc1 == MOD(ZBAR_MOD_AIM)) {
            if (fnc1_2ai < 100) {
                sa_text[sa_ntext++] = '0' + fnc1_2ai / 10;
                sa_text[sa_ntext++] = '0' + fnc1_2ai % 10;
            }
            else sa_text[sa_ntext++] = (char)(fnc1_2ai - 100);
        }
        eci = -1;
        enc_list[0] = sjis_cd;
        enc_list[1] = latin1_cd;
        enc_list[2] = utf8_cd;
        eci_cd = (iconv_t)-1;
        err = 0;
        for (j = 0; j < sa_size && !err; j++, sym = &(*sym)->next) {
            *sym = _zbar_image_scanner_alloc_sym(iscn, ZBAR_QRCODE, 0);
            (*sym)->datalen = sa_ntext;
            if (sa[j] < 0) {
                (*sym)->type = ZBAR_PARTIAL;
                for (j++; j < sa_size && sa[j] < 0; j++);
                if (j >= sa_size)break;
                sa_text[sa_ntext++] = '\0';
                (*sym)->datalen = sa_ntext;
                sym = &(*sym)->next;
                *sym = _zbar_image_scanner_alloc_sym(iscn, ZBAR_QRCODE, 0);
            }

            qrdataj = qrdata + sa[j];
            sym_add_point(*sym, qrdataj->bbox[0][0], qrdataj->bbox[0][1]);
            sym_add_point(*sym, qrdataj->bbox[2][0], qrdataj->bbox[2][1]);
            sym_add_point(*sym, qrdataj->bbox[3][0], qrdataj->bbox[3][1]);
            sym_add_point(*sym, qrdataj->bbox[1][0], qrdataj->bbox[1][1]);

            dir[0] = (qrdataj->bbox[0][0] - qrdataj->bbox[2][0] +
                qrdataj->bbox[1][0] - qrdataj->bbox[3][0]);
            dir[1] = (qrdataj->bbox[2][1] - qrdataj->bbox[0][1] +
                qrdataj->bbox[3][1] - qrdataj->bbox[1][1]);
            horiz = abs(dir[0]) > abs(dir[1]);
            (*sym)->orient = (zbar_orientation_t)(horiz + 2 * (dir[1 - horiz] < 0));

            for (k = 0; k < qrdataj->nentries && !err; k++) {
                size_t              inleft;
                size_t              outleft;
                char* in;
                char* out;
                entry = qrdataj->entries + k;
                switch (entry->mode) {
                case QR_MODE_NUM: {
                    if (sa_ctext - sa_ntext >= (size_t)entry->payload.data.len) {
                        memcpy(sa_text + sa_ntext, entry->payload.data.buf,
                            entry->payload.data.len * sizeof(*sa_text));
                        sa_ntext += entry->payload.data.len;
                    }
                    else err = 1;
                }break;
                case QR_MODE_ALNUM: {
                    char* p;
                    in = (char*)entry->payload.data.buf;
                    inleft = entry->payload.data.len;
                    if (fnc1)for (;;) {
                        size_t plen;
                        char   c;
                        p = (char*)memchr(in, '%', inleft * sizeof(*in));
                        if (p == nullptr)break;
                        plen = p - in;
                        if (sa_ctext - sa_ntext < plen + 1)break;
                        memcpy(sa_text + sa_ntext, in, plen * sizeof(*in));
                        sa_ntext += plen;
                        if (plen + 1 < inleft && p[1] == '%') {
                            c = '%';
                            plen++;
                            p++;
                        }
                        else c = 0x1D;
                        sa_text[sa_ntext++] = c;
                        inleft -= plen + 1;
                        in = p + 1;
                    }
                    else p = nullptr;
                    if (p != nullptr || sa_ctext - sa_ntext < inleft)err = 1;
                    else {
                        memcpy(sa_text + sa_ntext, in, inleft * sizeof(*sa_text));
                        sa_ntext += inleft;
                    }
                }break;
                case QR_MODE_BYTE:
                case QR_MODE_KANJI: {
                    in = (char*)entry->payload.data.buf;
                    inleft = entry->payload.data.len;
                    out = sa_text + sa_ntext;
                    outleft = sa_ctext - sa_ntext;
                    if (eci < 0) {
                        int ei;
                        if (has_kanji)enc_list_mtf(enc_list, sjis_cd);
                        else if (inleft >= 3 &&
                            in[0] == (char)0xEF && in[1] == (char)0xBB && in[2] == (char)0xBF) {
                            in += 3;
                            inleft -= 3;
                            err = utf8_cd == (iconv_t)-1 ||
                                iconv(utf8_cd, &in, &inleft, &out, &outleft) == (size_t)-1;
                            if (!err) {
                                sa_ntext = out - sa_text;
                                enc_list_mtf(enc_list, utf8_cd);
                                continue;
                            }
                            in = (char*)entry->payload.data.buf;
                            inleft = entry->payload.data.len;
                            out = sa_text + sa_ntext;
                            outleft = sa_ctext - sa_ntext;
                        }
                        else if (text_is_ascii((unsigned char*)in, inleft)) {
                            enc_list_mtf(enc_list, utf8_cd);
                        }
                        for (ei = 0; ei < 3; ei++)if (enc_list[ei] != (iconv_t)-1) {
                            if (ei < 2 && enc_list[ei] == latin1_cd &&
                                !text_is_latin1((unsigned char*)in, inleft)) {
                                int ej;
                                for (ej = ei + 1; ej < 3; ej++)enc_list[ej - 1] = enc_list[ej];
                                enc_list[2] = latin1_cd;
                            }
                            err = iconv(enc_list[ei], &in, &inleft, &out, &outleft) == (size_t)-1;
                            if (!err) {
                                sa_ntext = out - sa_text;
                                enc_list_mtf(enc_list, enc_list[ei]);
                                break;
                            }
                            in = (char*)entry->payload.data.buf;
                            inleft = entry->payload.data.len;
                            out = sa_text + sa_ntext;
                            outleft = sa_ctext - sa_ntext;
                        }
                    }
                    else {
                        err = eci_cd == (iconv_t)-1 ||
                            iconv(eci_cd, &in, &inleft, &out, &outleft) == (size_t)-1;
                        if (!err)sa_ntext = out - sa_text;
                    }
                }break;
                case QR_MODE_ECI: {
                    const char* enc;
                    char        buf[16];
                    unsigned    cur_eci;
                    cur_eci = entry->payload.eci;
                    if (cur_eci <= QR_ECI_ISO8859_16 && cur_eci != 14) {
                        if (cur_eci != QR_ECI_GLI0 && cur_eci != QR_ECI_CP437) {
                            enc = buf;
                        }
                        else enc = "CP437";
                    }
                    else if (cur_eci == QR_ECI_SJIS)enc = "SJIS";
                    else if (cur_eci == QR_ECI_UTF8)enc = "UTF-8";
                    else continue;
                    eci = cur_eci;
                    eci_cd = iconv_open("UTF-8", enc);
                }break;
                default:break;
                }
            }
            if (eci <= QR_ECI_GLI1) {
                eci = -1;
                if (eci_cd != (iconv_t)-1)iconv_close(eci_cd);
            }
        }
        if (eci_cd != (iconv_t)-1)iconv_close(eci_cd);
        if (!err) {
            zbar_symbol_t* sa_sym;
            sa_text[sa_ntext++] = '\0';
            if (sa_ctext + 1 > sa_ntext) {
                sa_text = (char*)realloc(sa_text, sa_ntext * sizeof(*sa_text));
            }

            if (sa_size == 1)
                sa_sym = syms;
            else {
                int xmin = img->width, xmax = -2;
                int ymin = img->height, ymax = -2;

                sa_sym = _zbar_image_scanner_alloc_sym(iscn, ZBAR_QRCODE, 0);
                sa_sym->syms = _zbar_symbol_set_create();
                sa_sym->syms->head = syms;

                for (; syms; syms = syms->next) {
                    int next;
                    _zbar_symbol_refcnt(syms, 1);
                    if (syms->type == ZBAR_PARTIAL)
                        sa_sym->type = ZBAR_PARTIAL;
                    else
                        for (j = 0; j < syms->npts; j++) {
                            int u = syms->pts[j].x;
                            if (xmin >= u) xmin = u - 1;
                            if (xmax <= u) xmax = u + 1;
                            u = syms->pts[j].y;
                            if (ymin >= u) ymin = u - 1;
                            if (ymax <= u) ymax = u + 1;
                        }
                    syms->data = sa_text + syms->datalen;
                    next = (syms->next) ? syms->next->datalen : sa_ntext;
                    assert(next > syms->datalen);
                    syms->datalen = next - syms->datalen - 1;
                }
                if (xmax >= -1) {
                    sym_add_point(sa_sym, xmin, ymin);
                    sym_add_point(sa_sym, xmin, ymax);
                    sym_add_point(sa_sym, xmax, ymax);
                    sym_add_point(sa_sym, xmax, ymin);
                }
            }
            sa_sym->data = sa_text;
            sa_sym->data_alloc = sa_ntext;
            sa_sym->datalen = sa_ntext - 1;
            sa_sym->modifiers = fnc1;

            _zbar_image_scanner_add_sym(iscn, sa_sym);
        }
        else {
            _zbar_image_scanner_recycle_syms(iscn, syms);
            free(sa_text);
        }
    }
    if (utf8_cd != (iconv_t)-1)iconv_close(utf8_cd);
    if (sjis_cd != (iconv_t)-1)iconv_close(sjis_cd);
    if (latin1_cd != (iconv_t)-1)iconv_close(latin1_cd);
    free(mark);
    return ntext;
}

void qr_code_data_list_clear(qr_code_data_list* _qrlist) {
    int i;
    for (i = 0; i < _qrlist->nqrdata; i++)qr_code_data_clear(_qrlist->qrdata + i);
    free(_qrlist->qrdata);
    qr_code_data_list_init(_qrlist);
}

int minizbar::_zbar_qr_decode(qr_reader* reader, zbar_image_scanner_t* iscn, zbar_image_t* img)
{
    int nqrdata = 0, ncenters;
    qr_finder_edge_pt* edge_pts = nullptr;
    qr_finder_center* centers = nullptr;

    if (reader->finder_lines[0].nlines < 9 ||
        reader->finder_lines[1].nlines < 9)
        return(0);
    ncenters = qr_finder_centers_locate(&centers, &edge_pts, reader, 0, 0);

    if (ncenters >= 3) {
        void* bin = qr_binarize((const unsigned char*)img->data, img->width, img->height);

        qr_code_data_list qrlist;
        qr_code_data_list_init(&qrlist);

        qr_reader_match_centers(reader, &qrlist, centers, ncenters,
            (const unsigned char*)bin, img->width, img->height);

        if (qrlist.nqrdata > 0)
            nqrdata = qr_code_data_list_extract_text(&qrlist, iscn, img);

        qr_code_data_list_clear(&qrlist);
        free(bin);
    }

    if (centers)
        free(centers);
    if (edge_pts)
        free(edge_pts);
    return(nqrdata);
}

int minizbar::_zbar_get_symbol_hash(zbar_symbol_type_t sym)
{
    static const signed char hash[0x20] = {
    0x00, 0x01, 0x10, 0x11,   -1, 0x11, 0x16, 0x0c,
0x05, 0x06, 0x08,   -1, 0x04, 0x03, 0x07, 0x12,
  -1,   -1,   -1,   -1,   -1,   -1,   -1, 0x02,
  -1, 0x00, 0x12, 0x0c, 0x0b, 0x1d, 0x0a, 0x00,
    };
    int g0 = hash[sym & 0x1f];
    int g1 = hash[~(sym >> 4) & 0x1f];
    assert(g0 >= 0 && g1 >= 0);
    if (g0 < 0 || g1 < 0)
        return(0);
    return((g0 + g1) & 0x1f);
}

void minizbar::cache_sym(zbar_image_scanner_t* iscn, zbar_symbol_t* sym)
{
    if (iscn->enable_cache) {
        uint32_t age, near_thresh, far_thresh, dup;
        zbar_symbol_t* entry = cache_lookup(iscn, sym);
        if (!entry) {
            entry = _zbar_image_scanner_alloc_sym(iscn, sym->type,
                sym->datalen + 1);
            entry->configs = sym->configs;
            entry->modifiers = sym->modifiers;
            memcpy(entry->data, sym->data, sym->datalen);
            entry->time = sym->time - CACHE_HYSTERESIS;
            entry->cache_count = 0;
            /* add to cache */
            entry->next = iscn->cache;
            iscn->cache = entry;
        }
        age = sym->time - entry->time;
        entry->time = sym->time;
        near_thresh = (age < CACHE_PROXIMITY);
        far_thresh = (age >= CACHE_HYSTERESIS);
        dup = (entry->cache_count >= 0);
        if ((!dup && !near_thresh) || far_thresh) {
            int type = sym->type;
            int h = _zbar_get_symbol_hash((zbar_symbol_type_t)type);
            entry->cache_count = -iscn->sym_configs[0][h];
        }
        else if (dup || near_thresh)
            entry->cache_count++;

        sym->cache_count = entry->cache_count;
    }
    else
        sym->cache_count = 0;
}

void minizbar::_zbar_image_scanner_add_sym(zbar_image_scanner_t* iscn, zbar_symbol_t* sym)
{
    zbar_symbol_set_t* syms;
    cache_sym(iscn, sym);

    syms = iscn->syms;
    if (sym->cache_count || !syms->tail) {
        sym->next = syms->head;
        syms->head = sym;
    }
    else {
        sym->next = syms->tail->next;
        syms->tail->next = sym;
    }

    if (!sym->cache_count)
        syms->nsyms++;
    else if (!syms->tail)
        syms->tail = sym;

    _zbar_symbol_refcnt(sym, 1);
}

inline zbar_symbol_t* minizbar::_zbar_image_scanner_alloc_sym(zbar_image_scanner_t* iscn, zbar_symbol_type_t type, int datalen)
{
    zbar_symbol_t* sym = nullptr;
    int i;
    for (i = 0; i < RECYCLE_BUCKETS - 1; i++)
        if (datalen <= 1 << (i * 2))
            break;

    for (; i > 0; i--)
        if ((sym = iscn->recycle[i].head)) {
            STAT(sym_recycle[i]);
            break;
        }

    if (sym) {
        iscn->recycle[i].head = sym->next;
        sym->next = nullptr;
        assert(iscn->recycle[i].nsyms);
        iscn->recycle[i].nsyms--;
    }
    else {
        sym = (zbar_symbol_t*)calloc(1, sizeof(zbar_symbol_t));
        STAT(sym_new);
    }

    sym->type = type;
    sym->quality = 1;
    sym->npts = 0;
    sym->orient = ZBAR_ORIENT_UNKNOWN;
    sym->cache_count = 0;
    sym->time = iscn->time;
    assert(!sym->syms);

    if (datalen > 0) {
        sym->datalen = datalen - 1;
        if (sym->data_alloc < datalen) {
            if (sym->data)
                free(sym->data);
            sym->data_alloc = datalen;
            sym->data = (char*)malloc(datalen);
        }
    }
    else {
        if (sym->data)
            free(sym->data);
        sym->data = nullptr;
        sym->datalen = sym->data_alloc = 0;
    }
    return(sym);
}

zbar_symbol_t* minizbar::cache_lookup(zbar_image_scanner_t* iscn, zbar_symbol_t* sym)
{
    zbar_symbol_t** entry = &iscn->cache;
    while (*entry) {
        if ((*entry)->type == sym->type &&
            (*entry)->datalen == sym->datalen &&
            !memcmp((*entry)->data, sym->data, sym->datalen))
            break;
        if ((sym->time - (*entry)->time) > CACHE_TIMEOUT) {
            zbar_symbol_t* next = (*entry)->next;
            (*entry)->next = nullptr;
            _zbar_image_scanner_recycle_syms(iscn, *entry);
            *entry = next;
        }
        else
            entry = &(*entry)->next;
    }
    return(*entry);
}

static inline void minizbar::quiet_border(zbar_image_scanner_t* iscn)
{
    zbar_scanner_t* scn = iscn->scn;
    zbar_scanner_flush(scn);
    zbar_scanner_flush(scn);
    zbar_scanner_new_scan(scn);
}

static inline unsigned minizbar::calc_thresh(zbar_scanner_t* scn)
{
    unsigned dx, thresh = scn->y1_thresh;
    unsigned long t;
    if ((thresh <= scn->y1_min_thresh) || !scn->width) {
        return(scn->y1_min_thresh);
    }
    dx = (scn->x << ZBAR_FIXED) - scn->last_edge;
    t = thresh * dx;
    t /= scn->width;
    t /= ZBAR_SCANNER_THRESH_FADE;
    if (thresh > t) {
        thresh -= t;
        if (thresh > scn->y1_min_thresh)
            return(thresh);
    }
    scn->y1_thresh = scn->y1_min_thresh;
    return(scn->y1_min_thresh);
}

zbar_symbol_type_t minizbar::zbar_scan_y(zbar_scanner_t* scn, int y)
{
    register int x = scn->x;
    register int y0_1 = scn->y0[(x - 1) & 3];
    register int y0_0 = y0_1;
    register int y0_2, y0_3, y1_1, y2_1, y2_2;
    zbar_symbol_type_t edge;
    if (x) {
        y0_0 += ((int)((y - y0_1) * EWMA_WEIGHT)) >> ZBAR_FIXED;
        scn->y0[x & 3] = y0_0;
    }
    else
        y0_0 = y0_1 = scn->y0[0] = scn->y0[1] = scn->y0[2] = scn->y0[3] = y;
    y0_2 = scn->y0[(x - 2) & 3];
    y0_3 = scn->y0[(x - 3) & 3];
    y1_1 = y0_1 - y0_2;
    {
        register int y1_2 = y0_2 - y0_3;
        if ((abs(y1_1) < abs(y1_2)) &&
            ((y1_1 >= 0) == (y1_2 >= 0)))
            y1_1 = y1_2;
    }
    y2_1 = y0_0 - (y0_1 * 2) + y0_2;
    y2_2 = y0_1 - (y0_2 * 2) + y0_3;

    edge = ZBAR_NONE;
    if ((!y2_1 ||
        ((y2_1 > 0) ? y2_2 < 0 : y2_2 > 0)) &&
        (calc_thresh(scn) <= abs(y1_1)))
    {
        char y1_rev = (scn->y1_sign > 0) ? y1_1 < 0 : y1_1 > 0;
        if (y1_rev)
            edge = process_edge(scn, y1_1);

        if (y1_rev || (abs(scn->y1_sign) < abs(y1_1))) {
            int d;
            scn->y1_sign = y1_1;
            scn->y1_thresh = (abs(y1_1) * THRESH_INIT + ROUND) >> ZBAR_FIXED;
            if (scn->y1_thresh < scn->y1_min_thresh)
                scn->y1_thresh = scn->y1_min_thresh;

            d = y2_1 - y2_2;
            scn->cur_edge = 1 << ZBAR_FIXED;
            if (!d)
                scn->cur_edge >>= 1;
            else if (y2_1)
                scn->cur_edge -= ((y2_1 << ZBAR_FIXED) + 1) / d;
            scn->cur_edge += x << ZBAR_FIXED;
        }
    }
    scn->x = x + 1;
    return(edge);
}

static inline char minizbar::release_lock(zbar_decoder_t* dcode, zbar_symbol_type_t req)
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
        memset(scn->decoder->w, 0, sizeof(scn->decoder->w));
        scn->decoder->lock = (zbar_symbol_type_t)0;
        scn->decoder->idx = 0;
        scn->decoder->s6 = 0;
        scn->decoder->qrf.s5 = 0;

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
        }
        else if (!strncmp(cfgstr, "qrcode", len)) {
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
    unsigned* config = &dcode->qrf.config;
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

//inline void minizbar::zbar_image_free_data(zbar_image_t* img)
//{
//    if (!img) {
//        return;
//    }
//    if (img->cleanup && img->data) {
//        if (img->cleanup != zbar_image_free_data) {
//            zbar_image_cleanup_handler_t* cleanup = img->cleanup;
//            img->cleanup = zbar_image_free_data;
//            cleanup(img);
//        }
//    }
//    else {
//        free((void*)img->data);
//    }
//    img->data = nullptr;
//}

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
    line->boffs = u - zbar_scanner_get_edge(iscn->scn, line->boffs, QR_FINDER_SUBPREC);
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
    return (scn);
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





