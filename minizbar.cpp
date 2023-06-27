#include "minizbar.h"
#include "utils.h"
#include <chrono>
using namespace minizbar;

zbar_image_scanner_t* minizbar::zbar_image_scanner_create()
{
    zbar_image_scanner_t* iscn = (zbar_image_scanner_t*)calloc(1, sizeof(zbar_image_scanner_t));
    if (!iscn)
        return(nullptr);
    iscn->dcode = zbar_decoder_create();
    iscn->scn = zbar_scanner_create(iscn->dcode);
    if (!iscn->dcode || !iscn->scn) {
        zbar_image_scanner_destroy(iscn);
        return(NULL);
    }
    iscn->dcode->userdata = iscn;
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

static inline unsigned long minizbar::zbar_fourcc_parse(const char* format)
{
    unsigned long fourcc = 0;
    if (format) {
        for (int i = 0; i < 4 && format[i]; ++i) {
            fourcc |= ((unsigned long)format[i]) << (i * 8);
        }
    }
    return (fourcc);
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
    if (sym != ZBAR_NONE || sym != ZBAR_QRCODE) {
        throw "No Support Other Type";
    }
    if ( sym == 0  && cfg == ZBAR_CFG_ENABLE) {
        iscn->ean_config = !!val;
        if (sym)
            return(0);
    }

    if (cfg < ZBAR_CFG_UNCERTAINTY)
        return(zbar_decoder_set_config(iscn->dcode, sym, cfg, val));
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

                    /* recycle */
                    *symp = sym->next;
                    syms->nsyms--;
                    sym->next = NULL;
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
            /* create container symbol for composite result */
            zbar_symbol_t* ean = NULL, * addon = NULL;
            for (symp = &syms->head; *symp; ) {
                zbar_symbol_t* sym = *symp;
                if (sym->type < ZBAR_COMPOSITE && sym->type > ZBAR_PARTIAL) {
                    /* move to composite */
                    *symp = sym->next;
                    syms->nsyms--;
                    sym->next = NULL;
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


