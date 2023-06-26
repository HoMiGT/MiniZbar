#include "minizbar.h"

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
