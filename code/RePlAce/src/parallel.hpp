#include "global.h"
#include "bin.h"
#include "charge.h"
#include "wlen.h"
#include "fft.h"
#include <omp.h>
#include <cmath>
#include <bits/stdc++.h>

void __wirelength__(circuit_t* circuit, FPOS *st, fpos2_t** poff, size_t* wrk_ld, double* time);//, cell_t* cells);
void __FILLgradient__(FPOS *dst, FPOS *wdst,
                          FPOS *pdst, FPOS *pdstl, int N,
                          prec *cellLambdaArr, bool onlyPreCon, bin_t* bins, cell_phy_t* ios, net_t* nets, double* time, cell_den_t* cels);
void __gradient__(FPOS *dst, FPOS *wdst,
                          FPOS *pdst, FPOS *pdstl, int N,
                          prec *cellLambdaArr, bool onlyPreCon, bin_t* bins, cell_phy_t* ios, net_t* nets, double* time, size_t* W1, cell_den_t* cels);
void __bin_update7_cGP2D(cell_den_t* cells, bin_t* bins, area_t* areas, float** localAr, float** localAr2, size_t* bound, double* time, double* time2);
void update_cells_frame(cell_den_t* cells) ;