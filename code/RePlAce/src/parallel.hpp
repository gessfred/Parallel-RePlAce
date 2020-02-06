#include "global.h"
#include "bin.h"
#include "charge.h"
#include "wlen.h"
#include "fft.h"
#include <omp.h>
#include <cmath>
#include <bits/stdc++.h>

/**
update_wirelength computes the weighted average wirelength and caches some subterms/subfactors
@param circuit - instance must contain
@param st - the current 
@param poff - pin offsets
@param wrkl_ld - number of nets that each thread must take so that each iterates on an equal number of pins
@param time - time to execute the function
*/
void update_wirelength(circuit_t* circuit, FPOS *st, fpos2_t** poff, size_t* wrk_ld, double* time);
/**
update_fill_gradient 
*/
void update_fill_gradient(FPOS *dst, FPOS *wdst,
                          FPOS *pdst, FPOS *pdstl, int N,
                          prec *cellLambdaArr, bool onlyPreCon, bin_t* bins, cell_phy_t* ios, net_t* nets, double* time, cell_den_t* cels);
void update_gradient(FPOS *dst, FPOS *wdst,
                          FPOS *pdst, FPOS *pdstl, int N,
                          prec *cellLambdaArr, bool onlyPreCon, bin_t* bins, cell_phy_t* ios, net_t* nets, double* time, size_t* W1, cell_den_t* cels);
/**
update_field_potential updates the field and potential
*/
void update_field_potential(cell_den_t* cells, bin_t* bins, area_t* areas, float** localAr, float** localAr2, size_t* bound, double* time, double* time2);
/**
update_cells_frame updates each cell's frame according to 
*/
void update_cells_frame(cell_den_t* cells);