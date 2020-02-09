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
*/
void update_wirelength(circuit_t* circuit, FPOS *st);
/**
update_fill_gradient 
*/
void update_fill_gradient(FPOS *dst, FPOS *wdst,
                          FPOS *pdst, FPOS *pdstl, int N,
                          prec *cellLambdaArr, bool onlyPreCon, circuit_t* circuit);
void update_gradient(FPOS *dst, FPOS *wdst,
                          FPOS *pdst, FPOS *pdstl, int N,
                          prec *cellLambdaArr, bool onlyPreCon, circuit_t* circuit);
/**
update_field_potential updates the field and potential
*/
void update_field_potential(circuit_t* circuit);
/**
update_cells_frame updates each cell's frame according to 
*/
void update_cells_frame(cell_den_t* cells);