///////////////////////////////////////////////////////////////////////////////
// Authors: Ilgweon Kang and Lutong Wang
//          (respective Ph.D. advisors: Chung-Kuan Cheng, Andrew B. Kahng),
//          based on Dr. Jingwei Lu with ePlace and ePlace-MS
//
//          Many subsequent improvements were made by Mingyu Woo
//          leading up to the initial release.
//
// BSD 3-Clause License
//
// Copyright (c) 2018, The Regents of the University of California
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
///////////////////////////////////////////////////////////////////////////////

#ifndef __PL_BIN__
#define __PL_BIN__

#include "global.h"

extern int tot_bin_cnt;
extern int bin_cnt;
extern prec bin_area;
extern prec inv_bin_area;
extern prec global_macro_area_scale;

extern BIN *bin_mat;
//extern prec *bin_share_x_st;
//extern prec *bin_share_y_st;
//extern prec *bin_share_z_st;
//extern prec *bin_share_st;


extern FPOS bin_org;
extern FPOS bin_stp;
extern FPOS inv_bin_stp;
extern FPOS half_bin_stp;
extern FPOS bin_stp_mGP2D;
extern FPOS bin_stp_cGP2D;
extern POS max_bin;

//extern BIN **bin_list;
//extern int bin_list_cnt;

//extern POS *bin_st;

struct BIN {
    FPOS e;
    ////igkang
    //  FPOS    e_local;
    POS p;
    FPOS pmin;
    FPOS pmax;
    FPOS center;
    prec cell_area;
    prec cell_area2;

    prec virt_area;
    long term_area; // mgwoo

    prec filler_area;
    prec phi;
    int flg;
    prec pl_area;
    prec macro_area;
    prec macro_area2;
    prec den;
    prec den2;
    prec no_mac_den;
    void dump(string a) {
        cout << a << endl;
        e.Dump("e");
        p.Dump("p");
        pmin.Dump("pmin");
        pmax.Dump("pmax");
        center.Dump("center");
        cout << endl;
    }
};
struct bin_t {
    fpos2_t pmax;
    fpos2_t pmin;
    fpos2_t e;
    prec cell_area;
    prec cell_area2;
    prec phi;
    inline void from(BIN* bin) {
        pmax.from(bin->pmax);
        pmin.from(bin->pmin);
        e.from(bin->e);
        cell_area = bin->cell_area;
        cell_area2 = bin->cell_area2;
        phi = bin->phi;
    }
    inline void to(BIN* bin) {
        bin->pmax.from(pmax);
        bin->pmin.from(pmin);
        bin->e.from(e);
        bin->cell_area = cell_area;
        bin->cell_area2 = cell_area2;
        bin->phi = phi;
    }
};

struct Bin_t {
    float* pmax_x;
    float* pmax_y;
    float* pmin_x;
    float* pmin_y;
    float* e_x;
    float* e_y;
    float* cell_area;
    float* cell_area2;
    float* phi;
    size_t size;

    inline void build(size_t N) {
        pmax_x = (float*)calloc(sizeof(float), N);
        pmax_y = (float*)calloc(sizeof(float), N);
        pmin_x = (float*)calloc(sizeof(float), N);
        pmin_y = (float*)calloc(sizeof(float), N);
        e_x = (float*)calloc(sizeof(float), N);
        e_y = (float*)calloc(sizeof(float), N);
        cell_area = (float*)calloc(sizeof(float), N);
        cell_area2 = (float*)calloc(sizeof(float), N);
        phi = (float*)calloc(sizeof(float), N);
        size = N;
    }

    inline void copy(bin_t* bins) {
        for(int i = 0; i < size; ++i) {
            pmax_x[i] = bins[i].pmax.x;
            pmax_y[i] = bins[i].pmax.y;
            pmin_x[i] = bins[i].pmin.x;
            pmin_y[i] = bins[i].pmin.y;
            e_x[i] = bins[i].e.x;
            e_y[i] = bins[i].e.y;
            cell_area[i] = bins[i].cell_area;
            cell_area2[i] = bins[i].cell_area2;
            phi[i] = bins[i].phi;
        }
    }

    inline void copyback(bin_t* bins) {
        for(int i = 0; i < size; ++i) {
            bins[i].pmax.x = pmax_x[i];
            bins[i].pmax.y = pmax_y[i];
            bins[i].pmin.x = pmin_x[i];
            bins[i].pmin.y = pmin_y[i];
            bins[i].e.x = e_x[i];
            bins[i].e.y = e_y[i];
            bins[i].cell_area = cell_area[i];
            bins[i].cell_area2 = cell_area2[i];
            bins[i].phi = phi[i];
        }
    }

    inline void destroy() {
        free(pmax_x);
        free(pmax_y);
        free(pmin_x);
        free(pmin_y);
        free(e_x);
        free(e_y);
        free(cell_area);
        free(cell_area2);
        free(phi);
        size = 0;
    }
};


//bins as a collection of area data
struct area_t {
    prec virt_area;
    long term_area;
    prec den;
    prec den2;
    pos_t p;
    inline void from(BIN* bin) {
        virt_area = bin->virt_area;
        term_area = bin->term_area;
        den = bin->den;
        den2 = bin->den2;
        p.set(bin->p.x, bin->p.y);
    }
    inline void to(BIN* bin) {
        bin->den = den;
        bin->den2 =den2;
        bin->p.from(p);
    }
};



int idx_in_bin_rect(POS *p, POS pmin, POS pmax);

void bin_init();
void bin_init_2D(int);

void UpdateTerminalArea( TIER* tier, FPOS* pmin, FPOS* pmax); 
void bin_update(int N);
void bin_update7(int N);

void bin_update7_cGP2D(double* time0, double* time1);
void Density(Cell_t* cells, Bin_t* bins, int* bound, float** localAr, float** localAr2);
void __density__(Cell_t* cells, Bin_t* bins, TIER *tier);
void __bin_update7_cGP2D(cell_den_t* cells, bin_t* bins, area_t* areas, float** localAr, float** localAr2, size_t* bound, double* time, double* time2);
void CellBox(Cell_t* cells);
void FFTSolve(bin_t* bins, area_t* areas);
void bin_update7_mGP2D();

void bin_delete(void);

void get_bin_grad(BIN **bin, int max_x, int max_y);

//int wthin(prec a, prec min_a, prec max_a);
//prec get_int_line(prec min_x1, prec max_x1, prec min_x2, prec max_x2);
//FPOS get_int_rgn(FPOS min1, FPOS max1, FPOS min2, FPOS max2);

void legal_bin_idx(POS *p);
POS get_bin_pt_from_point(FPOS p);

// calculating the da variable.
// using current da & half_size.
//
inline prec valid_coor2(prec da, prec half_size, int lab) {
    prec min_a = da - half_size;
    prec max_a = da + half_size;

    // according to X
    if(lab == 0) {
        if(min_a < place.org.x)
            da += place.org.x - min_a;
        if(max_a > place.end.x)
            da -= max_a - place.end.x;
    }
    // accroding to Y
    else if(lab == 1) {
        if(min_a < place.org.y)
            da += place.org.y - min_a;
        if(max_a > place.end.y)
            da -= max_a - place.end.y;
    }
    // mgwoo
    // accroding to Z
//    else if(lab == 2) {
//        if(min_a < place.org.z)
//            da += place.org.z - min_a;
//        if(max_a > place.end.z)
//            da -= max_a - place.end.z;
//    }
    return da;
}

inline FPOS valid_coor00(FPOS v, FPOS half_size) {
//    FPOS v1 = zeroFPoint;
    FPOS v1;
    v1.x = valid_coor2(v.x, half_size.x, 0);
    v1.y = valid_coor2(v.y, half_size.y, 1);
    return v1;
}

prec valid_coor(prec a, int xy);
prec valid_coor3(prec da, prec sz, int lab);

int is_IO_block(TERM *term);

void get_bins(FPOS center, CELLx *cell, POS *st, prec *share_st, int *bin_cnt);
prec get_bins_mac(FPOS center, MODULE *mac);
void fft_test(void);

#define DEN_SMOOTH_COF 5.0
enum { SIN_SMOOTH, LIN_SMOOTH };
#define SMOOTH_LAB LIN_SMOOTH /* SIN_SMOOTH   */

int find_next_dbin(int *x0, int *y0);

void loc_ext(int x0, int y0, int dir, int flg);
int find_next_loc_opt(POS *p);

int point_in_range(prec a, prec min_a, prec max_a);

FPOS valid_coor4(FPOS center, FPOS half_sz);

//#define GRAD_SCAL
//#define GRAD_SCAL_COMB

prec get_ovfl();

void bin_den_update(FPOS *x_st, int N, int iter, int lab);
bool get_common_rect(FPOS pmin1, FPOS pmax1, FPOS pmin2, FPOS pmax2, FPOS *p1,
                     FPOS *p2);
prec get_common_area(FPOS pmin1, FPOS pmax1, FPOS pmin2, FPOS pmax2);
prec get_common_vol(FPOS pmin1, FPOS pmax1, FPOS pmin2, FPOS pmax2);

prec get_den(prec area_num, prec area_denom);
prec get_den2(prec area_num, prec area_denom);

// return bin_mat's pointer
inline BIN *get_bin_from_idx(POS p) {
    int idx = 0;
    (flg_3dic == 1) ? idx = p.x *max_bin.y *max_bin.z + p.y *max_bin.z + p.z
                    : idx = p.x * max_bin.y + p.y;
    return &bin_mat[idx];
};

void add_net_to_bins(NET *ntp);
int point_in_rect(FPOS p0, FPOS pmin, FPOS pmax);

////// MACRO LG //////
void den_update();
void bin_den_mac_update(void);
void bin_cell_area_update(void);
void bin_clr(void);
prec get_mac_den(int idx);

POS get_bin_idx(int idx);

prec get_all_macro_den(void);

void bin_den_update_for_filler_adj(void);
//////-MACRO LG-//////

enum { NoneDpre, AreaDpre, DenDpre };
enum { NoneTemppre, AreaTemppre, DenTemppre };
#define CHARGE_PRE /* NoneDpre */ AreaDpre     /* DenDpre */
#define TEMP_PRE /* NoneTemppre */ AreaTemppre /* DenTemppre */

void den_comp(int cell_idx);
inline void den_comp_2d_mGP2D(CELLx *cell, TIER *tier);
inline void den_comp_2d_cGP2D(CELLx *cell, TIER *tier);
void den_comp_3d(int cell_idx);

// void    bin_zum_z ();
void bin_delete_mGP2D(void);
void bin_init_mGP2D(void);
void bin_init_cGP2D(void);

int compare(const void *p1, const void *p2);
void __FFTsolve__(cell_den_t* cell, bin_t* bins, area_t* areas);

#endif
