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


#ifndef __GLOBAL__
#define __GLOBAL__

#include <stdint.h>
#include <algorithm>
#include <climits>
#include <iostream>
#include <map>
#include <string>
#include <tuple>
#include <vector>
#include <iomanip>
#include <google/dense_hash_map>
#include <cfloat>
#include <ctime>
#include <cassert>
#include <atomic>
#include <chrono>

#include <bitset>
#include "data_structures.h"
#include "parallel_data_structures.hpp"
#define PI 3.141592653589793238462L
#define SQRT2 1.414213562373095048801L
#define INV_SQRT2 0.707106781186547524401L


// for PREC_MODE variable => required for different codes. 
#define IS_FLOAT 0
#define IS_DOUBLE 1

// precision settings
#define PREC_MODE IS_FLOAT

#if PREC_MODE == IS_FLOAT

typedef float prec;
#define PREC_MAX FLT_MAX 
#define PREC_MIN FLT_MIN

#elif PREC_MODE == IS_DOUBLE

typedef double prec;
#define PREC_MAX DBL_MAX 
#define PREC_MIN DBL_MIN

#endif

#define INT_CONVERT(a)          (int)(1.0*(a) + 0.5f)
#define INT_DOWN(a)             (int)(a)
#define UNSIGNED_CONVERT(a)     (unsigned)(1.0*(a) + 0.5f)

#define TSV_CAP 30.0  // fF
// data from M. Jung et al. "How to Reduce Power in 3D IC Designs:
// A Case Study with OpenSPARC T2 Core", CICC 2013, Section 2.B
#define MIV_CAP 0.1  // fF
// data from S. Panth et al. "Power-Performance Study of Block-Level
// Monolithc 3D-ICs Considering Inter-Tier Performance Variations",
// DAC 2014, Section 4.
#define ROW_CAP 0.3  // fF
// row_hei is 1.4um at 45nm tech-node with 70nm M1 half-pitch and
// 10 M1 tracks per placement row, unit capacitance per centimeter
// at intermediate routing layers are constantly 2pF/cm from ITRS
// INTC2 tabls between 2008 and 2013, which makes ROW_CAP roughly
// to be 0.3fF/row

#define BUF_SZ 511
#define Epsilon 1.0E-15
#define MIN_AREA3 /* 1.0 */ 1.0e-12
#define MIN_LEN 25.0 /* 10.0 */ /* 5.0 */ /* 1.0 */
#define LS_DEN
#define DetailPlace
#define SA_LG
#define FILLER_ADJ RandomAdj
#define TIER_DEP /* 64.0 */ /* 600.0 */ 1
#define TIER_Z0 0
#define MSH_Z_RES /* 8 */ 1
#define THETA_XY_3D_PLOT PI / 6.0
#define Z_SCAL 1.00

// GP_Dimension_one, 
// 1: true, 0: false
#define GP_DIM_ONE_FOR_IP 1
#define GP_DIM_ONE_FOR_GP1 1
#define GP_DIM_ONE_FOR_GP2 1
#define GP_DIM_ONE_FOR_MAC 1

#define GP_DIM_ONE_FOR_mGP2D 0
#define GP_DIM_ONE_FOR_GP3 0

#define GP_SCAL /* 1000.0 */ /* 1.0 */ /* 0.001 */ 0.0002
#define DEN_GRAD_SCALE 1.0 /* 0.1 */ /* 0.5 */ /* 0.25 */ /* 0.125 */
#define LAYER_ASSIGN_3DIC MIN_TIER_ORDER                  /* MAX_AREA_DIS_DIV */

#define time_since(start)  (std::chrono::duration<double>(std::chrono::steady_clock::now()-start)).count(); start = std::chrono::steady_clock::now()
///////////////////////////////////////////////////////////////////////////

using std::string;
using std::cout;
using std::endl;

using std::vector;
using std::pair;
using std::tuple;
using std::max;
using std::min;
using std::bitset;
using std::fixed;
using std::setprecision;
using std::map;
using std::to_string;

using google::dense_hash_map;
struct POS;
struct FPOS;
struct NET;
struct PIN;
struct CELLx;

struct pin_t;
struct pos2_t;
//struct net_t;
struct fpos2_t;
struct net_t;

inline void FPOS::from(fpos2_t fp) {
        x = fp.x;
        y = fp.y;
}

inline void POS::from(pos2_t p ) {
    x = p.x;
    y = p.y;
}

inline void FPOS::Set(POS p) {
    x = p.x;
    y = p.y;
    z = p.z;
}
       
inline void POS::Set(FPOS fp) { 
    x = (int)(fp.x + 0.5f);
    y = (int)(fp.y + 0.5f);
    z = (int)(fp.z + 0.5f);
}

inline void pin_t::from(PIN* pin) {
        expL.from(pin->e1);
        expR.from(pin->e2);
        coord.from(pin->fp);
        pinID = pin->pinIDinModule;
        __getMeta__();
        netID = pin->netID;
        moduleID = pin->moduleID;
    }
inline void pin_t::to(PIN* pin) {
    pin->e1.from(expL);
    pin->e2.from(expR);
    pin->fp.from(coord);
    pin->term = meta[0]; 
    pin->X_MIN = meta[1]; 
    pin->X_MAX = meta[2]; 
    pin->Y_MIN = meta[3]; 
    pin->Y_MAX = meta[4]; 
    pin->flg1.x = meta[5];
    pin->flg1.y = meta[6];
    pin->flg2.x = meta[7];
    pin->flg2.y = meta[8];
}

// for saving pinName
// If it's lied in the PIN structure, it'll enlarge the runtime
extern vector< vector<string> > mPinName;
extern vector< vector<string> > tPinName;

typedef struct CELLx CELL;

inline void fpos2_t::from(FPOS fp) {
        x = fp.x;
        y = fp.y;
    }

inline void cell_den_t::from(CELLx* cell) {
    min.from(cell->den_pmin);
    max.from(cell->den_pmax);
    //pmin.from(cell->pmin);
    //pmax.from(cell->pmax);
    size.from(cell->half_den_size);
    scale = cell->den_scal;
    type = cell->flg;
}
inline void cell_den_t::to(CELLx* cell) {
    cell->den_pmin.from(min);
    cell->den_pmax.from(max);
    cell->half_den_size.from(size);
    cell->flg = type;
}

class SHAPE {
   public:
//    char *name;
//    char *prefix;

    string name, instName;

    int idx;
    prec llx;
    prec lly;
    prec width;
    prec height;

    SHAPE(string _name, string _instName, int _idx, prec _llx, prec _lly,
          prec _width, prec _height)
        : name(_name),
          instName(_instName),
          idx(_idx),
          llx(_llx),
          lly(_lly),
          width(_width),
          height(_height){};

    void Dump() {
        printf("shape[%d]: name: %s, instName: %s, llx: %lf, lly: %lf, width: %lf, height: %lf\n",
               idx, name.c_str(), instName.c_str(), llx, lly, width, height);
        fflush(stdout); 
    }
};

struct T0 {
    int z;
    prec dis;
};

// for *.scl files
// ROW -> PLACE structure
struct PLACE {
    FPOS org;
    prec area;
    FPOS center;
    FPOS stp;
    FPOS cnt;
    FPOS end;
    void Dump(string a) {
        cout << a << endl;
        org.Dump("origin");
        cout << fixed <<setprecision(0)<<"area: " << area << endl;
        center.Dump("center");
        stp.Dump("stp");
        cnt.Dump("cnt");
        end.Dump("end");
        cout << endl;
    }
};

// for multi-rows
struct TIER {
    struct FPOS bin_org;
    struct FPOS inv_bin_stp;
    struct POS dim_bin;
    struct FPOS pmin;
    struct FPOS pmax;
    struct BIN *bin_mat;
    struct FPOS center;
    struct FPOS size;
    struct ROW *row_st;
    struct MODULE **modu_st;
    struct TERM **term_st;
    struct MODULE **mac_st;
    struct CELLx **cell_st;
    struct FPOS bin_stp;
    prec area;
    prec modu_area;
    prec term_area;
    prec virt_area;
    prec filler_area;
    prec pl_area;
    prec ws_area;
    prec bin_area;
    prec tot_bin_area;
    prec inv_bin_area;
    prec sum_ovf;
    int row_cnt;
    int row_term_cnt;
    int modu_cnt;
    int filler_cnt;
    int cell_cnt;
    int mac_cnt;
    int term_cnt;
    int tot_bin_cnt;
    prec temp_mac_area;
    struct FPOS bin_off;
    struct FPOS half_bin_stp;
    struct MODULE *max_mac;
    struct CELLx **cell_st_tmp;

    // routability
    struct FPOS tile_stp;
    prec tile_area;
    struct POS dim_tile;
    int tot_tile_cnt;
    struct FPOS half_tile_stp;
    struct FPOS inv_tile_stp;
    struct FPOS tile_org;
    struct FPOS tile_off;
    struct TILE *tile_mat;
};


inline void net_t::from(NET* net) {
    sumNumL.from(net->sum_num1);
    sumDenomL.from(net->sum_denom1);
    sumNumR.from(net->sum_num2);
    sumDenomR.from(net->sum_denom2);
    //min.set(net->min_x, net->min_y);
    //max.set(net->max_x, net->max_y);
    min.set(net->terminalMin.x, net->terminalMin.y);
    max.set(net->terminalMax.x, net->terminalMax.y);
    pinCNT = net->pinCNTinObject;
    pinArray = (pin_t*)calloc(sizeof(pin_t), pinCNT);
    for(int k = 0; k < pinCNT; ++k) {
        pinArray[k].from(net->pin[k]); //this copy is unnecessary, those values are not read where this type is used
    }
}
inline void net_t::copy(NET* net) {
    sumNumL.from(net->sum_num1);
    sumDenomL.from(net->sum_denom1);
    sumNumR.from(net->sum_num2);
    sumDenomR.from(net->sum_denom2);
    //min.set(net->min_x, net->min_y);
    //max.set(net->max_x, net->max_y);
    min.set(net->min_x, net->min_y);
    max.set(net->max_x, net->max_y);
    pinCNT = net->pinCNTinObject;
    pinArray = (pin_t*)calloc(sizeof(pin_t), pinCNT);
    for(int k = 0; k < pinCNT; ++k) {
        pinArray[k].from(net->pin[k]);
    }
}



inline void cell_phy_t::from(CELLx* cell, net_t* nets) {
    type = cell->flg;
    pinCNT = cell->pinCNTinObject;
    pinArrayPtr = (pin_t**)calloc(sizeof(pin_t*), pinCNT);
    for(int k = 0; k < pinCNT; ++k) {
        PIN* p = cell->pin[k];
        pinArrayPtr[k] = &nets[p->netID].pinArray[p->pinIDinNet];
    }
}

inline void bin_t::from(BIN* bin) {
    max.from(bin->pmax);
    min.from(bin->pmin);
    field.from(bin->e);
    cellArea = bin->cell_area;
    fillerArea = bin->cell_area2;
    potential = bin->phi;
}
inline void bin_t::to(BIN* bin) {
    bin->pmax.from(max);
    bin->pmin.from(min);
    bin->e.from(field);
    bin->cell_area = cellArea;
    bin->cell_area2 = fillerArea;
    bin->phi = potential;
}
inline void area_t::from(BIN* bin) {
    virtArea = bin->virt_area;
    terminArea = bin->term_area;
    binDensity = bin->den;
    fillerDensity = bin->den2;
    coord.set(bin->p.x, bin->p.y);
}
inline void area_t::to(BIN* bin) {
    bin->den = binDensity;
    bin->den2 = fillerDensity;
    bin->p.from(coord);
}

/*inline void NET::from(net_t net) {
        min_x = net.min.x;
        min_y = net.min.y;
        max_x = net.max.x;
        max_y = net.max.y;
        for(int k = 0; k < pinCNTinObject; ++k) {
            PIN* pin1 = pin[k];
            (*pin1).from(net.pin[k]);
        }
    }

*/
inline void NET::copy(net_t* net) {
    sum_num1.from(net->sumNumL);
    sum_denom1.from(net->sumDenomL);
    sum_num2.from(net->sumNumR);
    sum_denom2.from(net->sumDenomR);
    min_x = net->min.x;
    min_y = net->min.y;
    max_x = net->max.x;
    max_y = net->max.y;
    for(int k = 0; k < pinCNTinObject; ++k) {
        net->pinArray[k].to(pin[k]);
    }
    free(net->pinArray);
}

inline void Cell_t::copy(CELLx* origin) {
    for(size_t i = 0; i < size; ++i) {
        den_pmin_x[i] = origin[i].den_pmin.x;
        den_pmin_y[i] = origin[i].den_pmin.y;
        den_pmax_x[i] = origin[i].den_pmax.x;
        den_pmax_y[i] = origin[i].den_pmax.y;
        //pmin_x[i] = origin[i].pmin.x;
        //pmin_y[i] = origin[i].pmin.y;
        pmax_x[i] = origin[i].pmax.x;
        pmax_y[i] = origin[i].pmax.y;
        half_size_x[i] = origin[i].half_den_size.x;
        half_size_y[i] = origin[i].half_den_size.y;
        scale[i] = origin[i].den_scal;
        flg[i] = origin[i].flg;
    }
}

inline void Cell_t::copyback(CELLx* destination) {
    for(size_t i = 0; i < size; ++i) {
        destination[i].den_pmin.x = den_pmin_x[i];
        destination[i].den_pmin.y = den_pmin_y[i];
        destination[i].den_pmax.x = den_pmax_x[i];
        destination[i].den_pmax.y = den_pmax_y[i];
        //destination[i].pmin.x = pmin_x[i];
        //destination[i].pmin.y = pmin_y[i];
        destination[i].pmax.x = pmax_x[i];
        destination[i].pmax.y = pmax_y[i];
        destination[i].half_den_size.x = half_size_x[i];
        destination[i].half_den_size.y = half_size_y[i];
        destination[i].den_scal = scale[i];
        destination[i].flg = flg[i];
    }
}
/*
inline void Cell_t::Copy(cell_t* origin) {
    for(size_t i = 0; i < size; ++i) {
        den_pmin_x[i] = origin[i].min.x;
        den_pmin_y[i] = origin[i].min.y;
        den_pmax_x[i] = origin[i].max.x;
        den_pmax_y[i] = origin[i].max.y;
        //pmin_x[i] = origin[i].pmin.x;
        //pmin_y[i] = origin[i].pmin.y;
        //pmax_x[i] = origin[i].pmax.x;
        //pmax_y[i] = origin[i].pmax.y;
        half_size_x[i] = origin[i].size.x;
        half_size_y[i] = origin[i].size.y;
        scale[i] = origin[i].scale;
        flg[i] = origin[i].type;
    }
}

inline void Cell_t::Copyback(cell_t* destination) {
    for(size_t i = 0; i < size; ++i) {
        destination[i].den_pmin.x = den_pmin_x[i];
        destination[i].den_pmin.y = den_pmin_y[i];
        destination[i].den_pmax.x = den_pmax_x[i];
        destination[i].den_pmax.y = den_pmax_y[i];
        //destination[i].pmin.x = pmin_x[i];
        //destination[i].pmin.y = pmin_y[i];
        //destination[i].pmax.x = pmax_x[i];
        //destination[i].pmax.y = pmax_y[i];
        destination[i].half_size.x = half_size_x[i];
        destination[i].half_size.y = half_size_y[i];
        destination[i].scale = scale[i];
        destination[i].flg = flg[i];
    }
}
*/
struct cpu_t {
    double total;
    double ip;
    double tgp;
    double cgp;
    double ns;
    double wlen;
    double bins;
    double density;
    double fft;
    double wgrad;
    double pgrad;
    double pre;
    double cost;
};

enum { STDCELLonly, MIXED };
extern std::string outputPATH;
extern int pinCNT;
extern int moduleCNT;
extern int gcell_cnt;

enum { FillerCell, StdCell, Macro };
extern int terminalCNT;
extern int netCNT;
extern int numNonRectangularNodes;
extern int totalShapeCount;

extern int g_rrr;
extern int STAGE;

extern vector< SHAPE > shapeStor;
extern dense_hash_map< string, vector< int > > shapeMap;

// STAGE's possible inner variables
enum {
    INITIAL_PLACE,  // INITIAL_PLACE
    tGP,            // TRIAL
    mGP3D,          // MIXED_SIZE_3D_GLOBAL_PLACE
    mGP2D,          // MIXED_SIZE_2D_GLOBAL_PLACE
    mLG3D,          // MACRO_LEGALIZATION
    cGP3D,          // STDCELL_ONLY_3D_GLOBAL_PLACE
    cGP2D,          // STDCELL_ONLY_2D_GLOBAL_PLACE
    DETAIL_PLACE
};

// these variable is required to have detailPlacer settings
extern int dpMode;
enum { None, FastPlace, NTUplace3, NTUplace4h };
extern string dpLocation;
extern string dpFlag;

extern bool isOnlyLGinDP;
extern bool isSkipPlacement;
extern bool hasDensityDP ;
extern prec densityDP ;


extern int placementStdcellCNT;
extern int gfiller_cnt;
extern int placementMacroCNT;
enum { CDFT, RDFT, DDCT };
extern int msh_yz;
extern int INPUT_FLG;
enum { ISPD05, ISPD06, ISPD06m, ISPD, MMS, SB, ETC, IBM };
extern int gmov_mac_cnt;
extern int cGP3D_buf_iter;
extern int cGP2D_buf_iter;
extern int mGP3D_tot_iter;
extern int mGP2D_tot_iter;
extern int cGP3D_tot_iter;
extern int cGP2D_tot_iter;
extern int flg_3dic;
extern int flg_3dic_io;
extern int numLayer;
extern int GP_DIM_ONE;
enum { NoneAdj, RandomAdj, SmartAdj };
enum { zCenterPlace, zCenterTierZero, zCenterTierMax };

enum { global_router_based, prob_ripple_based };

enum { FastDP, NTUpl3, NTUpl4h };
enum { MAX_PCNT_ORDER, MIN_TIER_ORDER, MAX_AREA_DIS_DIV };
enum { secondOrder, thirdOrder };
enum {
    potnPhase1 = 1,
    potnPhase2 = 2,
    potnPhase3 = 3,
    potnPhase4 = 4,
    potnPhase5 = 5,
    potnPhase6 = 6,
    potnPhase7 = 7,
    potnPhase8 = 8
};

extern int potnPhaseDS;
extern int bloating_max_count;
extern int bloatCNT;
extern int trial_iterCNT;
extern int mGP3D_iterCNT;
extern int mGP2D_iterCNT;
extern int cGP3D_iterCNT;
extern int cGP2D_iterCNT;

// routability
extern prec row_height;   // lutong
extern prec gcell_x_cof;  // lutong
extern prec gcell_y_cof;  // lutong
extern prec gcell_cof;
extern prec eq_pitch;
extern prec max_inflation_ratio;
extern prec inflation_ratio_coef;
extern prec edgeadj_coef;
extern prec pincnt_coef;
extern prec gRoute_pitch_scal;
extern prec inflation_threshold;
extern prec total_inflation_ratio;
extern prec h_max_inflation_ratio;
extern prec v_max_inflation_ratio;
extern prec total_inflatedNewAreaDelta;
extern prec currTotalInflation;
extern int inflation_cnt;
extern int inflation_max_cnt;
extern int routepath_maxdist;
extern prec inflation_area_over_whitespace;
extern prec curr_filler_area;
extern prec adjust_ratio;
extern bool is_inflation_h;
extern bool flg_noroute;

extern prec ALPHAmGP;
extern prec ALPHAcGP;
extern prec ALPHA;
extern prec BETAmGP;
extern prec BETAcGP;
extern prec BETA;

extern prec dampParam;
extern prec global_ovfl;  ////////////
extern prec maxALPHA;
extern prec ExtraWSfor3D;
extern prec MaxExtraWSfor3D;
extern prec rowHeight;
extern prec SITE_SPA;
extern prec layout_area;
extern double tot_HPWL;
extern prec tx_HPWL;
extern prec ty_HPWL;
extern prec tz_HPWL;
extern prec tot_overlap;
extern prec total_std_area;
extern prec total_std_den;
extern prec total_modu_area;
extern prec inv_total_modu_area;
extern prec total_cell_area;
extern prec curr_cell_area;  // lutong


extern prec total_term_area;
extern prec total_move_available_area;
extern prec total_filler_area;
extern prec total_PL_area;
extern long total_termPL_area;
extern long total_WS_area;

extern prec curr_WS_area;  // lutong
extern prec filler_area;
extern prec target_cell_den;
extern prec target_cell_den_orig;  // lutong
extern prec total_macro_area;
extern prec ignoreEdgeRatio;
extern prec grad_stp;
extern prec gsum_phi;
extern prec gsum_ovfl;
extern prec gsum_ovf_area;
extern prec overflowMin;
extern prec overflowMin_initial;
extern prec mGP3D_opt_phi_cof;
extern prec mGP2D_opt_phi_cof;
extern prec cGP3D_opt_phi_cof;
extern prec cGP2D_opt_phi_cof;
extern prec inv_RAND_MAX;
extern prec theta;
extern prec dp_margin_per_tier;
extern prec stn_weight;  // lutong
extern prec opt_w_cof;   // lutong

extern unsigned extPt1_2ndOrder;
extern unsigned extPt2_2ndOrder;
extern unsigned extPt1_1stOrder;
extern unsigned extPt2_1stOrder;
extern unsigned extPt3_1stOrder;

extern char gbch_dir[BUF_SZ];
extern char gbch_aux[BUF_SZ];
extern char gbch[BUF_SZ];
extern char gGP_dir[BUF_SZ];
extern char gGP_pl[BUF_SZ];
extern char gIP_pl[BUF_SZ];
extern char gGP_pl_file[BUF_SZ];
extern char gmGP2D_pl[BUF_SZ];
extern char gGP3_pl[BUF_SZ];
extern char gLG_pl[BUF_SZ];
extern char gDP_log[BUF_SZ];
extern char gDP_pl[BUF_SZ];
extern char gDP_tmp[BUF_SZ];
extern char gDP2_pl[BUF_SZ];
extern char gDP3_pl[BUF_SZ];
extern char gGR_dir[BUF_SZ];
extern char gGR_log[BUF_SZ];
extern char gGR_tmp[BUF_SZ];
extern char gFinal_DP_pl[BUF_SZ];
extern char gTMP_bch_dir[BUF_SZ];
extern char gTMP_bch_aux[BUF_SZ];
extern char gTMP_bch_nodes[BUF_SZ];
extern char gTMP_bch_nets[BUF_SZ];
extern char gTMP_bch_wts[BUF_SZ];
extern char gTMP_bch_pl[BUF_SZ];
extern char gTMP_bch_scl[BUF_SZ];
extern char sing_fn_aux[BUF_SZ];
extern char sing_fn_nets[BUF_SZ];
extern char sing_fn_nodes[BUF_SZ];
extern char sing_fn_pl[BUF_SZ];
extern char sing_fn_wts[BUF_SZ];
extern char sing_fn_scl[BUF_SZ];
extern char fn_bch_IP[BUF_SZ];
extern char fn_bch_GP[BUF_SZ];
extern char fn_bch_GP2[BUF_SZ];
extern char fn_bch_GP3[BUF_SZ];
extern char fn_bch_mac_LG[BUF_SZ];
extern char fn_bch_DP[BUF_SZ];
extern char bench_aux[BUF_SZ];
extern char dir_bnd[BUF_SZ];
extern char global_router[1023];
extern char output_dir[BUF_SZ];
extern char currentDir[BUF_SZ];

extern string sourceCodeDir;

extern bool isBloatingBegin;
extern bool isRoutabilityInit;
extern bool isTrial;
extern bool isFirst_gp_opt;
extern bool DEN_ONLY_PRECON;
extern int orderHPWL;

extern vector< std::pair< int, prec > > inflationList;
extern vector< std::pair< prec, prec > > trial_HPWLs;
extern vector< prec > trial_POTNs;
extern vector< std::pair< prec, prec > > hpwlEPs_1stOrder;
extern vector< std::pair< prec, prec > > hpwlEPs_2ndOrder;
extern vector< prec > potnEPs;
extern std::map< string, vector< int > > routeBlockageNodes;
extern int nXgrids, nYgrids, nMetLayers;
extern vector< int > verticalCapacity;
extern vector< int > horizontalCapacity;
extern vector< prec > minWireWidth;
extern vector< prec > minWireSpacing;
extern vector< prec > viaSpacing;
extern vector< tuple< int, int, int, int, int, int, int > >
    edgeCapAdj;
extern prec gridLLx, gridLLy;
extern prec tileWidth, tileHeight;
extern prec blockagePorosity;

extern RECT cur_rect;
extern PIN *pinInstance;
extern MODULE *moduleInstance;
extern CELLx *gcell_st;
extern TERM *terminalInstance;
extern NET *netInstance;

// structure for *.scl
extern ROW *row_st;
extern int row_cnt;

// ??
extern PLACE *place_st;
extern PLACE *place_backup_st; // why?
extern int place_st_cnt;

// ??
extern PLACE place;
extern PLACE place_backup; //  why?

extern FPOS avgCellSize;
extern FPOS avgTerminalSize;

extern FPOS term_pmax;
extern FPOS term_pmin;

extern FPOS filler_size;

extern FPOS zeroFPoint;
extern POS zeroPoint;

extern POS msh;

// global Space Min/Max
extern FPOS gmin;
extern FPOS gmax;

extern FPOS gwid;
extern TIER *tier_st;
extern POS dim_bin;
extern POS dim_bin_mGP2D;
extern POS dim_bin_cGP2D;

extern FPOS grow_pmin;
extern FPOS grow_pmax;

///////////////////////////////////////////////////////////////////////////
/*  ARGUMENTS: main.cpp                                                  */
///////////////////////////////////////////////////////////////////////////
extern string bmFlagCMD;
extern string auxCMD;
extern string lefCMD;
extern string defCMD;
extern string verilogCMD;
extern string outputCMD;

extern string benchName;

extern int numThread;
enum class InputMode { bookshelf, lefdef };
extern InputMode inputMode;

extern string denCMD;
extern string bxMaxCMD;
extern string byMaxCMD;
extern string bzMaxCMD;
extern string overflowMinCMD;
extern string pcofmaxCMD;
extern string pcofminCMD;
extern string racntiCMD;       // lutong
extern string maxinflCMD;      // lutong
extern string inflcoefCMD;     // lutong
extern string filleriterCMD;
extern string refdWLCMD;
extern int conges_eval_methodCMD;  // grouter | prob
extern bool isPreplace;
extern bool isPlot;
extern bool plotCellCMD;
extern bool plotMacroCMD;
extern bool plotDensityCMD;
extern bool plotFieldCMD;
extern bool constraintDrivenCMD;
extern bool routabilityCMD;
extern bool lambda2CMD;
extern bool dynamicStepCMD;
extern bool thermalAwarePlaceCMD;
extern bool onlyGlobalPlaceCMD;
extern bool isARbyUserCMD;
extern bool stnCMD;  // lutong
extern bool trialRunCMD;
extern bool autoEvalRC_CMD;
extern     std::fstream fsWk;
extern cpu_t profile;

//////////////////////////////////////////////////////////////////////////
// Defined in main.cpp ///////////////////////////////////////////////////
void init();
void initialPlacement_main(void);
void trial_main(void);
void free_trial_mallocs(void);
void tmGP3DglobalPlacement_main(void);
void tmGP2DglobalPlacement_main(void);
void tcGP3DglobalPlacement_main(void);
void tcGP2DglobalPlacement_main(void);
void mGP3DglobalPlacement_main(void);
void mGP2DglobalPlacement_main(void);
void cGP3DglobalPlacement_main(void);
void cGP2DglobalPlacement_main(void);
void macroLegalization_main(void);
void calcRef_dWL(void);
void findMinHPWLduringTrial(void);
void findMaxHPWLduringTrial(void);
void get2ndOrderEP1(void);
void get2ndOrderEP2(void);
bool isLargestGapHPWL_1stEP(void);
void sort2ndOrderEPs(void);
void printTrend(void);
void get2ndOrderDiff_HPWL_LinearHPWLtrend(void);
void get1stOrderDiff_HPWL_LinearHPWLtrendforThird();
void get1stOrderDiff_HPWL_LinearHPWLtrendforSecond(void);
void get1stOrder_ExtremePointsforThird(void);
void get1stOrder_ExtremePointsforSecond(void);
prec get_abs(prec a);
void store2ndOrder_ExtremePoints(void);
void store_POTNs(void);
void store1stOrder_ExtremePointsforThird(void);
void store1stOrder_ExtremePointsforSecond(void);
void reassign_trial_2ndOrder_lastEP(prec);
void printEPs(void);

void printUsage(void);
void initArgument(int, char **);
void calcTSVweight(void);
bool argument(int, char **);
void printCMD(int, char **);
bool criticalArgumentError(void);

string getexepath();
int pos_eqv(struct POS p1, struct POS p2);
int prec_eqv(prec x, prec y);
int prec2int(prec a);
unsigned prec2unsigned(prec a);
int find_non_zero(prec *a, int cnt);
int prec_le(prec x, prec y);
int prec_ge(prec x, prec y);
int prec_lt(prec x, prec y);
int prec_gt(prec x, prec y);
int HPWL_count(void);
void overlap_count(int iter);
void update_net_by_pin(void);
void time_start(double *time_cost);
void time_end(double *time_cost);
void time_calc(double last_time, double *curr_time, double *time_cost);
void itoa(int n, char k[]);

inline int dge(prec a, prec b) {
    return (a > b || a == b) ? 1 : 0;
}
inline int dle(prec a, prec b) {
    return (a < b || a == b) ? 1 : 0;
}

void OR_opt(void);
void rdft2dsort(int, int, int, prec **);

//void call_DP(void);
//void call_FastDP(void);
//void call_NTUpl3(void);
//void call_NTUpl4h(void);
//void call_DP_StdCell_NTUpl4h(void);
void read_macro(char *fn);

FPOS fp_mul(struct FPOS a, struct FPOS b);
inline FPOS fp_add(struct FPOS a, struct FPOS b) {
    struct FPOS c = zeroFPoint;
    c.x = a.x + b.x;
    c.y = a.y + b.y;
    if(flg_3dic)
        c.z = a.z + b.z;
    return c;
}
FPOS fp_add_abs(struct FPOS a, struct FPOS b);
inline FPOS fp_scal(prec s, struct FPOS a) {
    struct FPOS c = a;
    c.x *= s;
    c.y *= s;
    if(flg_3dic)
        c.z *= s;
    return c;
}
FPOS fp_subt(struct FPOS a, struct FPOS b);
FPOS fp_subt_const(struct FPOS a, prec b);
prec fp_sum(struct FPOS a);
FPOS fp_exp(struct FPOS a);
prec fp_product(struct FPOS a);

int p_product(struct POS a);
int p_max(struct POS a);

FPOS fp_min2(struct FPOS a, struct FPOS b);
FPOS fp_max2(struct FPOS a, struct FPOS b);
FPOS fp_div(struct FPOS a, struct FPOS b);
FPOS fp_rand(void);
FPOS p2fp(struct POS a);
POS fp2p_floor(struct FPOS a);
POS fp2p_ceil(struct FPOS a);
FPOS fp_inv(struct FPOS a);

void tier_assign(int);
void tier_assign_with_macro(void);
void tier_assign_without_macro(void);
void find_close_tier(prec z, struct T0 *t0_st, int *z_st);
int prec_cmp(const void *a, const void *b);
int max_pinCNTinObject_cmp(const void *a, const void *b);
int min_tier_cmp(const void *a, const void *b);
int max_area_dis_div_cmp(const void *a, const void *b);

void call_FastDP_tier(char *tier_dir, char *tier_aux, char *tier_pl);
//void call_NTUpl3_tier(char *tier_dir, char *tier_aux, char *tier_pl);
//void call_NTUpl4h_tier(char *tier_dir, char *tier_aux, char *tier_pl);

void preprocess_SB_inputs(char *tier_dir);
void postprocess_SB_inputs(char *tier_dir);
void init_tier(void);
void tot_area_comp(void);
void calc_tier_WS(void);
void post_mac_tier(void);
void pre_mac_tier(void);
inline prec getStepSizefromEPs(prec hpwl, prec hpwlEP, prec hpwlSTD,
                               prec basePCOF, prec baseRange) {
    return min(
        basePCOF + baseRange,
        basePCOF +
            baseRange * (get_abs(hpwl - hpwlSTD) / get_abs(hpwlEP - hpwlSTD)));
}

// writing Bookshelf function 
void WriteBookshelf();
void CallDetailPlace();

// useful function

// return Common Area
// between Rectangle A and Rectangle B.
// type : casted long from prec
inline long lGetCommonAreaXY( FPOS aLL, FPOS aUR, FPOS bLL, FPOS bUR ) {
    long xLL = max( aLL.x, bLL.x ),
            yLL = max( aLL.y, bLL.y ),
            xUR = min( aUR.x, bUR.x ),
            yUR = min( aUR.y, bUR.y );

    if( xLL >= xUR || yLL >= yUR ) {
        return 0;
    }
    else {
        return (xUR - xLL) * (yUR - yLL);
    }
}

// return Common Area
// between Rectangle A and Rectangle B.
// type : integer
inline int iGetCommonAreaXY( POS aLL, POS aUR, POS bLL, POS bUR ) {
    int xLL = max( aLL.x, bLL.x ),
            yLL = max( aLL.y, bLL.y ),
            xUR = min( aUR.x, bUR.x ),
            yUR = min( aUR.y, bUR.y );

    if( xLL >= xUR || yLL >= yUR ) {
        return 0;
    }
    else {
        return (xUR - xLL) * (yUR - yLL);
    }
}

// return Common Area
// between Rectangle A and Rectangle B.
// type : prec  
inline prec pGetCommonAreaXY( FPOS aLL, FPOS aUR, FPOS bLL, FPOS bUR ) {
    prec xLL = max( aLL.x, bLL.x ),
            yLL = max( aLL.y, bLL.y ),
            xUR = min( aUR.x, bUR.x ),
            yUR = min( aUR.y, bUR.y );

    if( xLL >= xUR || yLL >= yUR ) {
        return 0;
    }
    else {
        return (xUR - xLL) * (yUR - yLL);
    }
}

#endif
