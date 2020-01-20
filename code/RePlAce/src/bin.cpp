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

#include <error.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctgmath>
#include <string>
#include <omp.h>

#include "bin.h"
#include "charge.h"
#include "fft.h"
#include "global.h"
#include "macro.h"
#include "opt.h"
#include "plot.h"
#include "wlen.h"

using std::to_string;

int tot_bin_cnt;
int bin_cnt;
prec bin_area;
prec inv_bin_area;
prec global_macro_area_scale;

BIN *bin_mat;
static POS *bin_st;
static prec *bin_share_x_st;
static prec *bin_share_y_st;
static prec *bin_share_z_st;
static prec *bin_share_st;

FPOS bin_org;
FPOS bin_stp;
FPOS inv_bin_stp;
FPOS half_bin_stp;
FPOS bin_stp_mGP2D;
FPOS bin_stp_cGP2D;
POS max_bin;

static BIN **bin_mat_st;
#define BNDBDY(i, x, m)\
    float b##i = INT_DOWN((cells->den_p##m##_##x[k] - bin_org.x) * stp.x);\
    if(b##i < 0) b##i = 0;\
    if(b##i > dim.x - 1) b##i = dim.x - 1;\
    cells->b##i##_##x[k] = b##i;

//BIN **bin_list;
//int bin_list_cnt;


// numLayer>1 : 3DIC cases.
//
// STAGE : mGP2D || cGP2D
// STAGE : INIT_PLACEMENT -> Totally not-working 
POS get_bin_idx(int idx) {
    POS p;
    if(numLayer > 1) {
        p.x = idx / (max_bin.y * max_bin.z);
        p.y = (idx % (max_bin.y * max_bin.z)) / max_bin.z;
        p.z = idx % max_bin.z;
    }
    else {
        if(STAGE == mGP2D) {
            cout << "this is Mixed GP" << endl;
            p.x = idx / dim_bin_mGP2D.y;
            p.y = idx % dim_bin_mGP2D.y;
        }
        else if(STAGE == cGP2D) {
            cout << "this is Original GP" << endl;
            p.x = idx / dim_bin_cGP2D.y;
            p.y = idx % dim_bin_cGP2D.y;
        }
    }

    return p;
}

// update below variable --> Totally not working
//
//
// max_bin -> tot_bin_cnt
// bin_stp -> bin_area -> inv_bin_area
//
// half_bin_stp, inv_bin_stp
//
//
// bin_mat, 
//
// bin_st,
// bin_share_x_st,
// bin_share_y_st,
// bin_share_z_st,  bin_share_st
//
// place_st


void bin_init() {
    POS p;
    bin_area = 0;

    max_bin = msh;
   
    // update bin_step 
    bin_stp.x = place.cnt.x / (prec)max_bin.x;
    bin_stp.y = place.cnt.y / (prec)max_bin.y;
    bin_stp.z = place.cnt.z / (prec)max_bin.z;

    half_bin_stp.x = bin_stp.x * 0.5;
    half_bin_stp.y = bin_stp.y * 0.5;
    half_bin_stp.z = bin_stp.z * 0.5;

    inv_bin_stp.x = 1.0 / (prec)bin_stp.x;
    inv_bin_stp.y = 1.0 / (prec)bin_stp.y;
    inv_bin_stp.z = 1.0 / (prec)bin_stp.z;


    if(flg_3dic) {
        printf("INFO:  Boundary Bins Max_b_x=%d, Max_b_y=%d, Max_b_z=%d\n",
               max_bin.x, max_bin.y, max_bin.z);
        printf("INFO:  Bin Step X=%.4lf, Y=%.4lf, Z=%.4lf\n\n", bin_stp.x,
               bin_stp.y, bin_stp.z);
    }
    else {
        printf("INFO:  Boundary Bins Max_b_x=%d, Max_b_y=%d\n", max_bin.x,
               max_bin.y);
        printf("INFO:  Bin Step X=%.4lf, Y=%.4lf\n\n", bin_stp.x, bin_stp.y);
    }

    bin_org = place.org;

    tot_bin_cnt = p_product(max_bin);
    bin_area = fp_product(bin_stp);
    inv_bin_area = 1.0 / bin_area;

    // update bin_mat's information
    bin_mat = (BIN *)mkl_malloc(sizeof(BIN) * tot_bin_cnt, 64);

    // for each allocated bin_mat
    for(int i = 0; i < tot_bin_cnt; i++) {
        // ???????? not working 
        p = get_bin_idx(i);
        exit(1);
//        cout << i << ", " ;
//        p.Dump();

        BIN *binp = &bin_mat[i];

        binp->den = 0;
        binp->e.SetZero();
        binp->phi = 0;

        binp->pmin = fp_add(bin_org, fp_mul(bin_stp, p2fp(p)));
        binp->pmax = fp_add(binp->pmin, bin_stp);
        binp->center = fp_add(binp->pmin, fp_scal(0.5, bin_stp));

        binp->p = p;

        binp->term_area = 0;
//        binp->term_area2 = 0;
        binp->flg = 0;
        binp->virt_area = 0;
        binp->pl_area = 0;

        for(int k = 0; k < place_st_cnt; k++) {
            PLACE *pl = &place_st[k];
            binp->pl_area += pGetCommonAreaXY(pl->org, pl->end, binp->pmin, binp->pmax);
        }

        binp->virt_area = (bin_area - binp->pl_area) * global_macro_area_scale;
//        binp->dump(to_string(i));
    }


    // calculate overlap area in terminalInstance.
    for(int i = 0; i < terminalCNT; i++) {
        TERM *curTerminal = &terminalInstance[i];

        // ignore for TerminalNI
        if(curTerminal->isTerminalNI) {
            continue;
        }
    
        for(int j = 0; j < place_st_cnt; j++) {
            PLACE *pl = &place_st[j];

            RECT rect;
            get_common_rect(curTerminal->pmin, curTerminal->pmax, pl->org, pl->end, 
                            &rect.pmin, &rect.pmax); 
            
            // if cannot find any common rectangles.
            if(rect.pmin.x < 0.0 || rect.pmin.y < 0.0 || rect.pmax.x < 0.0 ||
               rect.pmax.y < 0.0)
                continue;

            POS bin_pmin = get_bin_pt_from_point(rect.pmin);
            POS bin_pmax = get_bin_pt_from_point(rect.pmax);

            p.Set(-1);
    
            // bin_mat's information 
            while(idx_in_bin_rect(&p, bin_pmin, bin_pmax)) {
                BIN *binp = get_bin_from_idx(p);
                binp->term_area += pGetCommonAreaXY( rect.pmin, rect.pmax, binp->pmin, binp->pmax) 
                                   * global_macro_area_scale;
            }
        }
    }

    for(int i = 0; i < terminalCNT; i++) {
        TERM *curTerminal = &terminalInstance[i];

        if(is_IO_block(curTerminal)) {
            curTerminal->IO = 1;
        }
        else {
            curTerminal->IO = 0;
        }
    }

    bin_st = (POS *)mkl_malloc(sizeof(POS) * tot_bin_cnt, 64);
    bin_share_x_st = (prec *)mkl_malloc(sizeof(prec) * max_bin.x, 64);
    bin_share_y_st = (prec *)mkl_malloc(sizeof(prec) * max_bin.y, 64);
    bin_share_z_st = (prec *)mkl_malloc(sizeof(prec) * max_bin.z, 64);
    bin_share_st = (prec *)mkl_malloc(sizeof(prec) * tot_bin_cnt, 64);
}

//
// mainly update
//
// BIN matrix, 
// TIER, bin_pmin, bin_pmax
// TIER->bin_mat update
//
//
void bin_init_2D(int STAGE) {

    if(STAGE == cGP2D) {
        dim_bin_cGP2D.x = dim_bin_cGP2D.y = 0;
    }
    else if(STAGE == mGP2D) {
        dim_bin_mGP2D.x = dim_bin_mGP2D.y = 0;
    }

    for(int z = 0; z < numLayer; z++) {
        TIER *tier = &tier_st[z];
        prec avg_modu_area = 1.0 * tier->modu_area / tier->modu_cnt;
        prec ideal_bin_area = avg_modu_area / target_cell_den;
        int ideal_bin_cnt = INT_CONVERT(tier->area / ideal_bin_area);

        bool isUpdate = false;
        for(int i = 1; i <= 10; i++) {
            if((2 << i) * (2 << i) <= ideal_bin_cnt &&
               (2 << (i + 1)) * (2 << (i + 1)) > ideal_bin_cnt) {
                tier->dim_bin.x = tier->dim_bin.y = 2 << i;
                isUpdate = true;
                break;
            }
        }

        if(!isUpdate) {
            tier->dim_bin.x = tier->dim_bin.y = 1024;
        }

        if(STAGE == mGP2D) {
            if(dim_bin_mGP2D.x < tier->dim_bin.x)
                dim_bin_mGP2D.x = tier->dim_bin.x;
            if(dim_bin_mGP2D.y < tier->dim_bin.y)
                dim_bin_mGP2D.y = tier->dim_bin.y;
        }
        else if(STAGE == cGP2D) {
            if(dim_bin_cGP2D.x < tier->dim_bin.x)
                dim_bin_cGP2D.x = tier->dim_bin.x;
            if(dim_bin_cGP2D.y < tier->dim_bin.y)
                dim_bin_cGP2D.y = tier->dim_bin.y;
        }
    }

    if(STAGE == mGP2D) {
        // LW 05/30/17
        // IK 05/08/17
        // if (routabilityCMD == true) {
        //    if (placementMacroCNT == 0) {
        //        dim_bin_mGP2D.x *= 2;
        //        dim_bin_mGP2D.y *= 2;
        //    }
        //}
        bin_stp_mGP2D.x = place.cnt.x / (prec)dim_bin_mGP2D.x;
        bin_stp_mGP2D.y = place.cnt.y / (prec)dim_bin_mGP2D.y;
        bin_stp_mGP2D.z = TIER_DEP;
        printf("INFO:  dim_bin_mGP2D.(x,y) = (%d, %d)\n", dim_bin_mGP2D.x,
               dim_bin_mGP2D.y);
    }
    else if(STAGE == cGP2D) {
        // LW 05/30/17
        // IK 05/08/17
        // if (routabilityCMD == true && tier_st[0].modu_den < 0.32) {
        //    if (placementMacroCNT == 0) {
        //        dim_bin_cGP2D.x *= 2;
        //        dim_bin_cGP2D.y *= 2;
        //    }
        //}
        bin_stp_cGP2D.x = place.cnt.x / (prec)dim_bin_cGP2D.x;
        bin_stp_cGP2D.y = place.cnt.y / (prec)dim_bin_cGP2D.y;
        bin_stp_cGP2D.z = TIER_DEP;
        printf("INFO:  dim_bin_cGP2D.(x,y) = (%d, %d)\n", dim_bin_cGP2D.x,
               dim_bin_cGP2D.y);
    }
    
    bin_mat_st = (BIN **)mkl_malloc(sizeof(BIN *) * numLayer, 64);
    for(int z = 0; z < numLayer; z++) {
        TIER *tier = &tier_st[z];

        if(STAGE == cGP2D)
            tier->dim_bin = dim_bin_cGP2D;
        else if(STAGE == mGP2D)
            tier->dim_bin = dim_bin_mGP2D;

//        tier->dim_bin.Dump("tier->dim_bin");

        tier->tot_bin_cnt = tier->dim_bin.x * tier->dim_bin.y;

        tier->bin_stp.x = tier->size.x / tier->dim_bin.x;
        tier->bin_stp.y = tier->size.y / tier->dim_bin.y;
        tier->bin_stp.z = TIER_DEP;

        tier->half_bin_stp.x = 0.5 * tier->bin_stp.x;
        tier->half_bin_stp.y = 0.5 * tier->bin_stp.y;
        tier->half_bin_stp.z = 0.5 * tier->bin_stp.z;

        tier->inv_bin_stp.x = 1.0 / tier->bin_stp.x;
        tier->inv_bin_stp.y = 1.0 / tier->bin_stp.y;
        tier->inv_bin_stp.z = 1.0 / tier->bin_stp.z;

        tier->bin_off.SetZero();
        tier->bin_org.x = tier->pmin.x - tier->bin_off.x;
        tier->bin_org.y = tier->pmin.y - tier->bin_off.y;
        tier->bin_org.z = tier->pmin.z - tier->bin_off.z;

        tier->bin_area = tier->bin_stp.x * tier->bin_stp.y * tier->bin_stp.z;

        tier->inv_bin_area = 1.0 / tier->bin_area;
        tier->tot_bin_area = tier->bin_area * tier->tot_bin_cnt;

        bin_mat_st[z] = (BIN *)mkl_malloc(sizeof(BIN) * tier->tot_bin_cnt, 64);
        tier->bin_mat = bin_mat_st[z];

        // for each allocated bin_mat_st..
        for(int i = 0; i < tier->tot_bin_cnt; i++) {
            POS p( i / tier->dim_bin.y, i % tier->dim_bin.y, z);

            BIN *bp = &tier->bin_mat[i];

            bp->den = 0;
            bp->e = zeroFPoint;
            ////igkang
            // bp->e_local = zeroFPoint;
            bp->phi = 0;
            bp->pmin.x = tier->bin_org.x + (prec)p.x * tier->bin_stp.x;
            bp->pmin.y = tier->bin_org.y + (prec)p.y * tier->bin_stp.y;
            bp->pmin.z = (prec)z * TIER_DEP;

            bp->pmax.x = bp->pmin.x + tier->bin_stp.x;
            bp->pmax.y = bp->pmin.y + tier->bin_stp.y;
            bp->pmax.z = (prec)(z + 1) * TIER_DEP;

//            p.Dump("p");
//            bp->pmin.Dump("bp->pmin");
//            bp->pmax.Dump("bp->pmax");
//            cout << endl;

            bp->center.x = bp->pmin.x + 0.5 * tier->bin_stp.x;
            bp->center.y = bp->pmin.y + 0.5 * tier->bin_stp.y;
            bp->center.z = bp->pmin.z + 0.5 * tier->bin_stp.z;

            bp->p = p;

            bp->cell_area = 0;
            bp->cell_area2 = 0;

            bp->term_area = 0;
            bp->flg = 0;

            bp->virt_area = 0;

            bp->pl_area = 0;

            for(int k = 0; k < place_st_cnt; k++) {
                PLACE* pl = &place_st[k];
//                bp->pl_area += get_common_area(pl->org, pl->end, bp->pmin, bp->pmax);
                bp->pl_area += pGetCommonAreaXY(pl->org, pl->end, bp->pmin, bp->pmax);
            }

            bp->virt_area =
                (tier->bin_area - bp->pl_area) * global_macro_area_scale;
        }
    }


    // below step is for Update Terminal's Area.
    
    // for each terminal
    for(int i = 0; i < terminalCNT; i++) {
        TERM* term = &terminalInstance[i];

        // skip for IO pins.
        if(term->isTerminalNI) {
            // term->pmin.Dump("term->pmin");
            // term->pmax.Dump("term->pmax");
            continue;   
        }

        int z = term->tier;
        TIER *tier = &tier_st[z];

        // there is no shapes 
        if( shapeMap.find( term-> name ) == shapeMap.end() ) {
            UpdateTerminalArea(tier, &(term->pmin), &(term->pmax) ); 
        }
        // there are shapes
        else {
            for(auto& curIdx : shapeMap[term->name]) {
                prec llx = shapeStor[curIdx].llx, 
                     lly = shapeStor[curIdx].lly,
                     width = shapeStor[curIdx].width,
                     height = shapeStor[curIdx].height;
                FPOS tmpMin(llx, lly, 0), tmpMax(llx + width, lly + height, 1);

                UpdateTerminalArea(tier, &tmpMin, &tmpMax );
            }
        }
//        UpdateTerminalArea(tier, &(term->pmin), &(term->pmax));
    }

    if(STAGE == cGP2D) {
        // for each macro cell
        for(int i = 0; i < macro_cnt; i++) {
            MODULE* mac = macro_st[i];

            int z = mac->tier;
            TIER *tier = &tier_st[z];
            UpdateTerminalArea(tier, &(mac->pmin), &(mac->pmax));
        }
    }

    // check
//    TIER* tier = &tier_st[0];
//    for(int i=0; i<tier->tot_bin_cnt; i++) {
//        BIN* bp = &tier->bin_mat[i];
//        cout << bp->p.x << ", " << bp->p.y << ", " << bp->p.z << " :" << bp->term_area << endl;
//    }
}


// update bin_mat_st's term_area variable
void UpdateTerminalArea( TIER* tier, FPOS* pmin, FPOS* pmax ) {
    
    // for each place Instance
    for(int j = 0; j < place_st_cnt; j++) {
        PLACE* pl = &place_st[j];

        RECT rect;
        get_common_rect(*pmin, *pmax, pl->org, pl->end,
                &rect.pmin, &rect.pmax);

        // if cannot find any common rectangles.
        if(rect.pmin.x < 0.0 || rect.pmin.y < 0.0 || rect.pmax.x < 0.0 ||
                rect.pmax.y < 0.0)
            continue;

        POS bin_pmin( 
            INT_DOWN((rect.pmin.x - tier->bin_org.x) * tier->inv_bin_stp.x),
            INT_DOWN((rect.pmin.y - tier->bin_org.y) * tier->inv_bin_stp.y));
        
        POS bin_pmax( 
            INT_DOWN((rect.pmax.x - tier->bin_org.x) * tier->inv_bin_stp.x),
            INT_DOWN((rect.pmax.y - tier->bin_org.y) * tier->inv_bin_stp.y));

        bin_pmin.SetXYProjection( POS(0, 0), POS(tier->dim_bin.x-1, 
                    tier->dim_bin.y-1)); 
        bin_pmax.SetXYProjection( POS(0, 0), POS(tier->dim_bin.x-1, 
                    tier->dim_bin.y-1));

        for(int x = bin_pmin.x; x <= bin_pmax.x; x++) {
            for(int y = bin_pmin.y; y <= bin_pmax.y; y++) {
                BIN* curBin = &tier->bin_mat[x * tier->dim_bin.y + y];

                // due to precision
//                long terminalOverlap = lGetCommonAreaXY(rect.pmin, rect.pmax,
//                        curBin->pmin, curBin->pmax);
                prec terminalOverlap = pGetCommonAreaXY(rect.pmin, rect.pmax,
                        curBin->pmin, curBin->pmax);

                curBin->term_area +=
                    terminalOverlap * global_macro_area_scale;
            }
        }
    }
}


void bin_delete(void) {
    mkl_free(bin_mat);

    bin_mat = NULL;

    mkl_free(bin_st);
    mkl_free(bin_share_x_st);
    mkl_free(bin_share_y_st);
    mkl_free(bin_share_z_st);
    mkl_free(bin_share_st);
}

// min_a <= a && a <= max_a
//int wthin(prec a, prec min_a, prec max_a) {
//    return dle(a, max_a) && dge(a, min_a);
//}

/*
FPOS get_int_rgn(FPOS min1, FPOS max1, FPOS min2, FPOS max2) {
    FPOS a;

    a.x = get_int_line(min1.x, max1.x, min2.x, max2.x);
    a.y = get_int_line(min1.y, max1.y, min2.y, max2.y);

    if(flg_3dic) {
        a.z = get_int_line(min1.z, max1.z, min2.z, max2.z);
    }

    return a;
}

// get common Length 
prec get_int_line(prec aLeft, prec aRight, prec bLeft, prec bRight) {
    prec d12 = 0, d21 = 0;

    if(wthin(aLeft, bLeft, bRight) && wthin(aRight, bLeft, bRight)) {
        return aRight - aLeft;
    }
    else if(wthin(bLeft, aLeft, aRight) && wthin(bRight, aLeft, aRight)) {
        return bRight - bLeft;
    }
    else {
        d12 = aRight - bLeft;
        d21 = bRight - aLeft;

        if(dge(d12, 0) && dge(d21, 0))
            return min(d12, d21);
        else
            return 0;
    }
}
*/

// return BIN's point
POS get_bin_pt_from_point(FPOS fp) {
    POS p;
    p.x = INT_DOWN((fp.x - bin_org.x) / bin_stp.x);
    p.y = INT_DOWN((fp.y - bin_org.y) / bin_stp.y);

    if(flg_3dic)
        p.z = INT_DOWN((fp.z - bin_org.z) / bin_stp.z);

    // projection 
    legal_bin_idx(&p);

    return p;
}

// projection 
void legal_bin_idx(POS *p) {
    if(p->x < 0)
        p->x = 0;
    if(p->y < 0)
        p->y = 0;
    if(flg_3dic && p->z < 0)
        p->z = 0;

    if(p->x >= max_bin.x)
        p->x = max_bin.x - 1;
    if(p->y >= max_bin.y)
        p->y = max_bin.y - 1;
    if(flg_3dic && p->z >= max_bin.z)
        p->z = max_bin.z - 1;
}

int point_in_range(prec a, prec min_a, prec max_a) {
    return a >= min_a - Epsilon && a <= max_a + Epsilon;
}

FPOS valid_coor4(FPOS center, FPOS obj_size) {
    int i = 0;
    prec xyz_dis = 0;
    prec min_dis = dp_wlen_weight.x * place.cnt.x +
                   dp_wlen_weight.y * place.cnt.y +
                   dp_wlen_weight.z * place.cnt.z;
    POS flg = zeroPoint, obj_flg = zeroPoint;
    FPOS new_center = center;
    FPOS dorg = zeroFPoint, dend = zeroFPoint;
    FPOS dis = zeroFPoint, obj_dis = zeroFPoint;
    PLACE *place0 = NULL;
    FPOS half_size = zeroFPoint;

    half_size.x = obj_size.x * 0.5;
    half_size.y = obj_size.y * 0.5;
    half_size.z = obj_size.z * 0.5;

    for(i = 0; i < place_st_cnt; i++) {
        place0 = &place_st[i];

        flg.x = point_in_range(center.x, place0->org.x, place0->end.x);

        flg.y = point_in_range(center.y, place0->org.y, place0->end.y);

        flg.z = point_in_range(center.z, place0->org.z, place0->end.z);

        dorg.x = fabs(place0->org.x - center.x);
        dorg.y = fabs(place0->org.y - center.y);
        dorg.z = fabs(place0->org.z - center.z);

        dend.x = fabs(place0->end.x - center.x);
        dend.y = fabs(place0->end.y - center.y);
        dend.z = fabs(place0->end.z - center.z);

        if(flg.x && flg.y && flg.z) {
            if(dorg.x < dend.x && dorg.x < half_size.x) {
                new_center.x += half_size.x - dorg.x;
            }
            else if(dorg.x > dend.x && dend.x < half_size.x) {
                new_center.x -= half_size.x - dend.x;
            }

            if(dorg.y < dend.y && dorg.y < half_size.y) {
                new_center.y += half_size.y - dorg.y;
            }
            else if(dorg.y > dend.y && dend.y < half_size.y) {
                new_center.y -= half_size.y - dend.y;
            }

            if(dorg.z < dend.z && dorg.z < half_size.z) {
                new_center.z += half_size.z - dorg.z;
            }
            else if(dorg.z > dend.z && dend.z < half_size.z) {
                new_center.z -= half_size.z - dend.z;
            }

            return new_center;
        }
        else if(flg.x && flg.y) {
            if(dorg.x < dend.x && dorg.x < half_size.x) {
                dis.x = half_size.x - dorg.x;
            }
            else if(dorg.x > dend.x && dend.x < half_size.x) {
                dis.x = dend.x - half_size.x;
            }
            else {
                dis.x = 0;
            }

            if(dorg.y < dend.y && dorg.y < half_size.y) {
                dis.y = half_size.y - dorg.y;
            }
            else if(dorg.y > dend.y && dend.y < half_size.y) {
                dis.y = dend.y - half_size.y;
            }
            else {
                dis.y = 0;
            }

            if(center.z < place0->org.z) {
                dis.z = place0->org.z - center.z + half_size.z;
            }
            else {
                dis.z = place0->end.z - (center.z + half_size.z);
            }

            xyz_dis = fabs(dis.x) + fabs(dis.y) + fabs(dis.z);

            if(xyz_dis < min_dis) {
                min_dis = xyz_dis;
                obj_flg = flg;
                obj_dis = dis;
            }
        }
        else if(flg.x && flg.z) {
            if(dorg.x < dend.x && dorg.x < half_size.x) {
                dis.x = half_size.x - dorg.x;
            }
            else if(dorg.x > dend.x && dend.x < half_size.x) {
                dis.x = dend.x - half_size.x;
            }
            else {
                dis.x = 0;
            }

            if(dorg.z < dend.z && dorg.z < half_size.z) {
                dis.z = half_size.z - dorg.z;
            }
            else if(dorg.z > dend.z && dend.z < half_size.z) {
                dis.z = dend.z - half_size.z;
            }
            else {
                dis.z = 0;
            }

            if(center.y < place0->org.y) {
                dis.y = place0->org.y - center.y + half_size.y;
            }
            else {
                dis.y = place0->end.y - (center.y + half_size.y);
            }

            xyz_dis = fabs(dis.x) + fabs(dis.y) + fabs(dis.z);

            if(xyz_dis < min_dis) {
                min_dis = xyz_dis;
                obj_flg = flg;
                obj_dis = dis;
            }
        }
        else if(flg.y && flg.z) {
            if(dorg.y < dend.y && dorg.y < half_size.y) {
                dis.y = half_size.y - dorg.y;
            }
            else if(dorg.y > dend.y && dend.y < half_size.y) {
                dis.y = dend.y - half_size.y;
            }
            else {
                dis.y = 0;
            }

            if(dorg.z < dend.z && dorg.z < half_size.z) {
                dis.z = half_size.z - dorg.z;
            }
            else if(dorg.z > dend.z && dend.z < half_size.z) {
                dis.z = dend.z - half_size.z;
            }
            else {
                dis.z = 0;
            }

            if(center.x < place0->org.x) {
                dis.x = place0->org.x - center.x + half_size.x;
            }
            else {
                dis.x = place0->end.x - (center.x + half_size.x);
            }

            xyz_dis = fabs(dis.x) + fabs(dis.y) + fabs(dis.z);

            if(xyz_dis < min_dis) {
                min_dis = xyz_dis;
                obj_flg = flg;
                obj_dis = dis;
            }
        }
        else if(flg.x) {
            if(dorg.x < dend.x && dorg.x < half_size.x) {
                dis.x = half_size.x - dorg.x;
            }
            else if(dorg.x > dend.x && dend.x < half_size.x) {
                dis.x = dend.x - half_size.x;
            }
            else {
                dis.x = 0;
            }

            if(center.y < place0->org.y) {
                dis.y = place0->org.y - center.y + half_size.y;
            }
            else {
                dis.y = place0->end.y - (center.y + half_size.y);
            }

            if(center.z < place0->org.z) {
                dis.z = place0->org.z - center.z + half_size.z;
            }
            else {
                dis.z = place0->end.z - (center.z + half_size.z);
            }

            xyz_dis = fabs(dis.x) + fabs(dis.y) + fabs(dis.z);

            if(xyz_dis < min_dis) {
                min_dis = xyz_dis;
                obj_flg = flg;
                obj_dis = dis;
            }
        }
        else if(flg.y) {
            if(dorg.y < dend.y && dorg.y < half_size.y) {
                dis.y = half_size.y - dorg.y;
            }
            else if(dorg.y > dend.y && dend.y < half_size.y) {
                dis.y = dend.y - half_size.y;
            }
            else {
                dis.y = 0;
            }

            if(center.x < place0->org.x) {
                dis.x = place0->org.x - center.x + half_size.x;
            }
            else {
                dis.x = place0->end.x - (center.x + half_size.x);
            }

            if(center.z < place0->org.z) {
                dis.z = place0->org.z - center.z + half_size.z;
            }
            else {
                dis.z = place0->end.z - (center.z + half_size.z);
            }

            xyz_dis = fabs(dis.x) + fabs(dis.y) + fabs(dis.z);

            if(xyz_dis < min_dis) {
                min_dis = xyz_dis;
                obj_flg = flg;
                obj_dis = dis;
            }
        }
        else if(flg.z) {
            if(dorg.z < dend.z && dorg.z < half_size.z) {
                dis.z = half_size.z - dorg.z;
            }
            else if(dorg.z > dend.z && dend.z < half_size.z) {
                dis.z = dend.z - half_size.z;
            }
            else {
                dis.z = 0;
            }

            if(center.x < place0->org.x) {
                dis.x = place0->org.x - center.x + half_size.x;
            }
            else {
                dis.x = place0->end.x - (center.x + half_size.x);
            }

            if(center.y < place0->org.y) {
                dis.y = place0->org.y - center.y + half_size.y;
            }
            else {
                dis.y = place0->end.y - (center.y + half_size.y);
            }

            xyz_dis = fabs(dis.x) + fabs(dis.y) + fabs(dis.z);

            if(xyz_dis < min_dis) {
                min_dis = xyz_dis;
                obj_flg = flg;
                obj_dis = dis;
            }
        }
        else {
            if(center.x < place0->org.x) {
                dis.x = place0->org.x - center.x + half_size.x;
            }
            else {
                dis.x = place0->end.x - (center.x + half_size.x);
            }

            if(center.y < place0->org.y) {
                dis.y = place0->org.y - center.y + half_size.y;
            }
            else {
                dis.y = place0->end.y - (center.y + half_size.y);
            }

            if(center.z < place0->org.z) {
                dis.z = place0->org.z - center.z + half_size.z;
            }
            else {
                dis.z = place0->end.z - (center.z + half_size.z);
            }

            xyz_dis = fabs(dis.x) + fabs(dis.y) + fabs(dis.z);

            if(xyz_dis < min_dis) {
                min_dis = xyz_dis;
                obj_flg = flg;
                obj_dis = dis;
            }
        }
    }

    new_center.x += obj_dis.x;
    new_center.y += obj_dis.y;
    new_center.z += obj_dis.z;

    return new_center;
}

void bin_update(int N) {
    if(STAGE == mGP3D || STAGE == cGP3D || STAGE == mLG3D)
        return bin_update7(N);
    else if(STAGE==cGP2D) 
        return bin_update7_cGP2D(NULL, NULL);
    else if(STAGE == mGP2D)
        return bin_update7_mGP2D();
}

void bin_update7(int N) {
    BIN *bp = NULL;

    for(int i = 0; i < tot_bin_cnt; i++) {
        bp = &bin_mat[i];
        bp->cell_area = 0;
        bp->cell_area2 = 0;
        bp->den = 0;
    }

    for(int i = 0; i < N; i++) {
        if(STAGE == mGP3D) {
            den_comp_3d(i);
        }
        else if(STAGE == cGP3D) {
            if(gcell_st[i].flg != Macro)
                den_comp_3d(i);
        }
    }

    for(int i = 0; i < tot_bin_cnt; i++) {
        bp = &bin_mat[i];
        if(STAGE == mGP3D || STAGE == cGP3D || STAGE == mLG3D) {
            prec area_num2 = bp->cell_area + bp->virt_area + bp->term_area;
            prec area_num = area_num2 + bp->cell_area2;

            bp->den = area_num * inv_bin_area;
            bp->den2 = area_num2 * inv_bin_area;

            copy_den_to_fft_3D(bp->den, bp->p);
        }
    }

    charge_fft_call(1);

    // gsum_cell_phi   = 0;
    // gsum_term_phi   = 0;
    // gsum_virt_phi   = 0;
    gsum_ovfl = 0;

    int i = 0;
    prec ovf_area = 0;

    for(i = 0, bp = &bin_mat[0]; i < tot_bin_cnt; i++, bp++) {
        if(STAGE == mGP3D || STAGE == cGP3D) {
            POS p = bp->p;

            copy_e_from_fft_3D(&(bp->e), p);
            copy_phi_from_fft_3D(&(bp->phi), p);

            // gsum_cell_phi += bp->phi * bp->cell_area + bp->phi *
            // bp->cell_area2;
            // gsum_term_phi += bp->phi * bp->term_area;
            // gsum_virt_phi += bp->phi * bp->virt_area;

            prec ovf_cell_den = max((prec)0.0, bp->den2 - target_cell_den);
            prec ovf_cell_area = ovf_cell_den * bin_area;
            ovf_area += ovf_cell_area;
        }
    }

    // gsum_phi = gsum_cell_phi + gsum_term_phi + gsum_virt_phi;
    gsum_ovfl = ovf_area / total_modu_area;
}

int compare(const void *p1, const void *p2) {
    return (*(prec *)p2 - *(prec *)p1);
}

void bin_update7_mGP2D() {
    BIN *bp = NULL;
    CELLx *cell = NULL;
    TIER *tier = NULL;
    prec area_num = 0;
    prec area_num2 = 0;
    prec ovf_cell_den = 0;
    prec ovf_cell_area = 0;
    prec sum_ovf_area = 0;

    gsum_ovf_area = 0;
    gsum_phi = 0;

    for(int z = 0; z < numLayer; z++) {
        tier = &tier_st[z];

        for(int i = 0; i < tier->tot_bin_cnt; i++) {
            bp = &tier->bin_mat[i];
            bp->cell_area = 0;
            bp->cell_area2 = 0;
        }
        for(int i = 0; i < tier->cell_cnt; i++) {
            cell = tier->cell_st[i];
            den_comp_2d_mGP2D(cell, tier);
        }
        for(int i = 0; i < tier->tot_bin_cnt; i++) {
            bp = &tier->bin_mat[i];
            area_num2 = bp->cell_area + bp->virt_area + bp->term_area;
            area_num = area_num2 + bp->cell_area2;
            bp->den = area_num * tier->inv_bin_area;
            bp->den2 = area_num2 * tier->inv_bin_area;

            copy_den_to_fft_2D(bp->den, bp->p);
        }

        charge_fft_call(0);
        sum_ovf_area = 0;

        for(int i = 0; i < tier->tot_bin_cnt; i++) {
            bp = &tier->bin_mat[i];

            copy_e_from_fft_2D(&(bp->e), bp->p);
            copy_phi_from_fft_2D(&(bp->phi), bp->p);

            gsum_phi += bp->phi * bp->cell_area + bp->phi * bp->cell_area2 +
                        bp->phi * bp->term_area + bp->phi * bp->virt_area;
            // gsum_phi +=
            //  bp->phi * (bp->den - target_cell_den);

            ovf_cell_den = max((prec)0.0, bp->den2 - target_cell_den);
            ovf_cell_area = ovf_cell_den * tier->bin_area;
            sum_ovf_area += ovf_cell_area;
        }
        tier->sum_ovf = sum_ovf_area / tier->modu_area;
        gsum_ovf_area += sum_ovf_area;
    }
    gsum_ovfl = gsum_ovf_area / total_modu_area;
}

// 2D cGP2D
void __bin_update7_cGP2D(cell_t* cells, bin_t* bins, area_t* areas, float** localAr, float** localAr2, size_t* bound, double* time, double* time2) {
    gsum_ovf_area = 0;
    gsum_phi = 0;
    TIER* tier = &tier_st[0];
    const POS dim = tier->dim_bin;
    auto start = std::chrono::steady_clock::now();
    #pragma omp parallel num_threads(numThread) 
    {
        int tid = omp_get_thread_num();
        memset(localAr[tid], 0.0, tier->tot_bin_cnt*sizeof(float));
        memset(localAr2[tid], 0.0, tier->tot_bin_cnt*sizeof(float));
    }
    //#pragma omp parallel num_threads(numThread)
    {
        //int start = ((tid >= 1) ? bound[tid-1] : 0);
        //int end = (tid < numThread-1) ? bound[tid] : tier->cell_cnt;
        #pragma omp parallel for num_threads(numThread)
	for(int p = 0 ; p < tier->cell_cnt; ++p) {
	    int tid = omp_get_thread_num(); 
            cell_t* cell = &cells[p];
            pos_t b0, b1;
            fpos2_t den_pmin = cell->den_pmin, den_pmax = cell->den_pmax;
            b0.x = INT_DOWN((den_pmin.x - tier->bin_org.x) * tier->inv_bin_stp.x);
            b0.y = INT_DOWN((den_pmin.y - tier->bin_org.y) * tier->inv_bin_stp.y);
            b1.x = INT_DOWN((den_pmax.x - tier->bin_org.x) * tier->inv_bin_stp.x);
            b1.y = INT_DOWN((den_pmax.y - tier->bin_org.y) * tier->inv_bin_stp.y);
            if(b0.x < 0) b0.x = 0;
            if(b0.x > tier->dim_bin.x - 1) b0.x = tier->dim_bin.x - 1;
            if(b0.y < 0) b0.y = 0;
            if(b0.y > tier->dim_bin.y - 1) b0.y = tier->dim_bin.y - 1;
            if(b1.x < 0)  b1.x = 0;
            if(b1.x > tier->dim_bin.x - 1) b1.x = tier->dim_bin.x - 1;
            if(b1.y < 0) b1.y = 0;
            if(b1.y > tier->dim_bin.y - 1) b1.y = tier->dim_bin.y - 1;
            cell->b0 = b0;
            cell->b1 = b1;
            int idx = b0.x * tier->dim_bin.y + b0.y;

            int x = 0, y = 0;
            
            bin_t *bpx = NULL, *bpy = NULL;
            for(x = b0.x, bpx = &bins[idx]; x <= b1.x; x++, bpx += tier->dim_bin.y) {
                idx = x * dim.y + b0.y;
                prec max_x = min(bpx->pmax.x, den_pmax.x);
                prec min_x = max(bpx->pmin.x, den_pmin.x);
                for(y = b0.y, bpy = bpx; y <= b1.y; y++, bpy++, idx++) {
                    prec max_y = min(bpy->pmax.y, den_pmax.y);
                    prec min_y = max(bpy->pmin.y, den_pmin.y);
                    prec area_share = (max_x - min_x) * (max_y - min_y) * cell->scale;
                    switch(cell->flg) {
                        case FillerCell: localAr2[tid][idx] += area_share; break;
                        case Macro: localAr[tid][idx] += area_share; break;
                        default: localAr[tid][idx] += area_share; break;                    
                    }
                }
            }
        }
    }
    #pragma omp parallel for num_threads(numThread)
    for(int k = 0; k < tier->tot_bin_cnt; ++k) {
        for(int i = 0; i < numThread; ++i) {
            bins[k].cell_area += localAr[i][k];
        }
        for(int i = 0; i < numThread; ++i) {
            bins[k].cell_area2 += localAr2[i][k];
        }
    }
    *time += time_since(start);
    for(int i = 0; i < tier->tot_bin_cnt; i++) {
        bin_t* bp = &bins[i];
        area_t* a = &areas[i];
        prec area_num2 = bp->cell_area + a->virt_area + a->term_area;
        prec area_num = area_num2 + bp->cell_area2;
        a->den = area_num * tier->inv_bin_area;
        a->den2 = area_num2 * tier->inv_bin_area;
        __copy_den_to_fft_2D__(a->den, a->p);
    }
    charge_fft_call(0);
    prec sum_ovf_area = 0;
    for(int i = 0; i < tier->tot_bin_cnt; i++) {
        bin_t* bp = &bins[i];
        area_t* a = &areas[i];    
        __copy_e_from_fft_2D__(&(bp->e), a->p);
        __copy_phi_from_fft_2D__(&(bp->phi), a->p);
        gsum_phi += bp->phi * (bp->cell_area + bp->cell_area2 + a->term_area + a->virt_area);
        sum_ovf_area += max((prec)0.0, a->den2 - target_cell_den) * tier->bin_area;
    }
    
    tier->sum_ovf = sum_ovf_area / tier->modu_area;
    gsum_ovf_area += sum_ovf_area;
    gsum_ovfl = gsum_ovf_area / total_modu_area;
    *time2 += time_since(start);
}


// 2D cGP2D

void CellBox(Cell_t* cells) {
    TIER* tier = &tier_st[0];
    FPOS bin_org = tier->bin_org, stp = tier->inv_bin_stp;
    POS dim = tier->dim_bin;
    int N = tier->cell_cnt;
    for(int k = 0 ; k < N; ++k) {
        BNDBDY(0, x, min);
    }
    for(int k = 0; k < N; ++k) {
        BNDBDY(0, y, min);
    }
    for(int k = 0; k < N; ++k) {
        BNDBDY(1, x, max);
    }
    for(int k = 0; k < N; ++k) {
        BNDBDY(1, y, max);
    }
}

void Density(Cell_t* cells, Bin_t* bins, int* bound, float** localAr, float** localAr2) {
    const TIER* tier = &tier_st[0];
    const FPOS bin_org = tier->bin_org, stp = tier->inv_bin_stp;
    const POS dim = tier->dim_bin;
    auto start = std::chrono::system_clock::now();
    #pragma omp parallel num_threads(numThread) 
    {
        int tid = omp_get_thread_num();
        memset(localAr[tid], 0.0, tier->tot_bin_cnt*sizeof(float));
        memset(localAr2[tid], 0.0, tier->tot_bin_cnt*sizeof(float));
    }
    //TODO cast to (void*)
    auto time0 = std::chrono::system_clock::now();
	#pragma omp parallel num_threads(numThread)
    {
        int tid = omp_get_thread_num();
        int start = ((tid >= 1) ? bound[tid-1] : 0);
        int end = (tid < numThread-1) ? bound[tid] : tier->cell_cnt;
        for(int p = start ; p < end; ++p) {
            const pos_t b0 = {cells->b0_x[p],cells->b0_y[p]}, b1 = {cells->b1_x[p], cells->b1_y[p]};
            float scalar = cells->scale[p];
            char flg = cells->flg[p];
            int idx = b0.x * dim.y + b0.y;
            for(int x = b0.x; x <= b1.x; x++) {
                idx = x * dim.y + b0.y;
                const prec max_x = min(bins->pmax_x[idx], cells->den_pmax_x[p]);
                const prec min_x = max(bins->pmin_x[idx], cells->den_pmin_x[p]);
                prec dx = (max_x - min_x);
                for(int y = b0.y; y <= b1.y; y++, idx++) {
                    prec max_y = min(bins->pmax_y[idx], cells->den_pmax_y[p]);
                    prec min_y = max(bins->pmin_y[idx], cells->den_pmin_y[p]);
                    prec area_share = dx * (max_y - min_y) * scalar;
                    switch(flg) {
                        case FillerCell: localAr2[tid][idx] += area_share; break;
                        case Macro: localAr[tid][idx] += area_share; break;
                        default: localAr[tid][idx] += area_share; break;                    
                    }
                }
            }
        }
    }
    auto time1 = std::chrono::system_clock::now();
    #pragma omp parallel for num_threads(numThread)
    for(int k = 0; k < tier->tot_bin_cnt; ++k) {
        for(int i = 0; i < numThread; ++i) {
            bins->cell_area[k] += localAr[i][k];
        }
        for(int i = 0; i < numThread; ++i) {
            bins->cell_area2[k] += localAr2[i][k];
        }
    }
}

void FFTSolve(bin_t* bins, area_t* areas) {
    gsum_ovf_area = 0;
    gsum_phi = 0;
    TIER* tier = &tier_st[0];
    for(int i = 0; i < tier->tot_bin_cnt; i++) {
        bin_t* bp = &bins[i];
        area_t* a = &areas[i];
        prec area_num2 = bp->cell_area + a->virt_area + a->term_area;
        prec area_num = area_num2 + bp->cell_area2;
        a->den = area_num * tier->inv_bin_area;
        a->den2 = area_num2 * tier->inv_bin_area;
        __copy_den_to_fft_2D__(a->den, a->p);
    }
    charge_fft_call(0);
    prec sum_ovf_area = 0;
    for(int i = 0; i < tier->tot_bin_cnt; i++) {
        bin_t* bp = &bins[i];
        area_t* a = &areas[i];    
        __copy_e_from_fft_2D__(&(bp->e), a->p);
        __copy_phi_from_fft_2D__(&(bp->phi), a->p);
        gsum_phi += bp->phi * (bp->cell_area + bp->cell_area2 + a->term_area + a->virt_area);
        sum_ovf_area += max((prec)0.0, a->den2 - target_cell_den) * tier->bin_area;
    }
    
    tier->sum_ovf = sum_ovf_area / tier->modu_area;
    gsum_ovf_area += sum_ovf_area;
    gsum_ovfl = gsum_ovf_area / total_modu_area;
}



// 2D cGP2D
void bin_update7_cGP2D(double* time0, double* time1) {
    auto start = std::chrono::steady_clock::now();
    BIN *bp = NULL;

    gsum_ovf_area = 0;
    gsum_phi = 0;

    TIER* tier = &tier_st[0];
    bool timeon = false;
    double time = 0.0f;

    omp_set_num_threads(numThread);
    int i = 0;

    if( timeon ) {time_start(&time);}

    if(timeon) { time_end(&time); cout << "initialize: " << time << endl; time_start(&time);}

    // update cell_area & cell_area2
//    omp_set_num_threads(1);

//    CELL* cellInst = tier->cell_st;
//    cout << tier->cell_cnt << endl;
//    exit(0);
//#pragma omp parallel \
    default(none) \
    shared(tier) \
    private(i) 
//    {

//#pragma omp for


            for(int i = 0; i < tier->cell_cnt; i++) {
                CELLx* cell = tier->cell_st[i];
                POS b0, b1;
                b0.x = INT_DOWN((cell->den_pmin.x - tier->bin_org.x) * tier->inv_bin_stp.x);
                b0.y = INT_DOWN((cell->den_pmin.y - tier->bin_org.y) * tier->inv_bin_stp.y);
                b1.x = INT_DOWN((cell->den_pmax.x - tier->bin_org.x) * tier->inv_bin_stp.x);
                b1.y = INT_DOWN((cell->den_pmax.y - tier->bin_org.y) * tier->inv_bin_stp.y);
                if(b0.x < 0) b0.x = 0;
                if(b0.x > tier->dim_bin.x - 1) b0.x = tier->dim_bin.x - 1;
                if(b0.y < 0) b0.y = 0;
                if(b0.y > tier->dim_bin.y - 1) b0.y = tier->dim_bin.y - 1;
                if(b1.x < 0)  b1.x = 0;
                if(b1.x > tier->dim_bin.x - 1)  b1.x = tier->dim_bin.x - 1;
                if(b1.y < 0) b1.y = 0;
                if(b1.y > tier->dim_bin.y - 1) b1.y = tier->dim_bin.y - 1;

                int idx = b0.x * tier->dim_bin.y + b0.y;

                int x = 0, y = 0;
                BIN *bpx = NULL, *bpy = NULL;
                for(x = b0.x, bpx = &tier->bin_mat[idx]; x <= b1.x; x++, bpx += tier->dim_bin.y) {
                    prec max_x = min(bpx->pmax.x, cell->den_pmax.x);
                    prec min_x = max(bpx->pmin.x, cell->den_pmin.x);
                    for(y = b0.y, bpy = bpx; y <= b1.y; y++, bpy++) {
                        prec max_y = min(bpy->pmax.y, cell->den_pmax.y);
                        prec min_y = max(bpy->pmin.y, cell->den_pmin.y);
                        prec area_share = (max_x - min_x) * (max_y - min_y) * cell->den_scal;
                        switch(cell->flg) {
                            case FillerCell: bpy->cell_area2 += area_share; break;
                            case Macro: bpy->cell_area += area_share; break;
                            default: bpy->cell_area += area_share; break;
                        }
                    }
                }
            }
//}
if(time0)
    *time0 += time_since(start);
    if(timeon) { time_end(&time); cout << "fill cell_area: " << time << endl; time_start(&time);}

    omp_set_num_threads(numThread);
#pragma omp parallel \
    default(none) \
    shared(tier) \
    private(i) 
    {

#pragma omp for
    for(i = 0; i < tier->tot_bin_cnt; i++) {
        BIN* bp = &tier->bin_mat[i];

        prec area_num2 = bp->cell_area + bp->virt_area + bp->term_area;
        prec area_num = area_num2 + bp->cell_area2;

        bp->den = area_num * tier->inv_bin_area;
        bp->den2 = area_num2 * tier->inv_bin_area;

        copy_den_to_fft_2D(bp->den, bp->p);
    }
    }
    if(timeon) { time_end(&time); cout << "bin 1st update: " << time << endl; time_start(&time); }

    charge_fft_call(0);
    if(timeon) { time_end(&time); cout << "charge_fft_call: " << time << endl; time_start(&time); }

    prec sum_ovf_area = 0;
//#pragma omp parallel \
    default(none) \
    shared(tier, sum_ovf_area, gsum_phi, target_cell_den) \
    private(i) 
//    {
//        BIN *bp = NULL;
//#pragma omp for
    for(i = 0; i < tier->tot_bin_cnt; i++) {
        bp = &tier->bin_mat[i];

        copy_e_from_fft_2D(&(bp->e), bp->p);
        copy_phi_from_fft_2D(&(bp->phi), bp->p);

        gsum_phi += bp->phi * bp->cell_area + bp->phi * bp->cell_area2 +
            bp->phi * bp->term_area + bp->phi * bp->virt_area;
        
//        cout << i << " " << gsum_phi << endl;

        // gsum_phi +=
        //  bp->phi * (bp->den - target_cell_den);

//        prec ovf_cell_den = max((prec)0.0, bp->den2 - target_cell_den);
//        prec ovf_cell_area = ovf_cell_den * tier->bin_area;
//        sum_ovf_area += ovf_cell_area;
        
        sum_ovf_area += max((prec)0.0, bp->den2 - target_cell_den) * tier->bin_area;
    }

    if(timeon) { time_end(&time); cout << "bin final loop: " << time << endl; }

//    cout << "sum_ovf_area: " << sum_ovf_area << endl;

    tier->sum_ovf = sum_ovf_area / tier->modu_area;
    gsum_ovf_area += sum_ovf_area;
//    cout << "gsum_ovf_area: " << gsum_ovf_area << endl;
//    cout << "total_modu_area: " << total_modu_area << endl;

    gsum_ovfl = gsum_ovf_area / total_modu_area;
//    cout << "gsumovfl: " << gsum_ovfl << endl;
//    exit(1);
if(time1)
*time1 += time_since(start);
}

void get_term_den(prec *den) {
    int i = 0;
    prec sum_den = 0, term_den = 0, virt_den = 0;
    BIN *bp = NULL;

    for(i = 0; i < tot_bin_cnt; i++) {
        bp = &bin_mat[i];
        term_den = bp->term_area * bp->den;
        virt_den = bp->virt_area * bp->den;
        sum_den += term_den + virt_den;
    }

    *den = sum_den;
}

prec get_den2(prec area_num, prec area_denom) {
    return area_denom < MIN_AREA3 ? area_num / MIN_AREA3
                                  : area_num / area_denom;
}

int is_IO_block(TERM *term) {
    if(point_in_rect(term->center, place.org, place.end))
        return 0;
    else
        return 1;
}

prec get_bins_mac(FPOS center, MODULE *mac) {
    int x = 0, y = 0, z = 0, idx = 0;
    prec area_share = 0;
    prec min_x = 0, min_y = 0, min_z = 0;
    prec max_x = 0, max_y = 0, max_z = 0;
    POS b0 = zeroPoint, b1 = zeroPoint;
    FPOS pmin = zeroFPoint, pmax = zeroFPoint;
    BIN *bpx = NULL, *bpy = NULL, *bpz = NULL;
    prec cost = 0.0;

    pmin.x = center.x - mac->half_size.x;
    pmin.y = center.y - mac->half_size.y;
    pmin.z = center.z - mac->half_size.z;

    pmax.x = center.x + mac->half_size.x;
    pmax.y = center.y + mac->half_size.y;
    pmax.z = center.z + mac->half_size.z;

    b0.x = INT_DOWN ((pmin.x - bin_org.x) * inv_bin_stp.x);
    b0.y = INT_DOWN ((pmin.y - bin_org.y) * inv_bin_stp.y);
    b0.z = INT_DOWN ((pmin.z - bin_org.z) * inv_bin_stp.z);

    b1.x = INT_DOWN((pmax.x - bin_org.x) * inv_bin_stp.x);
    b1.y = INT_DOWN((pmax.y - bin_org.y) * inv_bin_stp.y);
    b1.z = INT_DOWN((pmax.z - bin_org.z) * inv_bin_stp.z);

    legal_bin_idx(&b0);
    legal_bin_idx(&b1);

    idx = b0.x * msh_yz + b0.y * msh.z + b0.z;

    for(x = b0.x, bpx = &bin_mat[idx]; x <= b1.x; x++, bpx += msh_yz) {
        // bp = & bin_mat[x_idx];

        max_x = min(bpx->pmax.x, pmax.x);
        min_x = max(bpx->pmin.x, pmin.x);

        for(y = b0.y, bpy = bpx; y <= b1.y; y++, bpy += msh.z) {
            max_y = min(bpy->pmax.y, pmax.y);
            min_y = max(bpy->pmin.y, pmin.y);

            for(z = b0.z, bpz = bpy; z <= b1.z; z++, bpz++) {
                max_z = min(bpz->pmax.z, pmax.z);
                min_z = max(bpz->pmin.z, pmin.z);

                area_share =
                    (max_x - min_x) * (max_y - min_y) * (max_z - min_z);
                // cell->den_scal;

                cost += area_share * bpz->no_mac_den;
            }
        }
    }

    return cost;
}

void get_bins(FPOS center, CELLx *cell, POS *st, prec *share_st, int *bin_cnt) {
    prec *x_st = bin_share_x_st;
    prec *y_st = bin_share_y_st;
    prec *z_st = bin_share_z_st;
    int x = 0, y = 0, z = 0;
    int bin_cnt0 = 0;
    prec sum_area = 0, diff_area = 0;
    prec x0 = 0, y0 = 0, z0 = 0;
    prec x1 = 0, y1 = 0, z1 = 0;
    prec min_x = 0, max_x = 0;
    prec min_y = 0, max_y = 0;
    prec min_z = 0, max_z = 0;
    prec area_share = 0;
    POS b0 = zeroPoint, b1 = zeroPoint;
    FPOS share = zeroFPoint;
    FPOS pmin = zeroFPoint, pmax = zeroFPoint;

    pmin.x = center.x - cell->half_den_size.x;
    pmin.y = center.y - cell->half_den_size.y;
    pmin.z = center.z - cell->half_den_size.z;

    pmax.x = center.x + cell->half_den_size.x;
    pmax.y = center.y + cell->half_den_size.y;
    pmax.z = center.z + cell->half_den_size.z;

    b0.x = INT_DOWN((pmin.x - bin_org.x) * inv_bin_stp.x);
    b0.y = INT_DOWN((pmin.y - bin_org.y) * inv_bin_stp.y);
    b0.z = INT_DOWN((pmin.z - bin_org.z) * inv_bin_stp.z);

    b1.x = INT_DOWN((pmax.x - bin_org.x) * inv_bin_stp.x);
    b1.y = INT_DOWN((pmax.y - bin_org.y) * inv_bin_stp.y);
    b1.z = INT_DOWN((pmax.z - bin_org.z) * inv_bin_stp.z);

    legal_bin_idx(&b0);
    legal_bin_idx(&b1);

    for(x = b0.x; x <= b1.x; x++) {
        x0 = (prec)x * bin_stp.x + bin_org.x;

        min_x = max(x0, pmin.x);

        x1 = (prec)(x + 1) * bin_stp.x + bin_org.x;

        max_x = min(x1, pmax.x);

        share.x = max_x - min_x;

        x_st[x - b0.x] = share.x;
    }

    for(y = b0.y; y <= b1.y; y++) {
        y0 = (prec)y * bin_stp.y + bin_org.y;

        min_y = max(y0, pmin.y);

        y1 = (prec)(y + 1) * bin_stp.y + bin_org.y;

        max_y = min(y1, pmax.y);

        share.y = max_y - min_y;

        y_st[y - b0.y] = share.y;
    }

    if(flg_3dic) {
        for(z = b0.z; z <= b1.z; z++) {
            z0 = (prec)z * bin_stp.z + bin_org.z;

            min_z = max(z0, pmin.z);

            z1 = (prec)(z + 1) * bin_stp.z + bin_org.z;

            max_z = min(z1, pmax.z);

            share.z = max_z - min_z;

            z_st[z - b0.z] = share.z;
        }
    }
    else {
        z_st[0] = 1.0;
        b0.z = b1.z = 0;
    }

    for(x = b0.x; x <= b1.x; x++) {
        for(y = b0.y; y <= b1.y; y++) {
            for(z = b0.z; z <= b1.z; z++) {
                share.x = x_st[x - b0.x];
                share.y = y_st[y - b0.y];
                share.z = z_st[z - b0.z];

                area_share = share.x * share.y * share.z * cell->den_scal;

                sum_area += area_share;

                st[bin_cnt0].x = x;
                st[bin_cnt0].y = y;
                st[bin_cnt0].z = z;

                share_st[bin_cnt0] = area_share;

                bin_cnt0++;
            }
        }
    }

    *bin_cnt = bin_cnt0;

    diff_area = fabs(sum_area - cell->area);

    if(diff_area > 0.01 * cell->area) {
        g_rrr++;
    }

    return;
}

prec get_ovfl() {
    prec ovfl1 = 0;
    prec ovf = 0;
    BIN *bp = NULL;

    for(int i = 0; i < tot_bin_cnt; i++) {
        bp = &bin_mat[i];
        ovf =
            max((prec)0.0,
                /* bp->den */ bp->den2 - target_cell_den);  // total_cell_den);
        ovf *= bin_area;

        ovfl1 += ovf;
    }
    return ovfl1 / total_modu_area;
}

bool get_common_rect(FPOS pmin1, FPOS pmax1, FPOS pmin2, FPOS pmax2, FPOS *p1,
                     FPOS *p2) {
    if((pmax1.x < pmin2.x || pmax1.x == pmin2.x) ||
       (pmin1.x > pmax2.x || pmin1.x == pmax2.x) ||
       (pmax1.y < pmin2.y || pmax1.y == pmin2.y) ||
       (pmin1.y > pmax2.y || pmin1.y == pmax2.y)) {
        p1->Set(-1);
        p2->Set(-1);
        return false;
    }
    else if(flg_3dic && (pmax1.z <= pmin2.z || pmin1.z >= pmax2.z)) {
        p1->Set(-1);
        p2->Set(-1);
        return false;
    }
    else {
        if((pmax1.x > pmin2.x || pmax1.x == pmin2.x) &&
           (pmax1.x < pmax2.x || pmax1.x == pmax2.x)) {
            if(pmin1.x > pmin2.x || pmin1.x == pmin2.x) {
                p1->x = pmin1.x;
                p2->x = pmax1.x;
            }
            else {
                p1->x = pmin2.x;
                p2->x = pmax1.x;
            }
        }
        else if(pmax1.x > pmax2.x || pmax1.x == pmax2.x) {
            if(pmin1.x > pmin2.x || pmin1.x == pmin2.x) {
                p1->x = pmin1.x;
                p2->x = pmax2.x;
            }
            else {
                p1->x = pmin2.x;
                p2->x = pmax2.x;
            }
        }
        else {
            error(0, 0, "get_common_area: eff_x error");
        }

        if((pmax1.y > pmin2.y || pmax1.y == pmin2.y) &&
           (pmax1.y < pmax2.y || pmax1.y == pmax2.y)) {
            if(pmin1.y > pmin2.y || pmin1.y == pmin2.y) {
                p1->y = pmin1.y;
                p2->y = pmax1.y;
            }
            else {
                p1->y = pmin2.y;
                p2->y = pmax1.y;
            }
        }
        else if(pmax1.y > pmax2.y || pmax1.y == pmax2.y) {
            if(pmin1.y > pmin2.y || pmin1.y == pmin2.y) {
                p1->y = pmin1.y;
                p2->y = pmax2.y;
            }
            else {
                p1->y = pmin2.y;
                p2->y = pmax2.y;
            }
        }
        else {
            error(0, 0, "get_common_area: eff_x error");
        }

        if(flg_3dic) {
            if((pmax1.z > pmin2.z || pmax1.z == pmin2.z) &&
               (pmax1.z < pmax2.z || pmax1.z == pmax2.z)) {
                if(pmin1.z > pmin2.z || pmin1.z == pmin2.z) {
                    p1->z = pmin1.z;
                    p2->z = pmax1.z;
                }
                else {
                    p1->z = pmin2.z;
                    p2->z = pmax1.z;
                }
            }
            else if(pmax1.z > pmax2.z || pmax1.z == pmax2.z) {
                if(pmin1.z > pmin2.z || pmin1.z == pmin2.z) {
                    p1->z = pmin1.z;
                    p2->z = pmax2.z;
                }
                else {
                    p1->z = pmin2.z;
                    p2->z = pmax2.z;
                }
            }
            else {
                error(0, 0, "get_common_area: eff_x error");
            }
        }
    }
    return true;
}


prec get_common_area(FPOS pmin1, FPOS pmax1, FPOS pmin2, FPOS pmax2) {
    prec eff_x = 0;
    prec eff_y = 0;
    prec eff_z = 0;

    if(pmax1.x <= pmin2.x || pmin1.x >= pmax2.x || pmax1.y <= pmin2.y ||
       pmin1.y >= pmax2.y)
        return 0;
    else if(flg_3dic && (pmin1.z >= pmax2.z || pmax1.z <= pmin2.z))
        return 0;
    else {
        if((pmax1.x > pmin2.x || pmax1.x == pmin2.x) &&
           (pmax1.x < pmax2.x || pmax1.x == pmax2.x)) {
            if(pmin1.x > pmin2.x || pmin1.x == pmin2.x)
                eff_x = pmax1.x - pmin1.x;
            else
                eff_x = pmax1.x - pmin2.x;
        }
        else if(pmax1.x > pmax2.x || pmax1.x == pmax2.x) {
            if(pmin1.x > pmin2.x || pmin1.x == pmin2.x)
                eff_x = pmax2.x - pmin1.x;
            else
                eff_x = pmax2.x - pmin2.x;
        }
        else {
            error(0, 0, "get_common_area: eff_x error");
            g_rrr++;
        }

        if((pmax1.y > pmin2.y || pmax1.y == pmin2.y) &&
           (pmax1.y < pmax2.y || pmax1.y == pmax2.y)) {
            if(pmin1.y > pmin2.y || pmin1.y == pmin2.y)
                eff_y = pmax1.y - pmin1.y;
            else
                eff_y = pmax1.y - pmin2.y;
        }
        else if(pmax1.y > pmax2.y || pmax1.y == pmax2.y) {
            if(pmin1.y > pmin2.y || pmin1.y == pmin2.y)
                eff_y = pmax2.y - pmin1.y;
            else
                eff_y = pmax2.y - pmin2.y;
        }
        else {
            error(0, 0, "get_common_area: eff_x error");
            g_rrr++;
        }

        if(!flg_3dic)
            return eff_x * eff_y;
        else {
            if((pmax1.z > pmin2.z || pmax1.z == pmin2.z) &&
               (pmax1.z < pmax2.z || pmax1.z == pmax2.z)) {
                if(pmin1.z > pmin2.z || pmin1.z == pmin2.z)
                    eff_z = pmax1.z - pmin1.z;
                else
                    eff_z = pmax1.z - pmin2.z;
            }
            else if(pmax1.z > pmax2.z || pmax1.z == pmax2.z) {
                if(pmin1.z > pmin2.z || pmin1.z == pmin2.z)
                    eff_z = pmax2.z - pmin1.z;
                else
                    eff_z = pmax2.z - pmin2.z;
            }
            else {
                error(0, 0, "get_common_area: eff_x error");
                g_rrr++;
            }
            return eff_x * eff_y * eff_z;
        }
    }
}

int point_in_rect(FPOS p0, FPOS pmin, FPOS pmax) {
    return p0.x >= pmin.x - Epsilon && p0.x <= pmax.x + Epsilon &&
           p0.y >= pmin.y - Epsilon && p0.y <= pmax.y + Epsilon;
}

/* it never used
void bin_den_update_for_filler_adj(void) {
    bin_clr();
    bin_cell_area_update();
}
*/

void den_update() {
    bin_clr();
    bin_cell_area_update();
    bin_den_mac_update();
}

void bin_den_mac_update(void) {
    for(int i = 0; i < tot_bin_cnt; i++) {
        BIN* bp = &bin_mat[i];
        prec area_num =
            bp->cell_area + bp->macro_area + bp->virt_area + bp->term_area;
        prec area_num2 = bp->cell_area + bp->virt_area + bp->term_area;
        prec area_denom = bin_area;
        bp->den = get_den2(area_num, area_denom);
        bp->no_mac_den = get_den2(area_num2, area_denom);
    }
}

// this function seems to have problem : calling bin_share_st
void bin_cell_area_update(void) {
    int i = 0, j = 0;
    prec area_share = 0;
    CELLx *cell = NULL;
    FPOS center = zeroFPoint;
    FPOS half_size = zeroFPoint;
    POS p = zeroPoint;

    for(i = 0; i < moduleCNT; i++) {
        center = moduleInstance[i].center;

        cell = &gcell_st[i];

        half_size = cell->half_size;

        get_bins(center, cell, bin_st, bin_share_st, &bin_cnt);

        for(j = 0; j < bin_cnt; j++) {
            p = bin_st[j];

            BIN* bp = get_bin_from_idx(p);

            area_share = bin_share_st[j];

            if(cell->flg == Macro) {
                bp->macro_area2 += area_share;
                bp->macro_area +=
                    area_share * global_macro_area_scale;  // mov-macro area
            }
            else if(cell->flg == StdCell)
                bp->cell_area += area_share;  // std area
            else if(cell->flg == FillerCell)
                bp->filler_area += area_share;  // filler area
            else {
                printf("error cell no group ID\n");
                exit(1);
            }
        }
    }
}

void bin_clr(void) {
    for(int i = 0; i < tot_bin_cnt; i++) {
        BIN* bp = &bin_mat[i];
        bp->cell_area = 0;
        bp->cell_area2 = 0;
        bp->filler_area = 0;
        bp->macro_area = 0;
        bp->macro_area2 = 0;
        bp->den = 0;
    }
}

prec get_all_macro_den(void) {
    MODULE *mac = NULL;
    MODULE *mdp = NULL;
    FPOS center = zeroFPoint;
    FPOS half_size = zeroFPoint;
    prec den_cost = 0;

    for(int i = 0; i < moduleCNT; i++) {
        mdp = &moduleInstance[i];

        if(mdp->flg != Macro)
            continue;

        center = mdp->center;
        mac = macro_st[mdp->mac_idx];
        half_size = mac->half_size;
        den_cost += get_bins_mac(center, mac);
    }
    return den_cost;
}

prec get_mac_den(int idx) {
    MODULE *mac = macro_st[idx];
    FPOS center = mac->center;
    prec den_cost = 0;

    den_cost = get_bins_mac(center, mac);
    return den_cost;
}

int idx_in_bin_rect(POS *p, POS pmin, POS pmax) {
    if(p->x < 0) {
        *p = pmin;
        return 1;
    }
    else {
        if(flg_3dic) {
            if(p->z < pmax.z) {
                p->z++;
                return 1;
            }
            else if(p->y < pmax.y) {
                p->z = pmin.z;
                p->y++;
                return 1;
            }
            else if(p->x < pmax.x) {
                p->z = pmin.z;
                p->y = pmin.y;
                p->x++;
                return 1;
            }
            else
                return 0;
        }
        else {
            if(p->y < pmax.y) {
                p->y++;
                return 1;
            }
            else if(p->x < pmax.x) {
                p->y = pmin.y;
                p->x++;
                return 1;
            }
            else
                return 0;
        }
    }
}

void den_comp_2d_mGP2D(CELLx *cell, TIER *tier) {
    int x = 0, y = 0;

    /*  CELLx*cell=&gcell_st[cell_idx]; */
    /*  FPOS pmin,pmax; */
    int idx = 0;
    prec area_share = 0;
    prec min_x = 0, min_y = 0;
    prec max_x = 0, max_y = 0;
    BIN *bpx = NULL, *bpy = NULL /* ,*bpz=NULL */;
    
    POS b0, b1;

    b0.x = INT_DOWN((cell->den_pmin.x - tier->bin_org.x) * tier->inv_bin_stp.x);
    b0.y = INT_DOWN((cell->den_pmin.y - tier->bin_org.y) * tier->inv_bin_stp.y);

    b1.x = INT_DOWN((cell->den_pmax.x - tier->bin_org.x) * tier->inv_bin_stp.x);
    b1.y = INT_DOWN((cell->den_pmax.y - tier->bin_org.y) * tier->inv_bin_stp.y);

    if(b0.x < 0)
        b0.x = 0;
    if(b0.x > tier->dim_bin.x - 1)
        b0.x = tier->dim_bin.x - 1;
    if(b0.y < 0)
        b0.y = 0;
    if(b0.y > tier->dim_bin.y - 1)
        b0.y = tier->dim_bin.y - 1;

    if(b1.x < 0)
        b1.x = 0;
    if(b1.x > tier->dim_bin.x - 1)
        b1.x = tier->dim_bin.x - 1;
    if(b1.y < 0)
        b1.y = 0;
    if(b1.y > tier->dim_bin.y - 1)
        b1.y = tier->dim_bin.y - 1;

    idx = b0.x * tier->dim_bin.y + b0.y;

    for(x = b0.x, bpx = &tier->bin_mat[idx]; x <= b1.x;
        x++, bpx += tier->dim_bin.y) {
        max_x = min(bpx->pmax.x, cell->den_pmax.x);
        min_x = max(bpx->pmin.x, cell->den_pmin.x);

        for(y = b0.y, bpy = bpx; y <= b1.y; y++, bpy++) {
            max_y = min(bpy->pmax.y, cell->den_pmax.y);
            min_y = max(bpy->pmin.y, cell->den_pmin.y);

            area_share = (max_x - min_x) * (max_y - min_y) * cell->size.z *
                         cell->den_scal;

            if(cell->flg == FillerCell) {
                bpy->cell_area2 += area_share;
            }
            else if(cell->flg == Macro) {
                bpy->cell_area += area_share * global_macro_area_scale;
            }
            else {
                bpy->cell_area += area_share;
            }
        }
    }
}

// calculate 
//
// tier->bin_mat->cell_area (Normal & Macro) and 
// tier->bin_mat->cell_area2(Filler Cell)
//
void den_comp_2d_cGP2D(CELLx *cell, TIER *tier) {

    POS b0, b1;
    b0.x = INT_DOWN((cell->den_pmin.x - tier->bin_org.x) * tier->inv_bin_stp.x);
    b0.y = INT_DOWN((cell->den_pmin.y - tier->bin_org.y) * tier->inv_bin_stp.y);
    b1.x = INT_DOWN((cell->den_pmax.x - tier->bin_org.x) * tier->inv_bin_stp.x);
    b1.y = INT_DOWN((cell->den_pmax.y - tier->bin_org.y) * tier->inv_bin_stp.y);
    //b0.Legalize(zeroFPoint, tier->dim_bin);
    //b1.Legalize(zeroFPoint, tier->dim_bin);

    int idx = b0.x * tier->dim_bin.y + b0.y;

    int x = 0, y = 0;
    BIN *bpx = NULL, *bpy = NULL;
    for(x = b0.x, bpx = &tier->bin_mat[idx]; x <= b1.x; x++, bpx += tier->dim_bin.y) {
        prec max_x = min(bpx->pmax.x, cell->den_pmax.x);
        prec min_x = max(bpx->pmin.x, cell->den_pmin.x);
        for(y = b0.y, bpy = bpx; y <= b1.y; y++, bpy++) {
            prec max_y = min(bpy->pmax.y, cell->den_pmax.y);
            prec min_y = max(bpy->pmin.y, cell->den_pmin.y);
            prec area_share = (max_x - min_x) * (max_y - min_y) * cell->den_scal;
            switch(cell->flg) {
                case FillerCell: bpy->cell_area2 += area_share; break;
                case Macro: bpy->cell_area += area_share; break;
                default: bpy->cell_area += area_share; break;
            }
        }
    }
}

void den_comp_3d(int cell_idx) {
    int x = 0, y = 0, z = 0;
    POS b0, b1;
    CELLx *cell = &gcell_st[cell_idx];
    prec area_share = 0;
    int idx = 0;
    prec min_x = 0, min_y = 0, min_z = 0;
    prec max_x = 0, max_y = 0, max_z = 0;
    BIN *bpx = NULL, *bpy = NULL, *bpz = NULL;

    b0.x = INT_DOWN((cell->den_pmin.x - bin_org.x) * inv_bin_stp.x);
    b0.y = INT_DOWN((cell->den_pmin.y - bin_org.y) * inv_bin_stp.y);
    b0.z = INT_DOWN((cell->den_pmin.z - bin_org.z) * inv_bin_stp.z);

    b1.x = INT_DOWN((cell->den_pmax.x - bin_org.x) * inv_bin_stp.x);
    b1.y = INT_DOWN((cell->den_pmax.y - bin_org.y) * inv_bin_stp.y);
    b1.z = INT_DOWN((cell->den_pmax.z - bin_org.z) * inv_bin_stp.z);

    if(b0.x < 0)
        b0.x = 0;
    if(b0.y < 0)
        b0.y = 0;
    if(b0.z < 0)
        b0.z = 0;

    if(b0.x >= max_bin.x)
        b0.x = max_bin.x - 1;
    if(b0.y >= max_bin.y)
        b0.y = max_bin.y - 1;
    if(b0.z >= max_bin.z)
        b0.z = max_bin.z - 1;

    if(b1.x < 0)
        b1.x = 0;
    if(b1.y < 0)
        b1.y = 0;
    if(b1.z < 0)
        b1.z = 0;

    if(b1.x >= max_bin.x)
        b1.x = max_bin.x - 1;
    if(b1.y >= max_bin.y)
        b1.y = max_bin.y - 1;
    if(b1.z >= max_bin.z)
        b1.z = max_bin.z - 1;

    idx = b0.x * msh_yz + b0.y * msh.z + b0.z;

    for(x = b0.x, bpx = &bin_mat[idx]; x <= b1.x; x++, bpx += msh_yz) {
        max_x = min(bpx->pmax.x, cell->den_pmax.x);
        min_x = max(bpx->pmin.x, cell->den_pmin.x);

        for(y = b0.y, bpy = bpx; y <= b1.y; y++, bpy += msh.z) {
            max_y = min(bpy->pmax.y, cell->den_pmax.y);
            min_y = max(bpy->pmin.y, cell->den_pmin.y);

            for(z = b0.z, bpz = bpy; z <= b1.z; z++, bpz++) {
                max_z = min(bpz->pmax.z, cell->den_pmax.z);
                min_z = max(bpz->pmin.z, cell->den_pmin.z);

                area_share = (max_x - min_x) * (max_y - min_y) *
                             (max_z - min_z) * cell->den_scal;

                if(cell_idx < moduleCNT) {
                    if(cell->flg == Macro)
                        bpz->cell_area += area_share * global_macro_area_scale;
                    else
                        bpz->cell_area += area_share;
                }
                else {
                    bpz->cell_area2 += area_share;
                }
            }
        }
    }

    return;
}

void bin_delete_mGP2D(void) {
    for(int currTier = 0; currTier < numLayer; currTier++) {
        TIER* tier = &tier_st[currTier];
        tier->bin_mat = NULL;
        tier->tot_bin_cnt = 0;
        // mkl_free (bin_mat_st[currTier]);
    }
    mkl_free(bin_mat_st);
    bin_mat_st = NULL;
}
