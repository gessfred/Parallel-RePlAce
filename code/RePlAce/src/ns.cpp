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
#include <sys/time.h>
#include <algorithm>
#include <cfloat>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <omp.h>
// #include        <ncurses.h>

#include "bookShelfIO.h"
#include "bin.h"
#include "charge.h"
#include "global.h"
#include "mkl.h"
#include "ns.h"
#include "opt.h"
#include "plot.h"
#include "wlen.h"
#include "parallel.hpp"
static int backtrack_cnt = 0;

void myNesterov::nesterov_opt() {
    auto nsS = std::chrono::steady_clock::now();
    int last_iter = 0;
    InitializationCommonVar();
    InitializationCellStatus();
    net_update(y_st);
    bin_update(N);
    InitializationCostFunctionGradient(&sum_wgrad, &sum_pgrad);
    InitializationCoefficients();
    if(DEN_ONLY_PRECON) {
        InitializationPrecondition_DEN_ONLY_PRECON();
    }
    else {
        InitializationPrecondition();
    }
    InitializationIter();
    z_init();
    if(dynamicStepCMD && isTrial == false && isFirst_gp_opt == true) {
        UPPER_PCOF = 1.0001;
        potnPhaseDS = potnPhase1;
    }
    if(isTrial == true)
        UPPER_PCOF = 1.05;
    net_update(z_st);
    bin_update(N);  // igkang
    it->tot_wwl = (1.0 - opt_w_cof) * it->tot_hpwl + opt_w_cof * it->tot_stnwl;
    getCostFuncGradient2(z_dst, z_wdst, z_pdst, z_pdstl, N, cellLambdaArr, false, false);
    a = 1.0;
    get_lc(y_st, y_dst, z_st, z_dst, &it0, N);
    it->alpha00 = it0.alpha00;
    if(isTrial) {
        initialOVFL = it->ovfl;
        trial_HPWLs.push_back(std::pair< prec, prec >(it->tot_hpwl, 0.0));
        trial_POTNs.push_back(it->potn);
    }
    PrintNesterovOptStatus(0);
    last_iter = (isPreplace) ? __optimize__() : DoNesterovOptimization();

    SummarizeNesterovOpt(last_iter);
    mkl_malloc_free();
    profile.ns += time_since(nsS);
}

void myNesterov::InitializationCommonVar() {
    N = gcell_cnt;
    N_org = moduleCNT;
    start_idx = 0;
    end_idx = N;

    if(dynamicStepCMD)
        max_iter = 6000;
    else
        max_iter = 2500;

    //debug : mgwoo
//    max_iter = 3;

    last_ra_iter = 0;
    a = 0;
    ab = 0;
    alpha = 0;
    cof = 0;
    initialOVFL = 0;
    alpha_pred = 0;
    alpha_new = 0;
    post_filler = 0;
    sum_wgrad = 0;
    sum_pgrad = 0;
    sum_tgrad = 0;
    u.SetZero();
    v.SetZero();
    half_densize.SetZero();
    wcof.SetZero();
    wpre.SetZero();
    charge_dpre.SetZero();
    temp_dpre.SetZero();
    pre.SetZero();

    iter_st =
        (struct ITER *)mkl_malloc(sizeof(struct ITER) * (max_iter + 1), 64);

    x_st = (struct FPOS *)mkl_malloc(sizeof(struct FPOS) * N, 64);

    y_st = (struct FPOS *)mkl_malloc(sizeof(struct FPOS) * N, 64);
    y_dst = (struct FPOS *)mkl_malloc(sizeof(struct FPOS) * N, 64);
    y_wdst = (struct FPOS *)mkl_malloc(sizeof(struct FPOS) * N, 64);
    y_pdst = (struct FPOS *)mkl_malloc(sizeof(struct FPOS) * N, 64);
    y_pdstl = (struct FPOS *)mkl_malloc(sizeof(struct FPOS) * N, 64);

    z_st = (struct FPOS *)mkl_malloc(sizeof(struct FPOS) * N, 64);
    z_dst = (struct FPOS *)mkl_malloc(sizeof(struct FPOS) * N, 64);
    z_wdst = (struct FPOS *)mkl_malloc(sizeof(struct FPOS) * N, 64);
    z_pdst = (struct FPOS *)mkl_malloc(sizeof(struct FPOS) * N, 64);
    z_pdstl = (struct FPOS *)mkl_malloc(sizeof(struct FPOS) * N, 64);

    x0_st = (struct FPOS *)mkl_malloc(sizeof(struct FPOS) * N, 64);

    y0_st = (struct FPOS *)mkl_malloc(sizeof(struct FPOS) * N, 64);
    y0_dst = (struct FPOS *)mkl_malloc(sizeof(struct FPOS) * N, 64);
    y0_wdst = (struct FPOS *)mkl_malloc(sizeof(struct FPOS) * N, 64);
    y0_pdst = (struct FPOS *)mkl_malloc(sizeof(struct FPOS) * N, 64);
    y0_pdstl = (struct FPOS *)mkl_malloc(sizeof(struct FPOS) * N, 64);

    memset(x_st, 0.0f, sizeof(prec) * N * 3);
    memset(y_st, 0.0f, sizeof(prec) * N * 3);
    memset(y_dst, 0.0f, sizeof(prec) * N * 3);
    memset(y_wdst, 0.0f, sizeof(prec) * N * 3);
    memset(y_pdst, 0.0f, sizeof(prec) * N * 3);
    memset(y_pdstl, 0.0f, sizeof(prec) * N * 3);

    memset(z_st, 0.0f, sizeof(prec) * N * 3);
    memset(z_dst, 0.0f, sizeof(prec) * N * 3);
    memset(z_wdst, 0.0f, sizeof(prec) * N * 3);
    memset(z_pdst, 0.0f, sizeof(prec) * N * 3);
    memset(z_pdstl, 0.0f, sizeof(prec) * N * 3);

    memset(x0_st, 0.0f, sizeof(prec) * N * 3);
    memset(y0_st, 0.0f, sizeof(prec) * N * 3);
    memset(y0_dst, 0.0f, sizeof(prec) * N * 3);
    memset(y0_wdst, 0.0f, sizeof(prec) * N * 3);
    memset(y0_pdst, 0.0f, sizeof(prec) * N * 3);
    memset(y0_pdstl, 0.0f, sizeof(prec) * N * 3);

    cellLambdaArr = (prec *)mkl_malloc(sizeof(prec) * N, 64);
    pcofArr = (prec *)mkl_malloc(sizeof(prec) * 100, 64);
    MIN_PRE = 1;
}

void myNesterov::InitializationCellStatus() {
    // OPT_INPUT == QWL_ISOL
    cg_input(x_st, N, OPT_INPUT);
    net_update(x_st);

    for(int i = 0; i < N; i++) {
        y_st[i] = y0_st[i] = x0_st[i] = x_st[i];
    }

    if(STAGE == mGP3D || STAGE == cGP3D)
        ShiftPL_SA(y_st, N);

    if(numLayer == 1) {
        if(placementMacroCNT == 0) {
            if(STAGE == cGP2D)
                ShiftPL_SA(y_st, N);
        }
        else {
            if(STAGE == mGP2D)
                ShiftPL_SA(y_st, N);
        }
    }
}

void myNesterov::InitializationCoefficients() {
    if(routabilityCMD)
        INIT_LAMBDA_COF_GP = 0.001;
    // if (routabilityCMD) INIT_LAMBDA_COF_GP = 0.1;
    else {
//        INIT_LAMBDA_COF_GP = 0.0001;
        INIT_LAMBDA_COF_GP = 0.00008;
    }

    if(STAGE == mGP3D) {
        opt_phi_cof = sum_wgrad / sum_pgrad * INIT_LAMBDA_COF_GP;
        opt_w_cof = 0;
        ALPHA = ALPHAmGP;  // 1E-12;
        BETA = BETAmGP;    // 1E-16;
        wcof = get_wlen_cof(gsum_ovfl);
        wlen_cof = fp_mul(base_wcof, wcof);
        wlen_cof_inv = fp_inv(wlen_cof);
    }
    else if(STAGE == mGP2D) {
        opt_phi_cof = sum_wgrad / sum_pgrad * INIT_LAMBDA_COF_GP;
        ALPHA = ALPHAmGP;
        if(INPUT_FLG == MMS &&
           (placementMacroCNT >= 3000 ||
            (!dynamicStepCMD &&
             (target_cell_den <= 0.60f || target_cell_den <= 0.60)))) {
            // Wants to have more inertia
            BETA = 1E-11;
            dampParam = 0.999900;
        }
        else {
            BETA = BETAmGP;
        }
        wcof = get_wlen_cof(gsum_ovfl);
        wlen_cof = fp_mul(base_wcof, wcof);
        wlen_cof_inv = fp_inv(wlen_cof);
    }
    else if(STAGE == cGP3D) {
        if(placementMacroCNT > 0) {
            cGP3D_buf_iter = mGP3D_tot_iter / 10;
            opt_phi_cof =
                mGP3D_opt_phi_cof / pow(UPPER_PCOF, (prec)cGP3D_buf_iter);
        }
        else if(placementMacroCNT == 0) {
            opt_phi_cof = sum_wgrad / sum_pgrad * INIT_LAMBDA_COF_GP;
        }
        ALPHA = ALPHAcGP;  // 1E-12;
        BETA = BETAcGP;    // 1E-16;
        wcof = get_wlen_cof(gsum_ovfl);
        wlen_cof = fp_mul(base_wcof, wcof);
        wlen_cof_inv = fp_inv(wlen_cof);
    }
    else if(STAGE == cGP2D) {
        if(placementMacroCNT > 0) {
            if(INPUT_FLG == MMS)
                cGP2D_buf_iter = mGP2D_tot_iter / 10;
            else
                cGP2D_buf_iter = mGP2D_tot_iter / 8;

            if(dynamicStepCMD && !constraintDrivenCMD) {
                UPPER_PCOF = 1.010;  // 1.010
                LOWER_PCOF = 0.60;   // 0.50
            }
            else if(dynamicStepCMD && constraintDrivenCMD) {
                UPPER_PCOF = 1.015;  // 1.015,1.012
                LOWER_PCOF = 0.90;   // 0.70, 0.60
            }

            opt_phi_cof =
                mGP2D_opt_phi_cof / pow(UPPER_PCOF, (prec)cGP2D_buf_iter);
        }
        else {
            opt_phi_cof = sum_wgrad / sum_pgrad * INIT_LAMBDA_COF_GP;
        }
        opt_w_cof = 0;
        if(placementMacroCNT == 0) {
            ALPHA = ALPHAcGP;  // 1E-15;
            BETA = BETAcGP;    // 1E-14;
        }
        else {
            ALPHA = ALPHAcGP;
            BETA = BETAcGP;
        }
        wcof = get_wlen_cof(gsum_ovfl);
        wlen_cof = fp_mul(base_wcof, wcof);
        wlen_cof_inv = fp_inv(wlen_cof);
        // IK
        if(!routabilityCMD)
            overflowMin = 0.07;
    }

    if(lambda2CMD == true) {
        opt_phi_cof_local = opt_phi_cof;
        for(int i = 0; i < N; i++) {
            cellLambdaArr[i] = opt_phi_cof_local;
        }
    }
}

void myNesterov::InitializationPrecondition() {
    if(lambda2CMD == false) {
        for(int i = 0; i < N; i++) {
            y_dst[i].x = y_wdst[i].x + opt_phi_cof * y_pdst[i].x;
            y_dst[i].y = y_wdst[i].y + opt_phi_cof * y_pdst[i].y;
            wlen_pre(i, &wpre);
            potn_pre(i, &charge_dpre);

            pre.x = wpre.x + opt_phi_cof * charge_dpre.x;
            pre.y = wpre.y + opt_phi_cof * charge_dpre.y;

            if(pre.x < MIN_PRE)
                pre.x = MIN_PRE;
            if(pre.y < MIN_PRE)
                pre.y = MIN_PRE;

            y_dst[i].x /= pre.x;
            y_dst[i].y /= pre.y;
        }
    }
    else if(lambda2CMD == true) {
        for(int i = 0; i < N; i++) {
            y_dst[i].x = y_wdst[i].x + opt_phi_cof * y_pdst[i].x +
                         cellLambdaArr[i] * y_pdstl[i].x;
            y_dst[i].y = y_wdst[i].y + opt_phi_cof * y_pdst[i].y +
                         cellLambdaArr[i] * y_pdstl[i].y;

            wlen_pre(i, &wpre);
            potn_pre(i, &charge_dpre);

            pre.x = wpre.x + opt_phi_cof * charge_dpre.x;
            pre.y = wpre.y + opt_phi_cof * charge_dpre.y;

            if(pre.x < MIN_PRE)
                pre.x = MIN_PRE;
            if(pre.y < MIN_PRE)
                pre.y = MIN_PRE;

            y_dst[i].x /= pre.x;
            y_dst[i].y /= pre.y;
        }
    }
    else {
    }
}

void myNesterov::InitializationPrecondition_DEN_ONLY_PRECON() {
    if(lambda2CMD == false) {
        for(int i = 0; i < N; i++) {
            y_dst[i].x = y_wdst[i].x + opt_phi_cof * y_pdst[i].x;
            y_dst[i].y = y_wdst[i].y + opt_phi_cof * y_pdst[i].y;

            wlen_pre(i, &wpre);
            potn_pre(i, &charge_dpre);
            pre = charge_dpre;
            if(pre.x < MIN_PRE)
                pre.x = MIN_PRE;
            if(pre.y < MIN_PRE)
                pre.y = MIN_PRE;
            y_dst[i].x /= pre.x;
            y_dst[i].y /= pre.y;
        }
    }
    else if(lambda2CMD == true) {
        for(int i = 0; i < N; i++) {
            y_dst[i].x = y_wdst[i].x + opt_phi_cof * y_pdst[i].x +
                         cellLambdaArr[i] * y_pdstl[i].x;
            y_dst[i].y = y_wdst[i].y + opt_phi_cof * y_pdst[i].y +
                         cellLambdaArr[i] * y_pdstl[i].y;
            wlen_pre(i, &wpre);
            potn_pre(i, &charge_dpre);
            pre = charge_dpre;
            if(pre.x < MIN_PRE)
                pre.x = MIN_PRE;
            if(pre.y < MIN_PRE)
                pre.y = MIN_PRE;
            y_dst[i].x /= pre.x;
            y_dst[i].y /= pre.y;
        }
    }
    else {
    }
}

void myNesterov::InitializationIter() {
    init_iter(iter_st, 0);
    it = &iter_st[0];
    time_calc(0, &it->cpu_curr, &it->cpu_cost);
    it->cpu_cost = 0;
    it->tot_wlen = get_wlen();
    it->wlen = total_wlen;
    it->tot_hpwl = GetHpwl();
    it->hpwl = total_hpwl;
    it->potn = gsum_phi;
    it->ovfl = gsum_ovfl;
    global_ovfl = gsum_ovfl;
    it->grad = get_norm(y_dst, N, 2.0);
}



void myNesterov::__post_filler__(int i, cell_den_t* cells, bin_t* bins, cell_phy_t* ios, net_t* mesh) {
        FILLER_PLACE = 0;
        if((isFirst_gp_opt == false) && post_filler == 0) {
            if(i < NUM_ITER_FILLER_PLACE) {
                FILLER_PLACE = 1;
                post_filler = 0;
                start_idx = moduleCNT;
                end_idx = N;
            }
            else {
                FILLER_PLACE = 0;
                post_filler = 1;
                start_idx = 0;
                end_idx = N;
//                __getCostFunc__(y_dst, y_wdst, y_pdst, y_pdstl, N, cellLambdaArr, DEN_ONLY_PRECON, FILLER_PLACE, cells, bins, ios, mesh);
                get_lc(y_st, y_dst, z_st, z_dst, &it0, N);
            }
        }
}



int myNesterov::__optimize__() {
    int i;
    prec minPotn = PREC_MAX;
    TIER* tier = &tier_st[0];
    temp_iter = 0;
    bool timeon = false;
    double time = 0.0f;
    circuit_t* circuit = circuit_t(numThread)
    .withNets(netInstance, netCNT)
    ->withCells(gcell_st, N)
    ->withModules(moduleInstance, moduleCNT)
    ->withBins(tier->bin_mat, tier->tot_bin_cnt);

    size_t* cellPinWk = schedule(refIo(circuit->cells, N), numThread);
	size_t* netPinWk = schedule(refNets(circuit->nets, netCNT), numThread);
    for(i = 0; i < max_iter; i++) {
        if( timeon ) time_start(&time);
        it = &iter_st[i + 1];
        init_iter(it, i + 1);
        FILLER_PLACE = 0;
        if((isFirst_gp_opt == false) && post_filler == 0) {
            if(i < NUM_ITER_FILLER_PLACE) {
                FILLER_PLACE = 1;
                post_filler = 0;
                start_idx = moduleCNT;
                end_idx = N;
            }
            else {
                FILLER_PLACE = 0;
                post_filler = 1;
                start_idx = 0;
                end_idx = N;
                exit(1);
                getCostFuncGradient3(y_dst, y_wdst, y_pdst, y_pdstl, N, cellLambdaArr);
                /*dCells.Copy(cells);
                __gradient__(y_dst, y_wdst, y_pdst, y_pdstl, N, cellLambdaArr, DEN_ONLY_PRECON, cells, bins, ios, mesh, &it->gradcost, cellPinWk, &dCells);   
                dCells.Copyback(cells);*/
                get_lc(y_st, y_dst, z_st, z_dst, &it0, N);
            }
        }
        it->lc = it0.lc;
        alpha = it->alpha00 = it0.alpha00;
        it->alpha00 = it0.alpha00 = alpha;
        ab = a;
        a = (1.0 + sqrt(4.0 * a * a + 1.0)) * 0.5;
        cof = (ab - 1.0) / a;
        alpha_pred = it->alpha00;
        backtrack_cnt = 0;
        while(1) {

            
            backtrack_cnt++;

            if( timeon ) { time_start(&time); };
            int j=0;
            FPOS u, v;
            for(j = start_idx; j < end_idx; j++) {
                half_densize = gcell_st[j].half_den_size;

                u.x = y_st[j].x + alpha_pred * y_dst[j].x;
                u.y = y_st[j].y + alpha_pred * y_dst[j].y;
                v.x = u.x + cof * (u.x - x_st[j].x);
                v.y = u.y + cof * (u.y - x_st[j].y);

                x0_st[j] = valid_coor00(u, half_densize);
                y0_st[j] = valid_coor00(v, half_densize);
            }
            auto t0 = std::chrono::steady_clock::now();
            update_wirelength(circuit, y0_st);
            profile.wlen += time_since(t0);
            for(int i = 0; i < tier->tot_bin_cnt; i++) {
                bin_t* bp = &circuit->bins[i];
                bp->cellArea = 0;
                bp->fillerArea = 0;
            }
            update_field_potential(circuit);
            t0 = std::chrono::steady_clock::now();
            if(FILLER_PLACE) update_fill_gradient(y0_dst, y0_wdst, y0_pdst, y0_pdstl, N, cellLambdaArr, DEN_ONLY_PRECON, circuit);
            else update_gradient(y0_dst, y0_wdst, y0_pdst, y0_pdstl, N, cellLambdaArr, DEN_ONLY_PRECON, circuit);   
            profile.cost += time_since(t0);
            get_lc(y_st, y_dst, y0_st, y0_dst, &it0, N);

            
            alpha_new = it0.alpha00;

            if(alpha_new > alpha_pred * 0.95 ||
               backtrack_cnt >= MAX_BKTRK_CNT) {
                alpha_pred = alpha_new;
                it->alpha00 = alpha_new;
                break;
            }
            else {
                alpha_pred = alpha_new;
            }
        }

        UpdateNesterovOptStatus();
        __UpdateNesterovIter__(circuit->nets, i + 1, it, &iter_st[i]);

        if(dynamicStepCMD && !isTrial) {
            if(isFirst_gp_opt)
                potnPhaseDS = definePOTNphase(it->potn);
            else
                potnPhaseDS = potnPhase8;

            stepSizeAdaptation_by2ndOrderEPs(it->tot_hpwl);
        }

        if(isTrial) {
            if(it->ovfl > initialOVFL / 2.5) {
                trial_HPWLs.push_back(std::pair< prec, prec >(it->tot_hpwl, 0.0));
                trial_POTNs.push_back(it->potn);
            }
            else {
                trial_iterCNT = trial_HPWLs.size();
                printf("\n");
                printf("INFO:    SUMMARY tGP\n");
                printf("INFO:    #iterations = %d\n", trial_iterCNT);
                return i;
            }
        }

        // Update ALPHA and BETA for Local Density Function
        if(constraintDrivenCMD == true)
            UpdateAlpha(it);
        if(lambda2CMD == true)
            UpdateBeta(it);
        // Elimination Condition
        if(it->tot_hpwl > 2000000000)
            exit(0);

        // Termination Condition 1
        if(it->ovfl <= overflowMin && i > 220)
            return i;

        // Termination Condition 2
        if(STAGE == cGP2D && i > 220) {
            if((it->ovfl <= 0.13f && dynamicStepCMD) || (it->ovfl <= 0.10f)) {
                if(minPotn * 1.01 <= it->potn && temp_iter == 0) {
                    temp_iter++;
                }
                if(minPotn > it->potn) {
                    minPotn = it->potn;
                    temp_iter = 0;
                }
            }
            if(temp_iter == 10) {
                return i;
            }
            if(temp_iter) {
                temp_iter++;
            }
        }
    }
    circuit->destroy(gcell_st, netInstance, tier->bin_mat);
    //Put back hotspots
    return -1;
}

int myNesterov::DoNesterovOptimization() {
    int i;
    prec minPotn = PREC_MAX;
    TIER* tier = &tier_st[0];
    temp_iter = 0;
    bool timeon = false;
    double time = 0.0f;

    for(i = 0; i < max_iter; i++) {
        if( timeon ) time_start(&time);
        if(isTrial == false && routabilityCMD == true &&
           isRoutabilityInit == false) {
            routability_init();
        }
        it = &iter_st[i + 1];
        init_iter(it, i + 1);
        FILLER_PLACE = 0;
        if((isFirst_gp_opt == false) && post_filler == 0) {
            if(i < NUM_ITER_FILLER_PLACE) {
                FILLER_PLACE = 1;
                post_filler = 0;
                start_idx = moduleCNT;
                end_idx = N;
            }
            else {
                FILLER_PLACE = 0;
                post_filler = 1;
                start_idx = 0;
                end_idx = N;
                getCostFuncGradient3(y_dst, y_wdst, y_pdst, y_pdstl, N,
                                     cellLambdaArr);
                get_lc(y_st, y_dst, z_st, z_dst, &it0, N);
            }
        }
        it->lc = it0.lc;
        alpha = it->alpha00 = it0.alpha00;
        it->alpha00 = it0.alpha00 = alpha;
            
        ab = a;
        a = (1.0 + sqrt(4.0 * a * a + 1.0)) * 0.5;
        cof = (ab - 1.0) / a;

        alpha_pred = it->alpha00;
        backtrack_cnt = 0;

        while(1) {
            backtrack_cnt++;

            int j=0;
            FPOS u, v;
            for(j = start_idx; j < end_idx; j++) {
                half_densize = gcell_st[j].half_den_size;

                u.x = y_st[j].x + alpha_pred * y_dst[j].x;
                u.y = y_st[j].y + alpha_pred * y_dst[j].y;
                v.x = u.x + cof * (u.x - x_st[j].x);
                v.y = u.y + cof * (u.y - x_st[j].y);

                x0_st[j] = valid_coor00(u, half_densize);
                y0_st[j] = valid_coor00(v, half_densize);
            }
            auto t0 = std::chrono::steady_clock::now();
            net_update(y0_st);
            profile.wlen += time_since(t0);
            TIER* tier = &tier_st[0];
            for(int i = 0; i < tier->tot_bin_cnt; i++) {
                BIN* bp = &tier->bin_mat[i];
                bp->cell_area = 0;
                bp->cell_area2 = 0;
            }
            bin_update7_cGP2D(&profile.density, &profile.fft);
            t0 = std::chrono::steady_clock::now();
            getCostFuncGradient3(y0_dst, y0_wdst, y0_pdst, y0_pdstl, N, cellLambdaArr);
            profile.cost += time_since(t0);

            get_lc(y_st, y_dst, y0_st, y0_dst, &it0, N);

            alpha_new = it0.alpha00;

            if(alpha_new > alpha_pred * 0.95 ||
               backtrack_cnt >= MAX_BKTRK_CNT) {
                alpha_pred = alpha_new;
                it->alpha00 = alpha_new;
                break;
            }
            else {
                alpha_pred = alpha_new;
            }
        }

        UpdateNesterovOptStatus();
        UpdateNesterovIter(i + 1, it, &iter_st[i]);

        if(dynamicStepCMD && !isTrial) {
            if(isFirst_gp_opt)
                potnPhaseDS = definePOTNphase(it->potn);
            else
                potnPhaseDS = potnPhase8;

            stepSizeAdaptation_by2ndOrderEPs(it->tot_hpwl);
        }

        if(isTrial) {
            if(it->ovfl > initialOVFL / 2.5) {
                trial_HPWLs.push_back(
                    std::pair< prec, prec >(it->tot_hpwl, 0.0));
                trial_POTNs.push_back(it->potn);
            }
            else {
                trial_iterCNT = trial_HPWLs.size();
                printf("\n");
                printf("INFO:    SUMMARY tGP\n");
                printf("INFO:    #iterations = %d\n", trial_iterCNT);
                return i;
            }
        }

        if(routabilityCMD == true && STAGE == cGP2D) {
            // LW mod 10/20/16 temp_con_orig = 0.12+
            // igkang
            // LW 05/30/17
            prec temp_con = 0.10 + bloating_max_count / 10.0 - bloatCNT / 10.0;
            if(DEN_ONLY_PRECON) {
                temp_con = 0.15 + bloating_max_count / 10.0 - bloatCNT / 10.0;
            }

            if((it->ovfl < temp_con) && (i - last_ra_iter > 10)) {
                // UPPER_PCOF = 1.01;
                if(bloatCNT < bloating_max_count) {
                    last_ra_iter = i;
                    cell_update(x_st, N_org);
                    modu_copy();
                    congEstimation(x_st);
                    // if (inflation_cnt == 0) calcCong_print_detail();
                    if(inflation_cnt % 2 == 0) {
                        is_inflation_h = true;
                    }
                    else {
                        is_inflation_h = false;
                    }
                    if(flg_noroute) {
                        inflation_cnt = 100;
                        if(inflation_cnt >= inflation_max_cnt) {
                            bloatCNT++;
                            inflation_cnt = 0;
                        }
                    }
                    else {
                        routability();
                        inflation_cnt++;
                        if(inflation_cnt >= inflation_max_cnt) {
                            bloatCNT++;
                            inflation_cnt = 0;
                        }
                        isBloatingBegin = true;
                        if(inflation_cnt == 1) {
                            before100iter_cof = pcofArr[(i + 2) % 100];
                        }
                        opt_phi_cof = before100iter_cof;
                    }
                }
            }
        }

        // Update ALPHA and BETA for Local Density Function
        if(constraintDrivenCMD == true)
            UpdateAlpha(it);
        if(lambda2CMD == true)
            UpdateBeta(it);
        // Elimination Condition
        if(it->tot_hpwl > 2000000000)
            exit(0);

        // Termination Condition 1
        if(it->ovfl <= overflowMin && i > 220)
            return i;

        // Termination Condition 2
        if(STAGE == cGP2D && i > 220) {
            if((it->ovfl <= 0.13f && dynamicStepCMD) || (it->ovfl <= 0.10f)) {
                if(minPotn * 1.01 <= it->potn && temp_iter == 0) {
                    temp_iter++;
                }
                if(minPotn > it->potn) {
                    minPotn = it->potn;
                    temp_iter = 0;
                }
            }
            if(temp_iter == 10) {
                return i;
            }
            if(temp_iter) {
                temp_iter++;
            }
        }
    }

    //Put back hotspots
    return -1;
}

void myNesterov::mkl_malloc_free() {
    mkl_free(iter_st);
    mkl_free(x_st);
    mkl_free(y_st);
    mkl_free(y_dst);
    mkl_free(y_wdst);
    mkl_free(y_pdst);
    mkl_free(y_pdstl);
    mkl_free(z_st);
    mkl_free(z_dst);
    mkl_free(z_wdst);
    mkl_free(z_pdst);
    mkl_free(z_pdstl);
    mkl_free(x0_st);
    mkl_free(y0_st);
    mkl_free(y0_dst);
    mkl_free(y0_wdst);
    mkl_free(y0_pdst);
    mkl_free(y0_pdstl);
    mkl_free(cellLambdaArr);
    mkl_free(pcofArr);
}

void myNesterov::SummarizeNesterovOpt(int last_index) {
    prec tot_hpwl_y;
    prec tot_hpwl_x;

    if(STAGE == mGP3D) {
        mGP3D_iterCNT = last_index + 1;
        hpwl_mGP3D = it->tot_hpwl;
        printf("\n");
        printf("INFO:    SUMMARY mGP3D\n");
        printf("INFO:    #iterations = %d\n", mGP3D_iterCNT);
        printf("INFO:        mGP3D ITERATION HPWL %.4lf\n", it->tot_hpwl);
        printf("INFO:        mGP3D ITERATION OVFL %.4lf\n", it->ovfl);
        printf("INFO:        mGP3D      POTENTIAL %.4lf\n", it->potn);
        mGP3D_tot_iter = last_index;
        mGP3D_opt_phi_cof = opt_phi_cof;
    }
    else if(STAGE == mGP2D) {
        mGP2D_iterCNT = last_index + 1;
        hpwl_mGP2D = it->tot_hpwl;
        printf("\n");
        printf("INFO:    SUMMARY mGP2D\n");
        printf("INFO:    #iterations = %d\n", mGP2D_iterCNT);
        printf("INFO:        mGP2D ITERATION HPWL %.4lf\n", it->tot_hpwl);
        printf("INFO:        mGP2D ITERATION OVFL %.4lf\n", it->ovfl);
        printf("INFO:        mGP2D      POTENTIAL %.4lf\n", it->potn);
        mGP2D_tot_iter = last_index;
        mGP2D_opt_phi_cof = opt_phi_cof;
    }
    else if(STAGE == cGP3D) {
        cGP3D_iterCNT = last_index + 1;
        hpwl_cGP3D = it->tot_hpwl;
        printf("\n");
        printf("INFO:    SUMMARY cGP3D\n");
        printf("INFO:    #iterations = %d\n", cGP3D_iterCNT);
        printf("INFO:        cGP3D ITERATION HPWL %.4lf\n", it->tot_hpwl);
        printf("INFO:        cGP3D ITERATION OVFL %.4lf\n", it->ovfl);
        printf("INFO:        cGP3D      POTENTIAL %.4lf\n", it->potn);
        cGP3D_tot_iter = last_index;
        cGP3D_opt_phi_cof = opt_phi_cof;
    }
    else if(STAGE == cGP2D) {
        cGP2D_iterCNT = last_index + 1;
        hpwl_cGP2D = it->tot_hpwl;
        printf("\n");
        printf("INFO:    SUMMARY cGP2D\n");
        printf("INFO:    #iterations = %d\n", cGP2D_iterCNT);
        printf("INFO:        cGP2D ITERATION HPWL %.4lf\n", it->tot_hpwl);
        printf("INFO:        cGP2D ITERATION OVFL %.4lf\n", it->ovfl);
        printf("INFO:        cGP2D      POTENTIAL %.4lf\n", it->potn);
        cGP2D_tot_iter = last_index;
        cGP2D_opt_phi_cof = opt_phi_cof;
    }
    else {
    }

    cell_update(x_st, N_org);

    net_update(y_st);
    tot_hpwl_y = GetHpwl();

    net_update(x_st);
    tot_hpwl_x = GetHpwl();

    printf("INFO:    TOTAL HPWL (U_k, V_k) = %.6lf, %.6lf\n", tot_hpwl_x,
           tot_hpwl_y);
}

void getCostFuncGradient3(struct FPOS *dst, struct FPOS *wdst,
                          struct FPOS *pdst, struct FPOS *pdstl, int N,
                          prec *cellLambdaArr) {

        getCostFuncGradient2(dst, wdst, pdst, pdstl, N, cellLambdaArr, DEN_ONLY_PRECON, FILLER_PLACE);

}


void getCostFuncGradient2(struct FPOS *dst, struct FPOS *wdst,
                          struct FPOS *pdst, struct FPOS *pdstl, int N,
                          prec *cellLambdaArr, bool onlyPreCon, bool filler) {
    CELL* cell = NULL;
    struct FPOS wgrad;
    struct FPOS pgrad;
    struct FPOS pgradl;
    struct FPOS wpre;
    struct FPOS charge_dpre;
    struct FPOS pre;
    for(int i = (filler ? moduleCNT : 0); i < N; i++) {
        cell = &gcell_st[i];
        if(filler) cellLambdaArr[i] *= dampParam;
        if(!filler && cell->flg == Macro && (STAGE == cGP3D || STAGE == cGP2D)) {
            wgrad.SetZero();
            pgrad.SetZero();
            pgradl.SetZero();
        }
        else if(filler && cell->flg != FillerCell) continue;
        else {
            if(!filler) wlen_grad(i, &wgrad);
            if(STAGE == mGP2D) {
                if(constraintDrivenCMD == false)
                    potn_grad_2D(i, &pgrad);
                else if(constraintDrivenCMD == true) {
                    exit(1);
                    potn_grad_2D(i, &pgrad);
                    potn_grad_2D_local(i, &pgradl, &cellLambdaArr[i]);
                }
            }
            else if(STAGE == cGP2D) {
                if(constraintDrivenCMD == false)
                    potn_grad_2D(i, &pgrad);
                else if(constraintDrivenCMD == true) {
                    exit(1);
                    potn_grad_2D(i, &pgrad);
                    potn_grad_2D_local(i, &pgradl, &cellLambdaArr[i]);
                }
            }
        }

        if(filler) wdst[i] = wgrad; else wdst[i].SetZero();
        pdst[i] = pgrad;

        dst[i].x = (filler ? 0 : wgrad.x) + opt_phi_cof * pgrad.x;
        dst[i].y = (filler ? 0 : wgrad.y) + opt_phi_cof * pgrad.y;

        if(lambda2CMD == true) {
            pdstl[i] = pgradl;
            dst[i].x += cellLambdaArr[i] * pgradl.x;
            dst[i].y += cellLambdaArr[i] * pgradl.y;
        }

        if(!filler) wlen_pre(i, &wpre);
        potn_pre(i, &charge_dpre);
        if(onlyPreCon) pre = charge_dpre;
        else { //?
          pre.x = (filler ? 0 : wpre.x) + opt_phi_cof * charge_dpre.x;
          pre.y = (filler ? 0 : wpre.y) + opt_phi_cof * charge_dpre.y;
        }
        if(pre.x < MIN_PRE)
            pre.x = MIN_PRE;
        if(pre.y < MIN_PRE)
            pre.y = MIN_PRE;

        dst[i].x /= pre.x;
        dst[i].y /= pre.y;
    }
}

void myNesterov::z_init() {
    FPOS half_densize;
    prec zx = 0.0f, zy = 0.0f, zz = 0.0f;
    prec coeffi = z_ref_alpha * GP_SCAL;

    for(int j = start_idx; j < end_idx; j++) {
        half_densize = gcell_st[j].half_den_size;

        if(GP_DIM_ONE) {
            zx = y_st[j].x + place_backup.cnt.x * coeffi * y_dst[j].x;
            zy = y_st[j].y + place_backup.cnt.y * coeffi * y_dst[j].y;
            zz = y_st[j].z + place_backup.cnt.z * coeffi * y_dst[j].z;
        }
        else {
            zx = y_st[j].x + z_ref_alpha * y_dst[j].x;
            zy = y_st[j].y + z_ref_alpha * y_dst[j].y;
            zz = y_st[j].z + z_ref_alpha * y_dst[j].z;
        }
        z_st[j].x = valid_coor2(zx, half_densize.x, 0);
        z_st[j].y = valid_coor2(zy, half_densize.y, 1);
        z_st[j].z = valid_coor2(zz, half_densize.z, 2);
    }
}

void myNesterov::UpdateNesterovOptStatus() {
    std::swap(z_st, y_st);
    std::swap(z_dst, y_dst);
    std::swap(z_wdst, y_wdst);
    std::swap(z_pdst, y_pdst);
    std::swap(z_pdstl, y_pdstl);
    std::swap(x_st, x0_st);
    std::swap(y_st, y0_st);
    std::swap(y_dst, y0_dst);
    std::swap(y_wdst, y0_wdst);
    std::swap(y_pdst, y0_pdst);
    std::swap(y_pdstl, y0_pdstl);
}

// this only calculate HPWL based on the stored value in NET's ure
prec __GetHpwl__(net_t* nets) {
    total_hpwl.SetZero();
    total_stnwl.SetZero();
    for(int i = 0; i < netCNT; i++) {
        net_t *curNet = &nets[i];
        if(curNet->pinCNT > 1){
        total_hpwl.x += curNet->max.x - curNet->min.x;
        total_hpwl.y += curNet->max.y - curNet->min.y;
        }
    }
    if(GP_DIM_ONE) {
        total_hpwl.x *= place_backup.cnt.x * GP_SCAL;
        total_hpwl.y *= place_backup.cnt.y * GP_SCAL;
    }
    return total_hpwl.x + total_hpwl.y;
}


void myNesterov::UpdateNesterovIter(int iter, struct ITER *it, struct ITER *last_it) {
    it->grad = get_norm(y_dst, N, 2.0);
    it->potn = gsum_phi;
    it->ovfl = gsum_ovfl;
    it->dis00 = get_dis(z_st, y_st, N);
    it->wcof = get_wlen_cof(it->ovfl);
    wlen_cof = fp_mul(base_wcof, it->wcof);
    wlen_cof_inv = fp_inv(wlen_cof);
    pcofArr[iter % 100] = opt_phi_cof;
    if(!FILLER_PLACE) {
        it->tot_hpwl = GetHpwl();
        it->tot_stnwl = total_stnwl.x + total_stnwl.y;  // lutong
        opt_w_cof = it->ovfl <= stn_weight ? 1 : (1.0 - it->ovfl) / (1.0 - stn_weight);
        it->tot_wwl = (1.0 - opt_w_cof) * it->tot_hpwl + opt_w_cof * it->tot_stnwl;  // lutong
        it->hpwl = total_hpwl;
        it->pcof = get_phi_cof1((it->tot_hpwl - last_it->tot_hpwl) / ref_dwl0);
        if(temp_iter < 1) {
            opt_phi_cof *= it->pcof;  // igkang
        }
        if(lambda2CMD == true) {
            opt_phi_cof_local = 0;
            for(int i = 0; i < N; i++) {
                opt_phi_cof_local += (cellLambdaArr[i]);
            }
            opt_phi_cof_local /= (prec)N;
        }
        it->wlen = total_wlen;
    }
    else {
        it->hpwl = last_it->hpwl;
        it->tot_hpwl = last_it->tot_hpwl;
        it->tot_stnwl = last_it->tot_stnwl;  // lutong
        it->tot_wwl = last_it->tot_wwl;      // lutong
        it->wlen.SetZero();
        it->tot_wlen = 0;
    }
    time_calc(last_it->cpu_curr, &it->cpu_curr, &it->cpu_cost);
    PrintNesterovOptStatus(iter);
    fflush(stdout);
}
void myNesterov::__UpdateNesterovIter__(net_t* nets, int iter, struct ITER *it, struct ITER *last_it) {
    it->grad = get_norm(y_dst, N, 2.0);
    it->potn = gsum_phi;
    it->ovfl = gsum_ovfl;
    it->dis00 = get_dis(z_st, y_st, N);
    it->wcof = get_wlen_cof(it->ovfl);
    wlen_cof = fp_mul(base_wcof, it->wcof);
    wlen_cof_inv = fp_inv(wlen_cof);
    pcofArr[iter % 100] = opt_phi_cof;
    if(!FILLER_PLACE) {
        it->tot_hpwl = __GetHpwl__(nets);
        it->tot_stnwl = total_stnwl.x + total_stnwl.y;  // lutong
        opt_w_cof = it->ovfl <= stn_weight ? 1 : (1.0 - it->ovfl) / (1.0 - stn_weight);
        it->tot_wwl = (1.0 - opt_w_cof) * it->tot_hpwl + opt_w_cof * it->tot_stnwl;  // lutong
        it->hpwl = total_hpwl;
        it->pcof = get_phi_cof1((it->tot_hpwl - last_it->tot_hpwl) / ref_dwl0);
        if(temp_iter < 1) {
            opt_phi_cof *= it->pcof;  // igkang
        }
        if(lambda2CMD == true) {
            opt_phi_cof_local = 0;
            for(int i = 0; i < N; i++) {
                opt_phi_cof_local += (cellLambdaArr[i]);
            }
            opt_phi_cof_local /= (prec)N;
        }
        it->wlen = total_wlen;
    }
    else {
        it->hpwl = last_it->hpwl;
        it->tot_hpwl = last_it->tot_hpwl;
        it->tot_stnwl = last_it->tot_stnwl;  // lutong
        it->tot_wwl = last_it->tot_wwl;      // lutong
        it->wlen.SetZero();
        it->tot_wlen = 0;
    }
    time_calc(last_it->cpu_curr, &it->cpu_curr, &it->cpu_cost);
    PrintNesterovOptStatus(iter);
    fflush(stdout);
}

/*void __lipschitz__(fpos2_t *y_st, fpos2_t *y_dst, fpos2_t *z_st,
             fpos2_t *z_dst, struct ITER *iter, int N) {
    prec yz_dis = 0;
    prec yz_dnm = 0;
    prec alpha;
    prec lc = 0;
    yz_dis = get_dis(y_st, z_st, N);
    yz_dnm = get_dis(y_dst, z_dst, N);
    lc = yz_dnm / yz_dis;
    alpha = 1.0 / lc;
    iter->lc = lc;
    iter->alpha00 = alpha;
}*/

void get_lc(struct FPOS *y_st, struct FPOS *y_dst, struct FPOS *z_st,
            struct FPOS *z_dst, struct ITER *iter, int N) {
    if(FILLER_PLACE)
        get_lc3_filler(y_st + moduleCNT, y_dst + moduleCNT, z_st + moduleCNT,
                       z_dst + moduleCNT, iter, gfiller_cnt);
    else
        get_lc3(y_st, y_dst, z_st, z_dst, iter, N);
}

void get_lc3(struct FPOS *y_st, struct FPOS *y_dst, struct FPOS *z_st,
             struct FPOS *z_dst, struct ITER *iter, int N) {
    prec yz_dis = 0;
    prec yz_dnm = 0;
    prec alpha;
    prec lc = 0;
    yz_dis = get_dis(y_st, z_st, N);
    yz_dnm = get_dis(y_dst, z_dst, N);
    lc = yz_dnm / yz_dis;
    alpha = 1.0 / lc;
    iter->lc = lc;
    iter->alpha00 = alpha;
}

void get_lc3_filler(struct FPOS *y_st, struct FPOS *y_dst, struct FPOS *z_st,
                    struct FPOS *z_dst, struct ITER *iter, int N) {
    prec yz_dis = 0;
    prec yz_dnm = 0;
    prec alpha;
    prec lc = 0;

    yz_dis = get_dis(y_st, z_st, N);
    yz_dnm = get_dis(y_dst, z_dst, N);

    lc = yz_dnm / yz_dis;
    alpha = 1.0 / lc;
    iter->lc = lc;
    iter->alpha00 = alpha;
}

void myNesterov::InitializationCostFunctionGradient(prec *sum_wgrad0,
                                                    prec *sum_pgrad0) {
    prec tmp_sum_wgrad = 0;
    prec tmp_sum_pgrad = 0;
    struct FPOS wgrad;
    struct FPOS pgrad;
    struct FPOS pgradl;
    struct CELLx *cell = NULL;

    for(int i = 0; i < N; i++) {
        cell = &gcell_st[i];
        wlen_grad(i, &wgrad);
        if(STAGE == cGP2D) {
            if(cell->flg == Macro) {
                wgrad.SetZero();
                pgrad.SetZero();
                pgradl.SetZero();
            }
            else {
                if(constraintDrivenCMD == false)
                    potn_grad_2D(i, &pgrad);
                else if(constraintDrivenCMD == true) {
                    if(lambda2CMD == false)
                        potn_grad_2D_local(i, &pgrad, &cellLambdaArr[i]);
                    else if(lambda2CMD == true) {
                        potn_grad_2D(i, &pgrad);
                        potn_grad_2D_local(i, &pgradl, &cellLambdaArr[i]);
                    }
                }
            }
        }
        else if(STAGE == mGP2D) {
            if(constraintDrivenCMD == false)
                potn_grad_2D(i, &pgrad);
            else if(constraintDrivenCMD == true) {
                if(lambda2CMD == false)
                    potn_grad_2D_local(i, &pgrad, &cellLambdaArr[i]);
                else if(lambda2CMD == true) {
                    potn_grad_2D(i, &pgrad);
                    potn_grad_2D_local(i, &pgradl, &cellLambdaArr[i]);
                }
            }
        }

        y_wdst[i] = wgrad;
        y_pdst[i] = pgrad;
        if(lambda2CMD == true)
            y_pdstl[i] = pgradl;

        tmp_sum_wgrad += fabs(wgrad.x) + fabs(wgrad.y);  
        tmp_sum_pgrad += fabs(pgrad.x) + fabs(pgrad.y);  
    }
    *sum_wgrad0 = tmp_sum_wgrad;
    *sum_pgrad0 = tmp_sum_pgrad;
}

void myNesterov::ShiftPL_SA(struct FPOS *y_st, int N) {
    if(FILLER_PLACE)
        ShiftPL_SA_sub(y_st + moduleCNT, gfiller_cnt);
    else
        ShiftPL_SA_sub(y_st, N);
}

void myNesterov::ShiftPL_SA_sub(struct FPOS *y_st, int N) {
    struct FPOS rnd;
    struct FPOS len;
    struct FPOS drnd;
    struct FPOS v;
    struct FPOS d;
    struct CELLx *cell = NULL;
    struct FPOS half_densize;
    struct FPOS center;

    prec ratio_var_pl = 0.025;

    for(int i = 0; i < N; i++) {
        cell = &gcell_st[i];
        if(STAGE == cGP3D) {
            if(cell->flg == Macro)
                continue;
        }
        half_densize = cell->half_den_size;

        center = y_st[i];

        len.x = 2.0 * ratio_var_pl * half_densize.x;
        len.y = 2.0 * ratio_var_pl * half_densize.y;

        rnd.x = (prec)rand();
        rnd.y = (prec)rand();

        drnd.x = inv_RAND_MAX * rnd.x - 0.5;
        drnd.y = inv_RAND_MAX * rnd.y - 0.5;

        d.x = drnd.x * len.x;
        d.y = drnd.y * len.y;

        v.x = center.x + d.x;
        v.y = center.y + d.y;

        center = valid_coor00(v, half_densize);

        y_st[i] = center;
    }
}

void myNesterov::UpdateAlpha(struct ITER *it) {
    if(ALPHA < maxALPHA) {
        if(it->pcof > LOWER_PCOF && it->pcof < UPPER_PCOF) {
            ALPHA *= it->pcof;
        }
        else if(it->pcof <= LOWER_PCOF) {
            ALPHA *= LOWER_PCOF;
        }
        else if(it->pcof >= UPPER_PCOF) {
            ALPHA *= UPPER_PCOF;
        }
    }
}

void myNesterov::UpdateBeta(struct ITER *it) {
    if(it->pcof > LOWER_PCOF && it->pcof < UPPER_PCOF) {
        BETA *= it->pcof;
    }
    else if(it->pcof <= LOWER_PCOF) {
        BETA *= LOWER_PCOF;
    }
    else if(it->pcof >= UPPER_PCOF) {
        BETA *= UPPER_PCOF;
    }
}

void myNesterov::PrintNesterovOptStatus(int iter) {
if(iter % 50 == 0 || iter > 2000) {
    printf("\n");
    printf("ITER: %d\n", iter);
    printf("    HPWL=%.6f\n", it->tot_hpwl);
    printf("    OVFL=%.6f\n", it->ovfl);
    printf("    HPWL=(%.6f, %.6f)\n", it->hpwl.x, it->hpwl.y);
    printf("    POTN=%.6E\n", it->potn);
    printf("    PHIC=%.6E\n", opt_phi_cof);
    printf("    GRAD=%.6E\n", it->grad);
    printf("    NuBT=%d\n", backtrack_cnt);
    printf("    CPU =%.3f  (wlen: %.4f) (grad: %.4f) (den: %.4f)\n", it->cpu_cost, it->wlcost, it->gradcost, it->dencost);
}
}

void myNesterov::PlotNesterov(int iter) {
    if(plotCellCMD) {
        cell_update(x_st, N);
        if(STAGE == mGP3D)
            plot("S1-cell-mGP3D-", iter, 1.0, 1);
        if(STAGE == mGP2D)
            plot("S2-cell-mGP2D-", iter, 1.0, 1);
        if(STAGE == cGP3D)
            plot("S4-cell-cGP3D-", iter, 1.0, 1);
        if(STAGE == cGP2D)
            plot("S5-cell-cGP2D-", iter, 1.0, 1);
    }
    if(plotDensityCMD) {
        if(STAGE == mGP3D)
            plot_bin("S1-den-mGP3D-", iter, 1.0, 0);
        if(STAGE == mGP2D)
            plot_bin("S2-den-mGP2D-", iter, 1.0, 0);
        if(STAGE == cGP3D)
            plot_bin("S4-den-cGP3D-", iter, 1.0, 0);
        if(STAGE == cGP2D)
            plot_bin("S5-den-cGP2D-", iter, 1.0, 0);
    }
    if(plotFieldCMD) {
        if(STAGE == mGP3D)
            plot_bin("S1-field-mGP3D-", iter, 1.0, 4);
        if(STAGE == mGP2D)
            plot_bin("S2-field-mGP2D-", iter, 1.0, 4);
        if(STAGE == cGP3D)
            plot_bin("S4-field-cGP3D-", iter, 1.0, 4);
        if(STAGE == cGP2D)
            plot_bin("S5-field-cGP2D-", iter, 1.0, 4);
    }
}
