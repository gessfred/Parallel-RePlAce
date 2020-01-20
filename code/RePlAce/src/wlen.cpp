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
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <omp.h>
#include <fstream>

#include "global.h"
#include "macro.h"
#include "mkl.h"
#include "opt.h"
#include "wlen.h"

int wcof_flg;
int MAX_EXP;
int NEG_MAX_EXP;
prec hpwl_mGP3D;
prec hpwl_mGP2D;
prec hpwl_cGP3D;
prec hpwl_cGP2D;
prec TSV_WEIGHT;

FPOS wcof00;
FPOS wcof00_dim1;
FPOS wcof00_org;
FPOS total_hpwl;
FPOS total_stnwl;  // lutong
FPOS total_wlen;
FPOS gp_wlen_weight;
FPOS dp_wlen_weight;
FPOS base_wcof;
FPOS wlen_cof;
FPOS wlen_cof_inv;
EXP_ST *exp_st;

void SetMAX_EXP_wlen() {
    MAX_EXP = 300;
    NEG_MAX_EXP = -300;
}

FPOS get_wlen_cof(prec ovf) {
    switch(wcof_flg) {
        case 1:
            return get_wlen_cof1(ovf);
            break;
        case 2:
            return get_wlen_cof2(ovf);
            break;
        default:
            return get_wlen_cof3(ovf);
            break;
    }
}

FPOS get_wlen_cof3(prec ovf) {
    FPOS cof;
    prec tmp = 0;
    if(ovf > 1.0) {
        cof.x = 1.0 / 10.0;
        cof.y = 1.0 / 10.0;
        cof.z = 1.0 / 10.0;
    }
    else if(ovf < 0.1) {
        cof.x = 1.0 / 0.01;
        cof.y = 1.0 / 0.01;
        cof.z = 1.0 / 0.01;
    }
    else {
        tmp = 1.0 / pow(10.0, (ovf - 0.1) * 10.0 / 3.0 - 2.0);
        cof.x = cof.y = cof.z = tmp;
    }
    return cof;
}

FPOS get_wlen_cof1(prec ovf) {
    FPOS cof;

    if(ovf > 1.0)
        cof.x = cof.y = cof.z = 1.0 / 10.0;
    else if(ovf < 0.1)
        cof.x = cof.y = cof.z = 1.0 / 0.01;
    else
        cof.x = cof.y = cof.z = 1.0 / pow(10.0, (ovf - 0.1) * 10.0 / 3.0 - 2.0);

    return cof;
}

FPOS get_wlen_cof2(prec ovf) {
    FPOS cof;
    prec tmp = 0.0;
    if(ovf > 1.0) {
        cof.x = 0.1;
        cof.y = 0.1;
        cof.z = 0.1;
    }
    else if(ovf < 0.1) {
        cof.x = 10.0;
        cof.y = 10.0;
        cof.z = 10.0;
    }
    else {
        tmp = 1.0 / pow(10.0, (ovf - 0.1) * 20 / 9.0 - 1.0);
        cof.x = cof.y = cof.z = tmp;
    }
    return cof;
}

void wlen_init(void) {
    int i = 0 /* ,cnt=exp_st_cnt */;
    /* prec interval = exp_interval ; */
    exp_st = (EXP_ST *)mkl_malloc(sizeof(EXP_ST) * exp_st_cnt, 64);
    for(i = 0; i < exp_st_cnt; i++) {
        exp_st[i].x = (prec)i * exp_interval - MAX_EXP;
        exp_st[i].val = exp(exp_st[i].x);
        if(i > 0) {
            exp_st[i - 1].y_h =
                (exp_st[i].val - exp_st[i - 1].val) / exp_interval;
        }
    }

    if(GP_DIM_ONE) {
        gp_wlen_weight.x = 1.0;
        gp_wlen_weight.y = place_backup.cnt.y / place_backup.cnt.x;
        // gp_wlen_weight.y = place.cnt.y / place.cnt.x;
        if(flg_3dic)
            gp_wlen_weight.z = TSV_WEIGHT;
        else
            gp_wlen_weight.z = 0.0;
    }
    else {
        gp_wlen_weight.x = gp_wlen_weight.y = 1.0;
        if(flg_3dic)
            gp_wlen_weight.z = TSV_WEIGHT;
        else
            gp_wlen_weight.z = 0.0;
    }

    dp_wlen_weight.x = 1.0;
    dp_wlen_weight.y = 1.0;
    if(GP_DIM_ONE)
        dp_wlen_weight.z = TSV_WEIGHT * place_backup.cnt.x / place_backup.cnt.z;
    else
        dp_wlen_weight.z = TSV_WEIGHT;
}

void wlen_init_mGP2D(void) {
    gp_wlen_weight.x = 1.0;
    gp_wlen_weight.y = 1.0;
    gp_wlen_weight.z = 0.0;

    dp_wlen_weight.x = 1.0;
    dp_wlen_weight.y = 1.0;
    dp_wlen_weight.z = 0.0;
}

void wlen_init_cGP2D(void) {
    gp_wlen_weight.x = 1.0;
    gp_wlen_weight.y = 1.0;
    gp_wlen_weight.z = 0.0;

    dp_wlen_weight.x = 1.0;
    dp_wlen_weight.y = 1.0;
    dp_wlen_weight.z = 0.0;
}

void wcof_init(FPOS bstp) {
    if(GP_DIM_ONE)
        wcof00 = wcof00_dim1;
    else
        wcof00 = wcof00_org;

    base_wcof.x = wcof00.x / (0.5 * (bstp.x + bstp.y));
    base_wcof.y = wcof00.y / (0.5 * (bstp.x + bstp.y));
    base_wcof.z = wcof00.z / (0.5 * (bstp.x + bstp.y));

    wlen_cof = fp_scal(0.1, base_wcof);
    wlen_cof_inv = fp_inv(wlen_cof);
}

prec get_wlen() {
    ////////////////////////////

    switch(WLEN_MODEL) {
        case WA:

            return get_wlen_wa();

            break;

        case LSE:

            return get_wlen_lse();

            break;
    }

    ////////////////////////////
    // return 0;
}

prec get_wlen_wa() {
    FPOS net_wlen;
    prec tot_wlen = 0;

    total_wlen.SetZero();

    for(int i = 0; i < netCNT; i++) {
        net_wlen = get_net_wlen_wa(&netInstance[i]);

        total_wlen.x += net_wlen.x;
        total_wlen.y += net_wlen.y;
        // total_wlen.z += net_wlen.z;
    }

    if(GP_DIM_ONE) {
        total_wlen.x *= place_backup.cnt.x * GP_SCAL;
        // * gp_wlen_weight.x +
        total_wlen.y *= place_backup.cnt.y * GP_SCAL;
        // * gp_wlen_weight.y +
        // total_wlen.z *= place_backup.cnt.z * GP_SCAL;
        // * gp_wlen_weight.z +
    }

    tot_wlen = total_wlen.x + total_wlen.y;  // +
    // total_wlen.z;

    return tot_wlen;
}

prec get_wlen_lse(void) {
    NET *net = NULL;
    FPOS net_wlen;
    prec tot_wlen = 0;

    total_wlen.SetZero();

    for(int i = 0; i < netCNT; i++) {
        net = &netInstance[i];
        if(net->pinCNTinObject <= 1)
            continue;
        net_wlen = get_net_wlen_lse(net);
        total_wlen.x += net_wlen.x;
        total_wlen.y += net_wlen.y;
        total_wlen.z += net_wlen.z;
    }

    if(GP_DIM_ONE) {
        total_wlen.x *= wlen_cof_inv.x * place_backup.cnt.x *
                        GP_SCAL;  // * gp_wlen_weight.x +
        total_wlen.y *= wlen_cof_inv.y * place_backup.cnt.y *
                        GP_SCAL;  // * gp_wlen_weight.y +
        total_wlen.z *= wlen_cof_inv.z * place_backup.cnt.z *
                        GP_SCAL;  // * gp_wlen_weight.z +
    }

    tot_wlen = total_wlen.x + total_wlen.y + total_wlen.z;

    return tot_wlen;
}

FPOS get_net_wlen_wa(NET *net) {
    FPOS net_wlen;
    FPOS sum_num1 = net->sum_num1;
    FPOS sum_num2 = net->sum_num2;
    FPOS sum_denom1 = net->sum_denom1;
    FPOS sum_denom2 = net->sum_denom2;

    if(net->pinCNTinObject <= 1)
        return zeroFPoint;

    net_wlen.x = sum_num1.x / sum_denom1.x - sum_num2.x / sum_denom2.x;

    net_wlen.y = sum_num1.y / sum_denom1.y - sum_num2.y / sum_denom2.y;

    return net_wlen;
}

FPOS get_net_wlen_lse(NET *net) {
    FPOS sum1, sum2;
    FPOS fp, wlen;
    PIN *pin = NULL;

    for(int i = 0; i < net->pinCNTinObject; i++) {
        pin = net->pin[i];
        fp = pin->fp;

        sum1.x += get_exp(wlen_cof.x * fp.x);
        sum1.y += get_exp(wlen_cof.y * fp.y);
        sum1.z += get_exp(wlen_cof.z * fp.z);

        sum2.x += get_exp(-1.0 * wlen_cof.x * fp.x);
        sum2.y += get_exp(-1.0 * wlen_cof.y * fp.y);
        sum2.z += get_exp(-1.0 * wlen_cof.z * fp.z);
    }

    wlen.x = log(sum1.x) + log(sum2.x);
    wlen.y = log(sum1.y) + log(sum2.y);
    wlen.z = log(sum1.z) + log(sum2.z);

    return wlen;
}

//
// this only calculate HPWL based on the stored value in NET's ure
prec GetHpwl() {
    total_hpwl.SetZero();
    total_stnwl.SetZero();

    for(int i = 0; i < netCNT; i++) {
        NET *curNet = &netInstance[i];
        if(curNet->pinCNTinObject <= 1)
            continue;

        total_hpwl.x += curNet->max_x - curNet->min_x;
        total_hpwl.y += curNet->max_y - curNet->min_y;
    }

    if(GP_DIM_ONE) {
        total_hpwl.x *= place_backup.cnt.x * GP_SCAL;
        total_hpwl.y *= place_backup.cnt.y * GP_SCAL;
    }

    return total_hpwl.x + total_hpwl.y;
}

//
// this calculate current NET's informations (min_xyz/max_xyz) & calculate HPWL
// simultaneously
prec UpdateNetAndGetHpwl() {
    FPOS pof, fp;
    total_hpwl.SetZero();

    for(int i = 0; i < netCNT; i++) {
        NET *curNet = &netInstance[i];

        curNet->min_x = curNet->terminalMin.x;
        curNet->min_y = curNet->terminalMin.y;
        curNet->min_z = curNet->terminalMin.z;

        curNet->max_x = curNet->terminalMax.x;
        curNet->max_y = curNet->terminalMax.y;
        curNet->max_z = curNet->terminalMax.z;

        // update min_xyz, max_xyz in NET's info
        for(int j = 0; j < curNet->pinCNTinObject; j++) {
            PIN *pin = curNet->pin[j];

            // only for modules
            if(pin->term) {
                continue;
            }

            MODULE *curModule = &moduleInstance[pin->moduleID];
            pin->fp.SetAdd(curModule->center,
                           curModule->pof[pin->pinIDinModule]);

            curNet->min_x = min(curNet->min_x, pin->fp.x);
            curNet->min_y = min(curNet->min_y, pin->fp.y);
            curNet->min_z = min(curNet->min_z, pin->fp.z);

            curNet->max_x = max(curNet->max_x, pin->fp.x);
            curNet->max_y = max(curNet->max_y, pin->fp.y);
            curNet->max_z = max(curNet->max_z, pin->fp.z);
        }

        if(curNet->pinCNTinObject <= 1) {
            continue;
        }

        // calculate HPWL
        total_hpwl.x += curNet->max_x - curNet->min_x;
        total_hpwl.y += curNet->max_y - curNet->min_y;
        total_hpwl.z += curNet->max_z - curNet->min_z;
    }

    return total_hpwl.x + total_hpwl.y + total_hpwl.z;
}

void wlen_grad2(int cell_idx, FPOS *grad2) {
    switch(WLEN_MODEL) {
        case LSE:
            wlen_grad2_lse(cell_idx, grad2);
            break;

        case WA:
            wlen_grad2_wa(grad2);
            break;
    }
    return;
}

void wlen_grad(int cell_idx, FPOS *grad) {
    grad->SetZero();
#ifdef NO_WLEN
    return;
#endif

    switch(WLEN_MODEL) {
        case LSE:
            wlen_grad_lse(cell_idx, grad);
            break;

        case WA:
            wlen_grad_wa(cell_idx, grad);
            break;
    }
    grad->x *= -1.0 * gp_wlen_weight.x;
    grad->y *= -1.0 * gp_wlen_weight.y;
}

void wlen_grad2_wa(FPOS *grad) {
    grad->x = 1.0;
    grad->y = 1.0;
    if(flg_3dic)
        grad->z = 1.0;
    return;
}

void wlen_grad2_lse(int cell_idx, FPOS *grad2) {
    FPOS net_grad2;
    CELLx *cell = &gcell_st[cell_idx];
    NET *net = NULL;
    PIN *pin = NULL;

    grad2->SetZero();

    for(int i = 0; i < cell->pinCNTinObject; i++) {
        pin = cell->pin[i];
        net = &netInstance[pin->netID];

        if(net->pinCNTinObject <= 1)
            continue;

        get_net_wlen_grad2_lse(net, pin, &net_grad2);

        grad2->x += net_grad2.x;
        grad2->y += net_grad2.y;
        grad2->z += net_grad2.z;
    }
}

void wlen_grad_lse(int cell_idx, FPOS *grad) {
    FPOS net_grad = zeroFPoint;
    CELLx *cell = &gcell_st[cell_idx];
    NET *net = NULL;
    PIN *pin = NULL;

    grad->SetZero();

    for(int i = 0; i < cell->pinCNTinObject; i++) {
        pin = cell->pin[i];
        net = &netInstance[pin->netID];

        if(net->pinCNTinObject <= 1)
            continue;

        get_net_wlen_grad_lse(net, pin, &net_grad);

        grad->x += net_grad.x;
        grad->y += net_grad.y;
        grad->z += net_grad.z;
    }
}

void __get_net_wlen_grad_wa__(FPOS obj, NET *net, PIN *pin, FPOS *grad) {
    FPOS grad_sum_num1 = zeroFPoint, grad_sum_num2 = zeroFPoint;
    FPOS grad_sum_denom1 = zeroFPoint, grad_sum_denom2 = zeroFPoint;
    FPOS grad1 = zeroFPoint;
    FPOS grad2 = zeroFPoint;
    FPOS e1 = pin->e1;
    FPOS e2 = pin->e2;
    POS flg1 = pin->flg1;
    POS flg2 = pin->flg2;
    FPOS sum_num1 = net->sum_num1;
    FPOS sum_num2 = net->sum_num2;
    FPOS sum_denom1 = net->sum_denom1;
    FPOS sum_denom2 = net->sum_denom2;

    if(flg1.x) {
        grad_sum_denom1.x = wlen_cof.x * e1.x;
        grad_sum_num1.x = e1.x + obj.x * grad_sum_denom1.x;
        grad1.x =
            (grad_sum_num1.x * sum_denom1.x - grad_sum_denom1.x * sum_num1.x) /
            (sum_denom1.x * sum_denom1.x);
    }

    if(flg1.y) {
        grad_sum_denom1.y = wlen_cof.y * e1.y;
        grad_sum_num1.y = e1.y + obj.y * grad_sum_denom1.y;
        grad1.y =
            (grad_sum_num1.y * sum_denom1.y - grad_sum_denom1.y * sum_num1.y) /
            (sum_denom1.y * sum_denom1.y);
    }

    if(flg2.x) {
        grad_sum_denom2.x = wlen_cof.x * e2.x;
        grad_sum_num2.x = e2.x - obj.x * grad_sum_denom2.x;
        grad2.x =
            (grad_sum_num2.x * sum_denom2.x + grad_sum_denom2.x * sum_num2.x) /
            (sum_denom2.x * sum_denom2.x);
    }

    if(flg2.y) {
        grad_sum_denom2.y = wlen_cof.y * e2.y;
        grad_sum_num2.y = e2.y - obj.y * grad_sum_denom2.y;
        grad2.y =
            (grad_sum_num2.y * sum_denom2.y + grad_sum_denom2.y * sum_num2.y) /
            (sum_denom2.y * sum_denom2.y);
    }

    grad->x = grad1.x - grad2.x;
    grad->y = grad1.y - grad2.y;
}

void __wlen_grad__(io_t* cell, field_t* nets, FPOS *grad) {
    fpos2_t net_grad;

    grad->SetZero();
    for(int i = 0; i < cell->pinCNT; i++) {
        charge_t* pin = cell->pins[i];
        field_t* net = &nets[pin->netID];
        if(net->pinCNT <= 1)
            continue;
        __dwlen__(pin->fp, net, pin, &net_grad);
        grad->x += net_grad.x;
        grad->y += net_grad.y;
    }
    grad->x *= -1.0 * gp_wlen_weight.x;
    grad->y *= -1.0 * gp_wlen_weight.y;
}

void wlen_grad_wa(int cell_idx, FPOS *grad) {
    CELLx *cell = &gcell_st[cell_idx];
    PIN *pin = NULL;
    NET *net = NULL;
    FPOS net_grad = zeroFPoint;

    grad->SetZero();

    for(int i = 0; i < cell->pinCNTinObject; i++) {
        pin = cell->pin[i];
        net = &netInstance[pin->netID];

        if(net->pinCNTinObject <= 1)
            continue;

        get_net_wlen_grad_wa(pin->fp, net, pin, &net_grad);
        grad->x += net_grad.x;
        grad->y += net_grad.y;
    }
}

void get_net_wlen_grad2_lse(NET *net, PIN *pin, FPOS *grad2) {
    POS flg1 = pin->flg1, flg2 = pin->flg2;
    FPOS e1 = pin->e1, e2 = pin->e2;
    FPOS grad2_1 = zeroFPoint, grad2_2 = zeroFPoint;
    FPOS sum_denom1 = net->sum_denom1;
    FPOS sum_denom2 = net->sum_denom2;

    if(flg1.x) {
        grad2_1.x =
            (e1.x) * (sum_denom1.x - e1.x) / (sum_denom1.x * sum_denom1.x);
    }

    if(flg2.x) {
        grad2_2.x =
            (e2.x) * (sum_denom2.x - e2.x) / (sum_denom2.x * sum_denom2.x);
    }

    grad2->x = grad2_1.x + grad2_2.x;

    if(flg1.y) {
        grad2_1.y =
            (e1.y) * (sum_denom1.y - e1.y) / (sum_denom1.y * sum_denom1.y);
    }

    if(flg2.y) {
        grad2_2.y =
            (e2.y) * (sum_denom2.y - e2.y) / (sum_denom2.y * sum_denom2.y);
    }

    grad2->y = grad2_1.y + grad2_2.y;

    if(flg_3dic) {
        if(flg1.z) {
            grad2_1.z =
                (e1.z) * (sum_denom1.z - e1.z) / (sum_denom1.z * sum_denom1.z);
        }

        if(flg2.z) {
            grad2_2.z =
                (e2.z) * (sum_denom2.z - e2.z) / (sum_denom2.z * sum_denom2.z);
        }

        grad2->z = (grad2_1.z + grad2_2.z);
    }

    return;
}

void get_net_wlen_grad_lse(NET *net, PIN *pin, FPOS *grad) {
    POS flg1 = pin->flg1, flg2 = pin->flg2;
    FPOS grad1 = zeroFPoint, grad2 = zeroFPoint;
    FPOS e1 = pin->e1, e2 = pin->e2;
    FPOS sum_denom1 = net->sum_denom1;
    FPOS sum_denom2 = net->sum_denom2;

    if(flg1.x) {
        grad1.x = e1.x / sum_denom1.x;
    }

    if(flg2.x) {
        grad2.x = e2.x / sum_denom2.x;
    }

    grad->x = grad1.x - grad2.x;

    if(flg1.y) {
        grad1.y = e1.y / sum_denom1.y;
    }

    if(flg2.y) {
        grad2.y = e2.y / sum_denom2.y;
    }

    grad->y = grad1.y - grad2.y;

    if(flg_3dic) {
        if(flg1.z) {
            grad1.z = e1.z / sum_denom1.z;
        }

        if(flg2.z) {
            grad2.z = e2.z / sum_denom2.z;
        }

        grad->z = (grad1.z - grad2.z);
    }

    return;
}

void __dwlen__(fpos2_t obj, field_t *net, charge_t *pin, fpos2_t *grad) {
    fpos2_t grad_sum_num1 , grad_sum_num2;
    fpos2_t grad_sum_denom1 , grad_sum_denom2;
    fpos2_t grad1 = {0, 0};
    fpos2_t grad2 = {0, 0};
    fpos2_t e1 = pin->e1;
    fpos2_t e2 = pin->e2;
    fpos2_t sum_num1 = net->sum_num1;
    fpos2_t sum_num2 = net->sum_num2;
    fpos2_t sum_denom1 = net->sum_denom1;
    fpos2_t sum_denom2 = net->sum_denom2;

    if(pin->meta[5]) {
        grad_sum_denom1.x = wlen_cof.x * e1.x;
        grad_sum_num1.x = e1.x + obj.x * grad_sum_denom1.x;
        grad1.x =
            (grad_sum_num1.x * sum_denom1.x - grad_sum_denom1.x * sum_num1.x) /
            (sum_denom1.x * sum_denom1.x);
    }

    if(pin->meta[6]) {
        grad_sum_denom1.y = wlen_cof.y * e1.y;
        grad_sum_num1.y = e1.y + obj.y * grad_sum_denom1.y;
        grad1.y =
            (grad_sum_num1.y * sum_denom1.y - grad_sum_denom1.y * sum_num1.y) /
            (sum_denom1.y * sum_denom1.y);
    }

    if(pin->meta[7]) {
        grad_sum_denom2.x = wlen_cof.x * e2.x;
        grad_sum_num2.x = e2.x - obj.x * grad_sum_denom2.x;
        grad2.x =
            (grad_sum_num2.x * sum_denom2.x + grad_sum_denom2.x * sum_num2.x) /
            (sum_denom2.x * sum_denom2.x);
    }

    if(pin->meta[8]) {
        grad_sum_denom2.y = wlen_cof.y * e2.y;
        grad_sum_num2.y = e2.y - obj.y * grad_sum_denom2.y;
        grad2.y =
            (grad_sum_num2.y * sum_denom2.y + grad_sum_denom2.y * sum_num2.y) /
            (sum_denom2.y * sum_denom2.y);
    }
    grad->x = grad1.x - grad2.x;
    grad->y = grad1.y - grad2.y;
}

inline void get_net_wlen_grad_wa(FPOS obj, NET *net, PIN *pin, FPOS *grad) {
    FPOS grad_sum_num1 = zeroFPoint, grad_sum_num2 = zeroFPoint;
    FPOS grad_sum_denom1 = zeroFPoint, grad_sum_denom2 = zeroFPoint;
    FPOS grad1 = zeroFPoint;
    FPOS grad2 = zeroFPoint;
    FPOS e1 = pin->e1;
    FPOS e2 = pin->e2;
    POS flg1 = pin->flg1;
    POS flg2 = pin->flg2;
    FPOS sum_num1 = net->sum_num1;
    FPOS sum_num2 = net->sum_num2;
    FPOS sum_denom1 = net->sum_denom1;
    FPOS sum_denom2 = net->sum_denom2;

    if(flg1.x) {
        grad_sum_denom1.x = wlen_cof.x * e1.x;
        grad_sum_num1.x = e1.x + obj.x * grad_sum_denom1.x;
        grad1.x =
            (grad_sum_num1.x * sum_denom1.x - grad_sum_denom1.x * sum_num1.x) /
            (sum_denom1.x * sum_denom1.x);
    }

    if(flg1.y) {
        grad_sum_denom1.y = wlen_cof.y * e1.y;
        grad_sum_num1.y = e1.y + obj.y * grad_sum_denom1.y;
        grad1.y =
            (grad_sum_num1.y * sum_denom1.y - grad_sum_denom1.y * sum_num1.y) /
            (sum_denom1.y * sum_denom1.y);
    }

    if(flg2.x) {
        grad_sum_denom2.x = wlen_cof.x * e2.x;
        grad_sum_num2.x = e2.x - obj.x * grad_sum_denom2.x;
        grad2.x =
            (grad_sum_num2.x * sum_denom2.x + grad_sum_denom2.x * sum_num2.x) /
            (sum_denom2.x * sum_denom2.x);
    }

    if(flg2.y) {
        grad_sum_denom2.y = wlen_cof.y * e2.y;
        grad_sum_num2.y = e2.y - obj.y * grad_sum_denom2.y;
        grad2.y =
            (grad_sum_num2.y * sum_denom2.y + grad_sum_denom2.y * sum_num2.y) /
            (sum_denom2.y * sum_denom2.y);
    }

    grad->x = grad1.x - grad2.x;
    grad->y = grad1.y - grad2.y;
}

void net_update_init(void) {
    for(int i = 0; i < netCNT; i++) {
        NET *net = &netInstance[i];

        net->terminalMin.Set(place.end);
        net->terminalMax.Set(place.org);

        bool first_term = true;

        for(int j = 0; j < net->pinCNTinObject; j++) {
            PIN *pin = net->pin[j];
            if(pin->term) {
                if(first_term) {
                    first_term = false;
                    net->terminalMin.Set(pin->fp);
                    net->terminalMax.Set(pin->fp);
                }
                else {
                    net->terminalMin.Min(pin->fp);
                    net->terminalMax.Max(pin->fp);
                }
            }
        }
    }
}

void net_update(FPOS *st) {
    switch(WLEN_MODEL) {
        case LSE:
            return net_update_lse(st);
            break;
        case WA:
            return net_update_wa(st);
            break;
    }
}

void net_update_lse(FPOS *st) {
    int i = 0, j = 0;
    CELLx *cell = NULL;
    NET *net = NULL;
    PIN *pin = NULL;
    MODULE *curModule = NULL;
    FPOS fp = zeroFPoint, pof = zeroFPoint, center = zeroFPoint;
    FPOS sum_denom1, sum_denom2;
    prec exp_val = 0;
    prec min_x = 0, min_y = 0;
    prec max_x = 0, max_y = 0;
    prec exp_min_x = 0, exp_min_y = 0;
    prec exp_max_x = 0, exp_max_y = 0;

    for(i = 0; i < gcell_cnt; i++) {
        cell = &gcell_st[i];

        cell->center = st[i];

        cell->den_pmin.x = cell->center.x - cell->half_den_size.x;
        cell->den_pmin.y = cell->center.y - cell->half_den_size.y;
        cell->den_pmax.x = cell->center.x + cell->half_den_size.x;
        cell->den_pmax.y = cell->center.y + cell->half_den_size.y;
    }

    for(i = 0; i < netCNT; i++) {
        net = &netInstance[i];

        net->min_x = net->terminalMin.x;
        net->min_y = net->terminalMin.y;
        // net->min_z = net->terminalMin.z;

        net->max_x = net->terminalMax.x;
        net->max_y = net->terminalMax.y;
        // net->max_z = net->terminalMax.z;

        for(j = 0; j < net->pinCNTinObject; j++) {
            pin = net->pin[j];

            if(!pin->term) {
                curModule = &moduleInstance[pin->moduleID];
                pof = curModule->pof[pin->pinIDinModule];
                center = st[pin->moduleID];
                fp.x = center.x + pof.x;
                fp.y = center.y + pof.y;
                // fp.z = center.z + pof.z ;
                pin->fp = fp;

                net->min_x = min(net->min_x, fp.x);
                net->min_y = min(net->min_y, fp.y);
                // net->min_z = min ( net->min_z , fp.z ) ;

                net->max_x = max(net->max_x, fp.x);
                net->max_y = max(net->max_y, fp.y);
                // net->max_z = max ( net->max_z , fp.z ) ;
            }
            else {
                continue;
            }
        }

        min_x = net->min_x;
        min_y = net->min_y;
        // min_z = net->min_z ;

        max_x = net->max_x;
        max_y = net->max_y;
        // max_z = net->max_z;

        sum_denom1 = zeroFPoint;
        sum_denom2 = zeroFPoint;

        for(j = 0; j < net->pinCNTinObject; j++) {
            pin = net->pin[j];

            if(!pin->term) {
#ifdef CELL_CENTER_WLEN_GRAD
                fp = st[pin->moduleID];
#else
                fp = pin->fp;
#endif
            }
            else {
                fp = pin->fp;
            }

            exp_max_x = (fp.x - max_x) * wlen_cof.x;

            if(fabs(exp_max_x) < MAX_EXP) {
                exp_val = get_exp(exp_max_x);
                sum_denom1.x += exp_val;
                pin->flg1.x = 1;
                pin->e1.x = exp_val;
            }
            else {
                exp_val = 0;
                pin->flg1.x = 0;
            }

            exp_min_x = (min_x - fp.x) * wlen_cof.x;

            if(fabs(exp_min_x) < MAX_EXP) {
                exp_val = get_exp(exp_min_x);
                sum_denom2.x += exp_val;
                pin->flg2.x = 1;
                pin->e2.x = exp_val;
            }
            else {
                pin->flg2.x = 0;
            }

            exp_max_y = (fp.y - max_y) * wlen_cof.y;

            if(fabs(exp_max_y) < MAX_EXP) {
                exp_val = get_exp(exp_max_y);
                sum_denom1.y += exp_val;
                pin->flg1.y = 1;
                pin->e1.y = exp_val;
            }
            else {
                pin->flg1.y = 0;
            }

            exp_min_y = (min_y - fp.y) * wlen_cof.y;

            if(fabs(exp_min_y) < MAX_EXP) {
                exp_val = get_exp(exp_min_y);
                sum_denom2.y += exp_val;
                pin->flg2.y = 1;
                pin->e2.y = exp_val;
            }
            else {
                pin->flg2.y = 0;
            }
        net->sum_denom1 = sum_denom1;
        net->sum_denom2 = sum_denom2;
    }
    }
}

prec net_update_hpwl_mac(void) {
    int i = 0, j = 0;
    prec hpwl = 0;
    NET *net = NULL;
    PIN *pin = NULL;
    MODULE *curModule = NULL;
    FPOS fp = zeroFPoint, pof = zeroFPoint, p0 = zeroFPoint;

    total_hpwl = zeroFPoint;

    for(i = 0; i < netCNT; i++) {
        net = &netInstance[i];

        net->min_x = net->terminalMin.x;
        net->min_y = net->terminalMin.y;
        net->min_z = net->terminalMin.z;

        net->max_x = net->terminalMax.x;
        net->max_y = net->terminalMax.y;
        net->max_z = net->terminalMax.z;

        for(j = 0; j < net->pinCNTinObject; j++) {
            pin = net->pin[j];
            if(!pin->term) {
                curModule = &moduleInstance[pin->moduleID];
                // cell = & gcell_st [pin->moduleID];
                p0 = curModule->center;
                pof = curModule->pof[pin->pinIDinModule];
                fp.x = p0.x + pof.x;
                fp.y = p0.y + pof.y;
                fp.z = p0.z + pof.z;
                pin->fp = fp;

                net->min_x = min(net->min_x, fp.x);
                net->min_y = min(net->min_y, fp.y);
                net->min_z = min(net->min_z, fp.z);

                net->max_x = max(net->max_x, fp.x);
                net->max_y = max(net->max_y, fp.y);
                net->max_z = max(net->max_z, fp.z);
            }
        }

        if(net->pinCNTinObject <= 1)
            continue;

        total_hpwl.x += (net->max_x - net->min_x);
        total_hpwl.y += (net->max_y - net->min_y);
        total_hpwl.z += (net->max_z - net->min_z);
    }

    /* total_hpwl_xy = total_hpwl.x + total_hpwl.y; */
    //  total_hpwl_xyz = total_hpwl.x + total_hpwl.y + total_hpwl.z;

    hpwl = total_hpwl.x * dp_wlen_weight.x + total_hpwl.y * dp_wlen_weight.y +
           total_hpwl.z * dp_wlen_weight.z;

    return hpwl;
}

void __wirelength__(field_t* nets, cell_t* cells, FPOS *st, fpos2_t** poff, size_t* wrk_ld, double* time) {
    for(int b = 0; b < gcell_cnt; b++) {
        CELL* cell = &gcell_st[b];
        cell->center = st[b];
        cells[b].den_pmin.x = cell->den_pmin.x = cell->center.x - cell->half_den_size.x;
        cells[b].den_pmin.y = cell->den_pmin.y = cell->center.y - cell->half_den_size.y;
        cells[b].den_pmax.x = cell->den_pmax.x = cell->center.x + cell->half_den_size.x;
        cells[b].den_pmax.y = cell->den_pmax.y = cell->center.y + cell->half_den_size.y;
    }
    for(int i = 0; i < netCNT; i++) {
        NET* net = &netInstance[i];
        field_t* net2 = &nets[i];
        net2->min.x = net->terminalMin.x;
        net2->min.y = net->terminalMin.y;
        net2->max.x = net->terminalMax.x;
        net2->max.y = net->terminalMax.y;
    }
    auto start = std::chrono::system_clock::now();
    #pragma omp parallel num_threads(numThread)
    {
        int tid=omp_get_thread_num();
        int start = ((tid >= 1) ? wrk_ld[tid - 1] : 0);
        int end = (tid < numThread-1) ? wrk_ld[tid] : netCNT;
        //#pragma omp parallel for num_threads(numThread)
        //for(int i = 0; i < netCNT; ++i) {
        for(int i = start; i < end; ++i) {
            field_t* net = &nets[i];
            //terminal min
            //net->min = net->terMin;
            //net->max = net->terMax;
            for(int j = 0; j < net->pinCNT; j++) {
                charge_t* pin = &net->pin[j];
                if(!pin->meta[0]) {
                    fpos2_t pof = poff[pin->moduleID][pin->idx];
                    FPOS center = st[pin->moduleID];
                    fpos2_t fp;
                    fp.x = center.x + pof.x;
                    fp.y = center.y + pof.y;
                    pin->fp = fp;
                    net->min.x = min(net->min.x, fp.x);
                    net->min.y = min(net->min.y, fp.y);
                    net->max.x = max(net->max.x, fp.x);
                    net->max.y = max(net->max.y, fp.y);
                }
            }
            prec min_x = net->min.x;
            prec min_y = net->min.y;
            prec max_x = net->max.x;
            prec max_y = net->max.y;
            fpos2_t sum_num1 = {0, 0}, sum_num2 = {0, 0};
            fpos2_t sum_denom1 = {0, 0}, sum_denom2 = {0, 0};
            for(int j = 0; j < net->pinCNT; j++) {
                charge_t* pin = &net->pin[j];
                fpos2_t fp = pin->fp;
                prec exp_max_x = (fp.x - max_x) * wlen_cof.x;
                prec exp_min_x = (min_x - fp.x) * wlen_cof.x;
                prec exp_max_y = (fp.y - max_y) * wlen_cof.y;
                prec exp_min_y = (min_y - fp.y) * wlen_cof.y;
                float e1x = get_exp(exp_max_x);
                float e2x = get_exp(exp_min_x);
                float e1y = get_exp(exp_max_y);
                float e2y = get_exp(exp_min_y);
                if(exp_max_x > NEG_MAX_EXP) {
                    pin->e1.x = e1x;
                    sum_num1.x += fp.x * e1x;
                    sum_denom1.x += e1x;
                    pin->meta[5] = 1;
                } else pin->meta[5] = 0;
                //TODO compute bus bandwidth
                if(exp_min_x > NEG_MAX_EXP) {
                    pin->e2.x = e2x;
                    sum_num2.x += fp.x * e2x;
                    sum_denom2.x += e2x;
                    pin->meta[7] = 1;
                }else pin->meta[7] = 0;
                
                if(exp_max_y > NEG_MAX_EXP) {
                    pin->e1.y = e1y;
                    sum_num1.y += fp.y * e1y;
                    sum_denom1.y += e1y;
                    pin->meta[6] =1;
                } else pin->meta[6] = 0;
                if(exp_min_y > NEG_MAX_EXP) {
                    pin->e2.y = e2y;
                    sum_num2.y += fp.y * e2y;
                    sum_denom2.y += e2y;
                    pin->meta[8] = 1;
                } else  pin->meta[8] = 0;
            }
            net->sum_num1 = sum_num1;
            net->sum_num2 = sum_num2;
            net->sum_denom1 = sum_denom1;
            net->sum_denom2 = sum_denom2;
        }
    }
    std::chrono::duration<double> diff = std::chrono::system_clock::now()-start;
    *time = diff.count();
}


void net_update_wa(FPOS *st) {

    int i = 0 ;

    bool timeon =false;
    double time = 0.0f;
    if(timeon) time_start(&time); 
    for(int b = 0; b < gcell_cnt; b++) {
        CELL* cell = &gcell_st[b];
        cell->center = st[b];
        cell->den_pmin.x = cell->center.x - cell->half_den_size.x;
        cell->den_pmin.y = cell->center.y - cell->half_den_size.y;
        cell->den_pmax.x = cell->center.x + cell->half_den_size.x;
        cell->den_pmax.y = cell->center.y + cell->half_den_size.y;
    }

    if(timeon) {time_end(&time); cout << "parallelTime : " << time << endl; };
    //#pragma omp parallel for num_threads(8)
    for(i = 0; i < netCNT; i++) {
        NET* net = &netInstance[i];
        net->min_x = net->terminalMin.x;
        net->min_y = net->terminalMin.y;
        net->max_x = net->terminalMax.x;
        net->max_y = net->terminalMax.y;
        for(int j = 0; j < net->pinCNTinObject; j++) {
            PIN* pin = net->pin[j];
            if(!pin->term) {
                MODULE* curModule = &moduleInstance[pin->moduleID];
                FPOS pof = curModule->pof[pin->pinIDinModule];
                FPOS center = st[pin->moduleID];
                FPOS fp;
                fp.x = center.x + pof.x;
                fp.y = center.y + pof.y;
                pin->fp = fp;

                net->min_x = min(net->min_x, fp.x);
                net->min_y = min(net->min_y, fp.y);
                net->max_x = max(net->max_x, fp.x);
                net->max_y = max(net->max_y, fp.y);
            }
            else {
                continue;
            }
        }
        prec min_x = net->min_x;
        prec min_y = net->min_y;
        prec max_x = net->max_x;
        prec max_y = net->max_y;

        FPOS sum_num1, sum_num2;
        FPOS sum_denom1, sum_denom2;
        for(int j = 0; j < net->pinCNTinObject; j++) {
            PIN* pin = net->pin[j];
            FPOS fp = pin->fp;
            prec exp_max_x = (fp.x - max_x) * wlen_cof.x;
            prec exp_min_x = (min_x - fp.x) * wlen_cof.x;
            prec exp_max_y = (fp.y - max_y) * wlen_cof.y;
            prec exp_min_y = (min_y - fp.y) * wlen_cof.y;

            if(exp_max_x > NEG_MAX_EXP) {
                pin->e1.x = get_exp(exp_max_x);
                sum_num1.x += fp.x * pin->e1.x;
                sum_denom1.x += pin->e1.x;
                pin->flg1.x = 1;
            }
            else {
                pin->flg1.x = 0;
            }

            if(exp_min_x > NEG_MAX_EXP) {
                pin->e2.x = get_exp(exp_min_x);
                sum_num2.x += fp.x * pin->e2.x;
                sum_denom2.x += pin->e2.x;
                pin->flg2.x = 1;
            }
            else {
                pin->flg2.x = 0;
            }

            if(exp_max_y > NEG_MAX_EXP) {
                pin->e1.y = get_exp(exp_max_y);
                sum_num1.y += fp.y * pin->e1.y;
                sum_denom1.y += pin->e1.y;
                pin->flg1.y = 1;
            }
            else {
                pin->flg1.y = 0;
            }

            if(exp_min_y > NEG_MAX_EXP) {
                pin->e2.y = get_exp(exp_min_y);
                sum_num2.y += fp.y * pin->e2.y;
                sum_denom2.y += pin->e2.y;
                pin->flg2.y = 1;
            }
            else {
                pin->flg2.y = 0;
            }
        }

        net->sum_num1 = sum_num1;
        net->sum_num2 = sum_num2;
        net->sum_denom1 = sum_denom1;
        net->sum_denom2 = sum_denom2;
    }
}

prec get_mac_hpwl(int idx) {
    MODULE *mac = macro_st[idx];
    PIN *pin = NULL;
    NET *net = NULL;
    int i = 0, j = 0;
    int moduleID = mac->idx;
    CELLx *cell = &gcell_st[moduleID];
    FPOS fp = zeroFPoint;

    mac_hpwl = zeroFPoint;
    mac_hpwl_xyz = 0;

    for(i = 0; i < cell->pinCNTinObject; i++) {
        pin = cell->pin[i];
        net = &netInstance[pin->netID];

        net->min_x = net->terminalMin.x;
        net->min_y = net->terminalMin.y;
        net->max_x = net->terminalMax.x;
        net->max_y = net->terminalMax.y;

        for(j = 0; j < net->pinCNTinObject; j++) {
            fp = net->pin[j]->fp;
            net->min_x = min(net->min_x, fp.x);
            net->min_y = min(net->min_y, fp.y);
            net->max_x = max(net->max_x, fp.x);
            net->max_y = max(net->max_y, fp.y);
        }

        if(net->pinCNTinObject <= 1)
            continue;

        mac_hpwl.x += (net->max_x - net->min_x);
        mac_hpwl.y += (net->max_y - net->min_y);
    }

    mac_hpwl_xyz =
        mac_hpwl.x * dp_wlen_weight.x + mac_hpwl.y * dp_wlen_weight.y;

    return mac_hpwl_xyz;
}

// update net->pin2->fp as updated version (modules' center + modules' offset)
void update_pin2(void) {
//    cout << "called update_pin2" << endl;
    for(int i = 0; i < netCNT; i++) {
        NET* net = &netInstance[i];

//        cout << net->pinCNTinObject2 << endl;
        for(int j = 0; j < net->pinCNTinObject2; j++) {
            PIN* pin = net->pin2[j];

            // if pin is terminal pin, then skip this procedure.
            if( pin -> term ) {
                continue;
            }

            MODULE* curModule = &moduleInstance[pin->moduleID];
            pin->fp.SetAdd( curModule->center, curModule->pof[pin->pinIDinModule] );
       
//            if( i == 0 ) {
//                cout << j << endl;
//                curModule->center.Dump("curModule->center");
//                curModule->pof[pin->pinIDinModule].Dump("pof");
//                pin->fp.Dump("Added pin->fp");
//            }

            // original code 
            // FPOS pof = curModule->pof[pin->pinIDinModule];
            // FPOS center = curModule->center;

            // fp = pin->fp;
            // fp.x = center.x + pof.x;
            // fp.y = center.y + pof.y;
            // fp.z = center.z + pof.z;
            // pin->fp = fp;
        }
    }
}

prec get_modu_hpwl(void) {
    int i = 0, j = 0;
    NET *net = NULL;
    FPOS fp;
    PIN *pin = NULL;

    tot_HPWL = tx_HPWL = ty_HPWL = tz_HPWL = 0;

    for(i = 0; i < netCNT; i++) {
        net = &netInstance[i];

        net->min_x = net->terminalMin.x;
        net->min_y = net->terminalMin.y;
        // net->min_z = net->terminalMin.z;

        net->max_x = net->terminalMax.x;
        net->max_y = net->terminalMax.y;
        // net->max_z = net->terminalMax.z;

        for(j = 0; j < net->pinCNTinObject2; j++) {
            pin = net->pin2[j];

            if(!pin->term) {
                fp = pin->fp;

                net->min_x = min(net->min_x, fp.x);
                net->min_y = min(net->min_y, fp.y);
                // net->min_z = min ( net->min_z , fp.z ) ;

                net->max_x = max(net->max_x, fp.x);
                net->max_y = max(net->max_y, fp.y);
                // net->max_z = max ( net->max_z , fp.z ) ;
            }
        }

        if(net->pinCNTinObject2 <= 1)
            continue;

        tx_HPWL += net->max_x - net->min_x;
        ty_HPWL += net->max_y - net->min_y;
        // tz_HPWL += net->max_z - net->min_z;
    }

    tot_HPWL = tx_HPWL + ty_HPWL;  // + tz_HPWL;

    return tot_HPWL;
}

int HPWL_count() {
    NET *curNet = NULL;

    tot_HPWL = 0;
    tx_HPWL = 0;
    ty_HPWL = 0;

    for(int i = 0; i < netCNT; i++) {
        curNet = &netInstance[i];
        tx_HPWL += (curNet->max_x - curNet->min_x);
        ty_HPWL += (curNet->max_y - curNet->min_y);
        // printf("%.16lf %.16lf\n",(curNet->max_x - curNet->min_x),
        //        (curNet->max_y - curNet->min_y));

        if(tx_HPWL < 0 || ty_HPWL < 0) {
            printf("ERROR\n");
            cout << tx_HPWL << " " << ty_HPWL << endl;
            g_rrr++;
            exit(1);
        }
    }
    // exit(0);
    tot_HPWL = tx_HPWL + ty_HPWL;
    return 0;
}
