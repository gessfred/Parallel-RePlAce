#include "parallel.hpp"

void update_wirelength(circuit_t* circuit, FPOS *st) {
    size_t* wrk_ld = circuit->constPinsPerNet;
    net_t* nets = circuit->nets;
    cell_den_t* cells = circuit->rects;
    for(int b = 0; b < gcell_cnt; b++) {
        CELL* cell = &gcell_st[b];
        cell->center = st[b];
        cells[b].min.x = cell->den_pmin.x = cell->center.x - cell->half_den_size.x;
        cells[b].min.y = cell->den_pmin.y = cell->center.y - cell->half_den_size.y;
        cells[b].max.x = cell->den_pmax.x = cell->center.x + cell->half_den_size.x;
        cells[b].max.y = cell->den_pmax.y = cell->center.y + cell->half_den_size.y;
    }
    for(int i = 0; i < netCNT; i++) {
        NET* net = &netInstance[i];
        net_t* net2 = &nets[i];
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
            net_t* net = &nets[i];
            //terminal min
            //net->min = net->terMin;
            //net->max = net->terMax;
            for(int j = 0; j < net->pinCNT; j++) {
                pin_t* pin = &net->pinArray[j];
                if(!pin->meta[0]) {
                    fpos2_t pof = circuit->pinOffsets[pin->moduleID][pin->pinID];
                    FPOS center = st[pin->moduleID];
                    fpos2_t coord;
                    coord.x = center.x + pof.x;
                    coord.y = center.y + pof.y;
                    pin->coord = coord;
                    net->min.x = min(net->min.x, coord.x);
                    net->min.y = min(net->min.y, coord.y);
                    net->max.x = max(net->max.x, coord.x);
                    net->max.y = max(net->max.y, coord.y);
                }
            }
            prec min_x = net->min.x;
            prec min_y = net->min.y;
            prec max_x = net->max.x;
            prec max_y = net->max.y;
            fpos2_t sum_num1 = {0, 0}, sum_num2 = {0, 0};
            fpos2_t sum_denom1 = {0, 0}, sum_denom2 = {0, 0};
            for(int j = 0; j < net->pinCNT; j++) {
                pin_t* pin = &net->pinArray[j];
                fpos2_t fp = pin->coord;
                prec exp_max_x = (fp.x - max_x) * wlen_cof.x;
                prec exp_min_x = (min_x - fp.x) * wlen_cof.x;
                prec exp_max_y = (fp.y - max_y) * wlen_cof.y;
                prec exp_min_y = (min_y - fp.y) * wlen_cof.y;
                float e1x = std::exp(exp_max_x);
                float e2x = std::exp(exp_min_x);
                float e1y = std::exp(exp_max_y);
                float e2y = std::exp(exp_min_y);
                if(exp_max_x > NEG_MAX_EXP) {
                    pin->expL.x = e1x;
                    sum_num1.x += fp.x * e1x;
                    sum_denom1.x += e1x;
                    pin->meta[5] = 1;
                } else pin->meta[5] = 0;
                //TODO compute bus bandwidth
                if(exp_min_x > NEG_MAX_EXP) {
                    pin->expR.x = e2x;
                    sum_num2.x += fp.x * e2x;
                    sum_denom2.x += e2x;
                    pin->meta[7] = 1;
                }else pin->meta[7] = 0;
                
                if(exp_max_y > NEG_MAX_EXP) {
                    pin->expL.y = e1y;
                    sum_num1.y += fp.y * e1y;
                    sum_denom1.y += e1y;
                    pin->meta[6] =1;
                } else pin->meta[6] = 0;
                if(exp_min_y > NEG_MAX_EXP) {
                    pin->expR.y = e2y;
                    sum_num2.y += fp.y * e2y;
                    sum_denom2.y += e2y;
                    pin->meta[8] = 1;
                } else  pin->meta[8] = 0;
            }
            net->sumNumL = sum_num1;
            net->sumNumR = sum_num2;
            net->sumDenomL = sum_denom1;
            net->sumDenomR = sum_denom2;
        }
    }
    std::chrono::duration<double> diff = std::chrono::system_clock::now()-start;
    //*time = diff.count();
}

//area share of a cell
void put_potential_gradient(cell_den_t* cells, bin_t* bins, FPOS* grad, int i) {
        grad->SetZero();
    cell_den_t cell = cells[i];
    POS b0, b1;
    
    TIER *tier = &tier_st[0];

    b0.x = cell.binStart.x;
    b0.y = cell.binStart.y;
    b1.x = cell.binEnd.x;
    b1.y = cell.binEnd.y;

    bin_t *bpx = NULL, *bpy = NULL;
    int x = 0, y = 0;
    int idx = b0.x * tier->dim_bin.y + b0.y;

    for(x = b0.x, bpx = &bins[idx]; x <= b1.x; x++, bpx += tier->dim_bin.y) {

        prec max_x = min(bpx->max.x, cell.max.x);
        prec min_x = max(bpx->min.x, cell.min.x);

        for(y = b0.y, bpy = bpx; y <= b1.y; y++, bpy++) {
            prec max_y = min(bpy->max.y, cell.max.y);
            prec min_y = max(bpy->min.y, cell.min.y);
            prec area_share = (max_x - min_x) * (max_y - min_y) * cell.scale;
            grad->x += area_share * bpy->field.x;
            grad->y += area_share * bpy->field.y;
        }
    }
}

void get_wirelength_grad(fpos2_t obj, net_t *net, pin_t *pin, fpos2_t *grad) {
    fpos2_t grad_sum_num1 , grad_sum_num2;
    fpos2_t grad_sum_denom1 , grad_sum_denom2;
    fpos2_t grad1 = {0, 0};
    fpos2_t grad2 = {0, 0};
    fpos2_t e1 = pin->expL;
    fpos2_t e2 = pin->expR;
    fpos2_t sum_num1 = net->sumNumL;
    fpos2_t sum_num2 = net->sumNumR;
    fpos2_t sum_denom1 = net->sumDenomL;
    fpos2_t sum_denom2 = net->sumDenomR;

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

void put_wirelength_grad(cell_phy_t* cell, net_t* nets, FPOS *grad) {
    fpos2_t net_grad;

    grad->SetZero();
    for(int i = 0; i < cell->pinCNT; i++) {
        pin_t* pin = cell->pinArrayPtr[i];
        net_t* net = &nets[pin->netID];
        if(net->pinCNT <= 1)
            continue;
        get_wirelength_grad(pin->coord, net, pin, &net_grad);
        grad->x += net_grad.x;
        grad->y += net_grad.y;
    }
    grad->x *= -1.0 * gp_wlen_weight.x;
    grad->y *= -1.0 * gp_wlen_weight.y;
}

//filler=true
void update_fill_gradient(FPOS *dst, FPOS *wdst,
                          FPOS *pdst, FPOS *pdstl, int N,
                          prec *cellLambdaArr, bool onlyPreCon, circuit_t* circuit) {
    bin_t* bins = circuit->bins;
    cell_phy_t* ios = circuit->cells;
    net_t* nets = circuit->nets;
    cell_den_t* cels = circuit->rects;
    struct FPOS wgrad, pgrad, pgradl, wpre, charge_dpre, pre;
    #pragma omp parallel for num_threads(numThread)
    for(int i = moduleCNT; i < N; i++) {
        if(cels[i].type != FillerCell) continue;
        if(STAGE == cGP2D) {
            if(constraintDrivenCMD == false) put_potential_gradient(cels, bins, &pgrad, i);
        }
        wdst[i].SetZero();
        pdst[i] = pgrad;
        dst[i].x = opt_phi_cof * pgrad.x;
        dst[i].y = opt_phi_cof * pgrad.y;
    }
    for(int i = moduleCNT; i < N; i++) {
        potn_pre(i, &charge_dpre);
        if(onlyPreCon) pre = charge_dpre;
        else { //?
          pre.x = opt_phi_cof * charge_dpre.x;
          pre.y = opt_phi_cof * charge_dpre.y;
        }
        if(pre.x < MIN_PRE) pre.x = MIN_PRE;
        if(pre.y < MIN_PRE) pre.y = MIN_PRE;
        dst[i].x /= pre.x;
        dst[i].y /= pre.y;
    }
}


void update_gradient(FPOS *dst, FPOS *wdst,
                          FPOS *pdst, FPOS *pdstl, int N,
                          prec *cellLambdaArr, bool onlyPreCon, circuit_t* circuit) {
    //double dampParam = circuit->params->dampParam;
    //double MIN_PRE = circuit->params->MIN_PRE;
    //double opt_phi_cof = circuit->params->opt_phi_cof;
    bin_t* bins = circuit->bins;
    cell_phy_t* ios = circuit->cells;
    net_t* nets = circuit->nets;
    cell_den_t* cels = circuit->rects;
    size_t* W1 = circuit->constPinsPerCell;

    struct FPOS wpre, charge_dpre, pre;
    auto t0 = std::chrono::steady_clock::now();
    #pragma omp parallel num_threads(numThread)
    {
        int tid = omp_get_thread_num();
        int start = ((tid >= 1) ? W1[tid-1] : 0);
        int end = (tid < numThread-1) ? W1[tid] : N;
        for(int i = start; i < end; ++i) {
            FPOS wgrad;
            if(STAGE == cGP2D && cels[i].type == Macro) wgrad.SetZero();
            else put_wirelength_grad(&ios[i], nets, &wgrad);
            wdst[i] = wgrad;
            dst[i].x = wgrad.x;
            dst[i].y = wgrad.y;
        }
    }
    profile.wgrad += time_since(t0);
    //time_start(time);
    t0 = std::chrono::steady_clock::now();
    #pragma omp parallel for num_threads(numThread)
    for(int i = 0; i < N; i++) {
        FPOS pgrad;
        if(STAGE == cGP2D) { //look to delete this test
            if(cels[i].type == Macro) pgrad.SetZero();
            else put_potential_gradient(cels, bins, &pgrad, i);
        }
        pdst[i] = pgrad;
        dst[i].x += opt_phi_cof * pgrad.x;
        dst[i].y += opt_phi_cof * pgrad.y;
    }
    profile.pgrad += time_since(t0);
    //time_end(time);
    t0 = std::chrono::steady_clock::now();
    for(int i = 0; i < N; i++) {
        cellLambdaArr[i] *= dampParam;
        wlen_pre(i, &wpre);
        potn_pre(i, &charge_dpre);
        if(onlyPreCon) pre = charge_dpre;
        else { //?
          pre.x = wpre.x + opt_phi_cof * charge_dpre.x;
          pre.y = wpre.y + opt_phi_cof * charge_dpre.y;
        }
        if(pre.x < MIN_PRE) pre.x = MIN_PRE;
        if(pre.y < MIN_PRE) pre.y = MIN_PRE;
        dst[i].x /= pre.x;
        dst[i].y /= pre.y;
    }
    profile.pre += time_since(t0);
}

void update_cells_frame(cell_den_t* cells) {
    TIER* tier = &tier_st[0];
    for( int i = 0; i < tier->cell_cnt; ++i ) {
        cell_den_t* cell = &cells[i];
        pos2_t b0, b1;
        fpos2_t den_pmin = cell->min, den_pmax = cell->max;
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
        cell->binStart = b0;
        cell->binEnd = b1;
    }
}

//MIN_PRE, opt_phi_cof, dampParam

// 2D cGP2D
void update_field_potential(circuit_t* circuit) {
    float** localAr = circuit->cellAreas;
    float** localAr2 = circuit->fillerAreas;
    size_t* bound = circuit->constPinsPerCell;
    cell_den_t* cells = circuit->rects;
    bin_t* bins = circuit->bins;
    area_t* areas = circuit->areas;
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
            cell_den_t* cell = &cells[p];
            pos2_t b0, b1;
            fpos2_t den_pmin = cell->min, den_pmax = cell->max;
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
            cell->binStart = b0;
            cell->binEnd = b1;
            int idx = b0.x * tier->dim_bin.y + b0.y;

            int x = 0, y = 0;
            
            bin_t *bpx = NULL, *bpy = NULL;
            for(x = b0.x, bpx = &bins[idx]; x <= b1.x; x++, bpx += tier->dim_bin.y) {
                idx = x * dim.y + b0.y;
                prec max_x = min(bpx->max.x, den_pmax.x);
                prec min_x = max(bpx->min.x, den_pmin.x);
                for(y = b0.y, bpy = bpx; y <= b1.y; y++, bpy++, idx++) {
                    prec max_y = min(bpy->max.y, den_pmax.y);
                    prec min_y = max(bpy->min.y, den_pmin.y);
                    prec area_share = (max_x - min_x) * (max_y - min_y) * cell->scale;
                    switch(cell->type) {
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
            bins[k].cellArea += localAr[i][k];
        }
        for(int i = 0; i < numThread; ++i) {
            bins[k].fillerArea += localAr2[i][k];
        }
    }
    for(int i = 0; i < tier->tot_bin_cnt; i++) {
        bin_t* bp = &bins[i];
        area_t* a = &areas[i];
        prec area_num2 = bp->cellArea + a->virtArea + a->terminArea;
        prec area_num = area_num2 + bp->fillerArea;
        a->binDensity = area_num * tier->inv_bin_area;
        a->fillerDensity = area_num2 * tier->inv_bin_area;
        __copy_den_to_fft_2D__(a->binDensity, a->coord);
    }
    charge_fft_call(0);
    prec sum_ovf_area = 0;
    for(int i = 0; i < tier->tot_bin_cnt; i++) {
        bin_t* bp = &bins[i];
        area_t* a = &areas[i];    
        __copy_e_from_fft_2D__(&(bp->field), a->coord);
        __copy_phi_from_fft_2D__(&(bp->potential), a->coord);
        gsum_phi += bp->potential * (bp->cellArea + bp->fillerArea + a->terminArea + a->virtArea);
        sum_ovf_area += max((prec)0.0, a->fillerDensity - target_cell_den) * tier->bin_area;
    }
    
    tier->sum_ovf = sum_ovf_area / tier->modu_area;
    gsum_ovf_area += sum_ovf_area;
    gsum_ovfl = gsum_ovf_area / total_modu_area;
}
