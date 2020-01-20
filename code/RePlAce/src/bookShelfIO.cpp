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

#include <cstdio>
#include <cstdlib>
#include <cstring>
//#include        <tgmath.h>
#include <climits>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <mkl.h>

#include "bookShelfIO.h"
#include "bin.h"
#include "global.h"
#include "initPlacement.h"
#include "macro.h"

// this is shared with lefdefIO.cpp
FPOS grow_pmin, grow_pmax;
FPOS terminal_pmin, terminal_pmax;

int numNonRectangularNodes;
int totalShapeCount;


static FPOS shrunk_ratio;

static char error_text[BUFFERSIZE];
static FPOS terminal_size_max;
static FPOS module_size_max;
static POS max_mac_dim;

static dense_hash_map< string, NODES > nodesMap;

void ParseBookShelf() {

    // extract directory
    char dir[511] = {0, };
    extract_dir(bench_aux, dir);

    char fn[511] = {0, };
    char ffn[511] = {0, };
    char sfx[511] = {0, };
    char fn_st[511][511] = {0, };


    // parsing *.aux
    char buf[511] = {0, };
    int buf_size = 511;

    FILE *fp_aux = fopen(bench_aux, "r");

    if( fp_aux == NULL ) {
        cout << "ERROR: Cannot open " << bench_aux << " file." << endl;
        exit(0);
    }
    fgets(buf, buf_size, fp_aux);


    char *token = NULL;
    token = strtok(buf, " \t\n");
    token = strtok(NULL, " \t\n");

    int cnt = 0;
    while((token = strtok(NULL, " \t\n")) != NULL) {
        strcpy(fn_st[cnt], token);
        cnt++;
    }

    // cout << "aux count: " << cnt << endl;
    for(int i = 0; i < cnt; i++) {
        strcpy(fn, fn_st[i]);
        extract_sfx(fn, sfx);
        strcpy(ffn, dir);
        strcat(ffn, fn);

        if(!strcmp(sfx, "nodes")) {
            read_nodes_3D(ffn);
        }
        else if(!strcmp(sfx, "shapes")) {
            read_shapes_3D(ffn);
        }
        else if(!strcmp(sfx, "route")) {
            read_routes_3D(ffn);
        }
        else if(!strcmp(sfx, "nets")) {
            read_nets_3D(ffn);
        }
        else if(!strcmp(sfx, "wts")) {
            // is currently do nothings..
        }
        else if(!strcmp(sfx, "pl")) {
            read_pl2(ffn);
        }
        else if(!strcmp(sfx, "scl")) {
            read_scl(ffn);
        }
    }
    POS tier_min, tier_max;
    int tier_row_cnt = 0;

    get_mms_3d_dim(&tier_min, &tier_max, &tier_row_cnt);
    transform_3d(&tier_min, &tier_max, tier_row_cnt);
    post_read_3d();

}

void transform_3d(POS *tier_min, POS *tier_max, int tier_row_cnt) {
    TERM *curTerminal = NULL;
    TIER *tier = NULL;
    MODULE *curModule = NULL;

    int tot_row_cnt = tier_row_cnt * numLayer;
    // prec site_wid    = 0.0;
    // prec site_spa    = 0.0;
    // prec ori         = 0.0;
    // prec sym         = 0.0;

    for(int i = 0; i < terminalCNT; i++) {
        curTerminal = &terminalInstance[i];

        // lutong
        // if (curTerminal->center.x > grow_pmin.x)
        // curTerminal->center.x
        // = grow_pmin.x + (curTerminal->center.x - grow_pmin.x) *
        // shrunk_ratio.x;
        // if (curTerminal->center.y > grow_pmin.y)
        // curTerminal->center.y
        // = grow_pmin.y + (curTerminal->center.y - grow_pmin.y) *
        // shrunk_ratio.y;

        curTerminal->tier = 0;

        curTerminal->center.z = 0.5 * TIER_DEP;

        curTerminal->pmin.x = curTerminal->center.x - 0.5 * curTerminal->size.x;
        curTerminal->pmin.y = curTerminal->center.y - 0.5 * curTerminal->size.y;
        curTerminal->pmin.z = curTerminal->center.z - 0.5 * curTerminal->size.z;
        curTerminal->pmax.x = curTerminal->center.x + 0.5 * curTerminal->size.x;
        curTerminal->pmax.y = curTerminal->center.y + 0.5 * curTerminal->size.y;
        curTerminal->pmax.z = curTerminal->center.z + 0.5 * curTerminal->size.z;
    }
   
   /* 
    for(int i=0; i<terminalCNT; i++) {
        if(terminalInstance[i].isTerminalNI) {
            continue;
        }   

        if(shapeMap.find( terminalInstance[i].name ) != shapeMap.end() ) {
            continue;
        }

        TERM* curTerminalA = &terminalInstance[i];

        for(int j=i+1; j<terminalCNT; j++) {
            if( terminalInstance[j].isTerminalNI) {
                continue;
            }

            if(shapeMap.find( terminalInstance[j].name ) != shapeMap.end() ) {
                continue;
            }

            TERM* curTerminalB = &terminalInstance[j];

            prec commonArea = pGetCommonAreaXY( curTerminalA->pmin, curTerminalA->pmax, 
                              curTerminalB->pmin, curTerminalB->pmax );
            if(commonArea > 0) {
                cout << "Warning!! " << curTerminalA->name << " , " << curTerminalB->name << " : " << commonArea << endl;
            }
        }
    }
    cout << "end" << endl;
    */

/*
    prec core_min_x = row_st[0].pmin.x;
    prec core_min_y = row_st[0].pmin.y;
    prec core_max_x = row_st[0].pmax.x;
    prec core_max_y = row_st[0].pmax.y;

    for(int i = 0; i < row_cnt; i++) {
        core_min_x = min(core_min_x, row_st[i].pmin.x);
        core_min_y = min(core_min_y, row_st[i].pmin.y);
        core_max_x = max(core_max_x, row_st[i].pmax.x);
        core_max_y = max(core_max_y, row_st[i].pmax.y);
    }
*/
    // 
    // mgwoo
    //
    assert( tot_row_cnt == row_cnt );

//    row_st = (ROW *)realloc(row_st, sizeof(struct ROW) * tot_row_cnt);
    tier_st = (TIER *)mkl_malloc(sizeof(struct TIER) * numLayer, 64);

    placementMacroCNT = 0;

    for(int i = 0; i < moduleCNT; i++) {
        curModule = &moduleInstance[i];

        if(curModule->size.y - rowHeight > Epsilon) {
            placementMacroCNT++;
            curModule->flg = Macro;
        }
        else if(curModule->size.y - rowHeight < -Epsilon) {
            printf("Cell %s   %d height ERROR: %f\n", curModule->name, i,
                   curModule->size.y);
        }
        else {
            curModule->flg = StdCell;
        }
    }

    for(int i = 0; i < numLayer; i++) {
        tier = &tier_st[i];
        tier->row_st = &(row_st[tier_row_cnt * i]);
        tier->row_cnt = tier_row_cnt;
        tier->term_cnt = 0;

        if(terminalCNT > 0) {
            tier->term_st =
                (TERM **)malloc(sizeof(TERM *) * terminalCNT);
        }
        else {
            tier->term_st = NULL;
        }

        tier->modu_cnt = 0;
        tier->modu_st = NULL;
        tier->mac_cnt = 0;
        tier->mac_st = NULL;
        tier->term_area = 0;
        tier->virt_area = 0;
    }

//    int curr_row_hei = tier_min->y;
    
//    FPOS pmin;
//    FPOS pmax;
//    FPOS size;

    // site_wid= row_st[0].site_wid;
    // site_spa= row_st[0].site_spa;
    // ori     = row_st[0].ori;
    // sym     = row_st[0].sym;
//    pmin = row_st[0].pmin;
//    pmax = row_st[0].pmax;
//    size = row_st[0].size;

//    for(int i = 0; i < tier_row_cnt; i++) {
        // for (int j=0; j<numLayer; j++) {
        // row = &row_st[i + tier_row_cnt*j];
        // row->site_wid= site_wid; // lutong
        // row->site_spa= site_spa;
        // row->ori    = ori;
        // row->sym    = sym;
        // row->size.x = tier_max->x - tier_min->x;
        // row->size.y = rowHeight;
        // row->size.z = TIER_DEP;
        // row->pmin.x = tier_min->x;
        // row->pmin.y = curr_row_hei;
        // row->pmin.z = (prec)j * TIER_DEP;
        // row->pmax.x = row->pmin.x + row->size.x;
        // row->pmax.y = row->pmin.y + row->size.y;
        // row->pmax.z = row->pmin.z + row->size.z;
        // row->x_cnt  = tier_max->x - tier_min->x;
        //}
//        curr_row_hei += rowHeight;
//    }


    // update global_variable
    row_cnt = tot_row_cnt;

    for(int i = 0; i < numLayer; i++) {
        tier = &tier_st[i];

        tier->pmin.x = (prec)tier_min->x;
        tier->pmin.y = (prec)tier_min->y;
        tier->pmin.z = (prec)i * TIER_DEP;

        tier->pmax.x = (prec)tier_max->x;
        tier->pmax.y = (prec)tier_max->y;
        tier->pmax.z = (prec)(i + 1) * TIER_DEP;

        tier->size.x = tier->pmax.x - tier->pmin.x;
        tier->size.y = tier->pmax.y - tier->pmin.y;
        tier->size.z = tier->pmax.z - tier->pmin.z;

        tier->center.x = (tier->pmin.x + tier->pmax.x) * 0.5;
        tier->center.y = (tier->pmin.y + tier->pmax.y) * 0.5;
        tier->center.z = ((prec)i + 0.5) * TIER_DEP;

        tier->area = tier->size.x * tier->size.y * tier->size.z;
    }

    // TIER's term_area update!
    // mgwoo
    for(int i = 0; i < terminalCNT; i++) {
        curTerminal = &terminalInstance[i];
        tier = &tier_st[curTerminal->tier];

        tier->term_st[tier->term_cnt++] = curTerminal;
       
        // skip for IO pin informations
//        if( curTerminal->isTerminalNI ) {
//            continue;
//        }

        // if there is no shape's information
        if( shapeMap.find( curTerminal->name ) == shapeMap.end() ) {
            tier->term_area += get_common_area(curTerminal->pmin, curTerminal->pmax,
                    tier->pmin, tier->pmax);
        }
        // if there exist shape's information
        else {
            for(auto& curIdx : shapeMap[curTerminal->name]) {
                prec llx = shapeStor[curIdx].llx, 
                     lly = shapeStor[curIdx].lly,
                     width = shapeStor[curIdx].width,
                     height = shapeStor[curIdx].height;
                FPOS tmpMin(llx, lly, 0), tmpMax(llx + width, lly + height, 1);
                
                tier->term_area += get_common_area(tmpMin, tmpMax,
                        tier->pmin, tier->pmax);
//                prec commonArea = get_common_area(tmpMin, tmpMax,
//                        tier->pmin, tier->pmax);
            }
        }
    }

    for(int i = 0; i < numLayer; i++) {
        tier = &tier_st[i];
        if(tier->term_cnt > 0) {
            tier->term_st = (TERM **)realloc(
                tier->term_st, sizeof(struct TERM *) * tier->term_cnt);
        }
        else if(tier->term_st) {
            free(tier->term_st);
        }
    }
}

//
// update: total_macro_area, 
//          total_std_area, 
//          total_PL_area, 
//          gmov_macro_cnt,
//
//          placementStdcellCNT, placementMacroCNT
//          
//          gmin, gmax, gwid // global minimum, global maximum, global width (plotting)
//
//  malloc : place_st, place_backup_st // place_st_cnt
//
void post_read_3d(void) {
    // int tmpMultiHeightCellCNT = 0;
    struct MODULE *curModule = NULL;
    struct TERM *curTerminal = NULL;

    placementStdcellCNT = 0;
    total_std_area = 0;
    placementMacroCNT = 0;
    total_macro_area = 0;

    if(INPUT_FLG == ISPD) {
        for(int i = 0; i < moduleCNT; i++) {
            curModule = &moduleInstance[i];
            if(get_abs(curModule->size.x - 504) < Epsilon &&
               get_abs(curModule->size.y - 216) < Epsilon) {
                // stupid 504 x 216 cell overlap problem (w/ terminal) ISPD 2005
                // b3 case.
                curModule->flg = StdCell;
                placementStdcellCNT++;
                // tmpMultiHeightCellCNT++;
                // CAUTION!!
                total_macro_area += curModule->area;
            }
            else if(curModule->size.y - 2 * rowHeight > Epsilon) {
                // if (curModule->size.y - rowHeight > 24.0 + Epsilon) {
                // //lutong
                curModule->flg = Macro;
                placementMacroCNT++;
                total_macro_area += curModule->area;
            }
            else if(curModule->size.x - 10 * rowHeight > Epsilon &&
                    curModule->size.y - rowHeight > Epsilon) {
                curModule->flg = Macro;
                placementMacroCNT++;
                total_macro_area += curModule->area;
            }
            else if(curModule->size.y - rowHeight < -Epsilon) {
                printf("*** ERROR:  Cell \"%d\" Height ERROR \n", i);
            }
            else {
                curModule->flg = StdCell;
                placementStdcellCNT++;
                // if (get_abs(curModule->size.y - 2*rowHeight) < Epsilon)
                // tmpMultiHeightCellCNT++;
                // stupid precision reproduction bug
                // for mental health, you should delete the next 5 lines except
                // for 4th line.
                if(curModule->size.y - rowHeight > Epsilon) {
                    total_macro_area += curModule->area;
                }
                else {
                    total_std_area += curModule->area;
                }
                continue;
            }
        }
    }
    else {
        for(int i = 0; i < moduleCNT; i++) {
            curModule = &moduleInstance[i];
            if(curModule->size.y - rowHeight > Epsilon) {
                // if (curModule->size.y - rowHeight > 24.0 + Epsilon) {
                // //lutong
                curModule->flg = Macro;
                placementMacroCNT++;
                total_macro_area += curModule->area;
            }
            else if(curModule->size.y - rowHeight < -Epsilon) {
                printf("*** ERROR:  Cell \"%d\" Height ERROR \n", i);
            }
            else {
                curModule->flg = StdCell;
                placementStdcellCNT++;
                total_std_area += curModule->area;
                continue;
            }
        }
    }

    printf("INFO:  #STDCELL=%d,  #movableMACRO=%d,  #fixedMACRO=%d\n",
           placementStdcellCNT /* - tmpMultiHeightCellCNT*/,
           placementMacroCNT /* + tmpMultiHeightCellCNT*/, terminalCNT);

    gmov_mac_cnt = placementMacroCNT;

    
    // initial malloc : row_cnt
    place_st_cnt = 0;
    place_st = (PLACE *)malloc(sizeof(struct PLACE) * row_cnt);

    ROW *last_row = NULL;
    PLACE *curPlace = NULL;

    // 
    // Assume, height is equally distributed.
    //
    int place_hei = 0;
    for(int i = 0; i < row_cnt; i++) {
        ROW* row = &row_st[i];
        if(i == 0) {
            place_hei = row->size.y;
        }
        else if(row->pmin.x != last_row->pmin.x ||
                row->pmax.x != last_row->pmax.x ||

                row->pmin.y != last_row->pmax.y ||

                row->pmin.z != last_row->pmin.z ||
                row->pmax.z != last_row->pmax.z) {
            curPlace = &place_st[place_st_cnt++];

            curPlace->org.x = last_row->pmin.x;
            curPlace->org.y = last_row->pmax.y - place_hei;
            curPlace->org.z = last_row->pmin.z;

            curPlace->end.x = last_row->pmax.x;
            curPlace->end.y = last_row->pmax.y;
            curPlace->end.z = last_row->pmax.z;

            curPlace->stp.x = SITE_SPA;
            curPlace->stp.y = rowHeight;
            curPlace->stp.z = TIER_DEP;

            curPlace->cnt.x = curPlace->end.x - curPlace->org.x;
            curPlace->cnt.y = curPlace->end.y - curPlace->org.y;
            curPlace->cnt.z = curPlace->end.z - curPlace->org.z;

            curPlace->area = curPlace->cnt.x * curPlace->cnt.y * curPlace->cnt.z;

            curPlace->center.x = 0.5 * (curPlace->org.x + curPlace->end.x);
            curPlace->center.y = 0.5 * (curPlace->org.y + curPlace->end.y);
            curPlace->center.z = 0.5 * (curPlace->org.z + curPlace->end.z);

            place_hei = row->size.y;
        }
        else {
            place_hei += row->size.y;
        }
        last_row = row;
    }

    curPlace = &place_st[place_st_cnt++];

    curPlace->org.x = last_row->pmin.x;
    curPlace->org.y = last_row->pmax.y - place_hei;
    curPlace->org.z = last_row->pmin.z;

    curPlace->end.x = last_row->pmax.x;
    curPlace->end.y = last_row->pmax.y;
    curPlace->end.z = last_row->pmax.z;

    curPlace->stp.x = SITE_SPA;
    curPlace->stp.y = rowHeight;
    curPlace->stp.z = TIER_DEP;

    curPlace->cnt.x = curPlace->end.x - curPlace->org.x;
    curPlace->cnt.y = curPlace->end.y - curPlace->org.y;
    curPlace->cnt.z = curPlace->end.z - curPlace->org.z;

    curPlace->area = curPlace->cnt.x * curPlace->cnt.y * curPlace->cnt.z;

    curPlace->center.x = 0.5 * (curPlace->org.x + curPlace->end.x);
    curPlace->center.y = 0.5 * (curPlace->org.y + curPlace->end.y);
    curPlace->center.z = 0.5 * (curPlace->org.z + curPlace->end.z);

//    cout << "current place_st_cnt Final: " << place_st_cnt << endl;
//    cout << "previous row_cnt: " << row_cnt << endl;

//    for(int i=0; i<place_st_cnt; i++) {
//        place_st[i].dump(to_string(i));
//    }

    // second re-malloc : upto place_st_cnt;
    place_st = (PLACE *)realloc(place_st, sizeof(struct PLACE) * place_st_cnt);
    place_backup_st =
        (PLACE *)realloc(place_backup_st, sizeof(struct PLACE) * place_st_cnt);

    for(int i = 0; i < place_st_cnt; i++) {
        place_backup_st[i] = place_st[i];
    }

    printf("INFO:  #ROW=%d,  ROW HEIGHT=%d, TIER DEPTH=%d\n", row_cnt,
           (int)(rowHeight + 0.5), (int)(TIER_DEP + 0.5));

    // global variable 'place' update
    place.stp.x = SITE_SPA;
    place.stp.y = 1.0; // rowHeight....................????
    place.stp.z = TIER_DEP;

    place.org = place_st[0].org;
    place.end = place_st[0].end;

    total_PL_area = 0;
    for(int i = 0; i < place_st_cnt; i++) {
        curPlace = &place_st[i];

        place.org.x = min(place.org.x, curPlace->org.x);
        place.org.y = min(place.org.y, curPlace->org.y);
        place.org.z = min(place.org.z, curPlace->org.z);

        place.end.x = max(place.end.x, curPlace->end.x);
        place.end.y = max(place.end.y, curPlace->end.y);
        place.end.z = max(place.end.z, curPlace->end.z);

        total_PL_area += curPlace->area;
    }

    // global min & global max variable setting
    gmin = place.org;
    gmax = place.end;

    for(int i = 0; i < terminalCNT; i++) {
        curTerminal = &terminalInstance[i];

        gmin.x = min(gmin.x, curTerminal->pmin.x);
        gmin.y = min(gmin.y, curTerminal->pmin.y);
        gmin.z = min(gmin.z, curTerminal->pmin.z);

        gmax.x = max(gmax.x, curTerminal->pmax.x);
        gmax.y = max(gmax.y, curTerminal->pmax.y);
        gmax.z = max(gmax.z, curTerminal->pmax.z);
    }
    
    gwid.x = gmax.x - gmin.x;
    gwid.y = gmax.y - gmin.y;
    gwid.z = gmax.z - gmin.z;

/*
// setting additional margin 
// for plotting.. here It is not required

    gmin.Dump("gmin first");

    prec lx_mg = place.org.x - gmin.x, 
         ly_mg = place.org.y - gmin.y,
         rx_mg = gmax.x - place.end.x,
         ry_mg = gmax.y - place.end.y;
    
//    cout << "lx_mg: " << lx_mg << ", ly_mg: " << ly_mg << endl;
//    cout << "rx_mg: " << rx_mg << ", ry_mg: " << ry_mg << endl;

    prec max_mg = PREC_MIN;
    max_mg = max(max_mg, lx_mg);
    max_mg = max(max_mg, ly_mg);
    max_mg = max(max_mg, rx_mg);
    max_mg = max(max_mg, ry_mg);
    
//    cout << "max_mg: " << max_mg << endl;

    gmin.x = place.org.x - max_mg;
    gmin.y = place.org.y - max_mg;
    gmin.z = place.org.z;

    gmax.x = place.end.x + max_mg;
    gmax.y = place.end.y + max_mg;
    gmax.z = place.end.z;
*/

    printf("INFO:  GLOBAL.SPACE.MIN = (%f, %f, %f)\n", gmin.x, gmin.y, gmin.z);
    printf("INFO:  GLOBAL.SPACE.MAX = (%f, %f, %f)\n", gmax.x, gmax.y, gmax.z);

//    gwid.Dump("gwid");

    printf("INFO:  PLACE.ORIGIN = (%d, %d, %d)\n", (int)(place.org.x + 0.5),
           (int)(place.org.y + 0.5), (int)(place.org.z + 0.5));
    printf("INFO:  PLACE.END = (%d, %d, %d)\n\n", (int)(place.end.x + 0.5),
           (int)(place.end.y + 0.5), (int)(place.end.z + 0.5));

    place.cnt.x = place.end.x - place.org.x;
    place.cnt.y = place.end.y - place.org.y;
    place.cnt.z = place.end.z - place.org.z;

    place.center.x = 0.5 * (place.org.x + place.end.x);
    place.center.y = 0.5 * (place.org.y + place.end.y);
    place.center.z = 0.5 * (place.org.z + place.end.z);

    place.area = place.cnt.x * place.cnt.y * place.cnt.z;
//    place.dump("global variable place");
    
}
/*
void post_read_2d(void) {
    struct PLACE *place0 = NULL;
    struct MODULE *curModule = NULL;
    struct TERM *curTerminal = NULL;
    struct ROW *row = NULL, *last_row = NULL;

    total_std_area = 0;
    placementStdcellCNT = 0;
    placementMacroCNT = 0;
    total_macro_area = 0;

    for(int i = 0; i < moduleCNT; i++) {
        curModule = &moduleInstance[i];
        if(curModule->size.y - (prec)rowHeight > Epsilon) {
            curModule->flg = Macro;
            placementMacroCNT++;
            total_macro_area += curModule->area;
        }
        else if(curModule->size.y - (prec)rowHeight < -Epsilon) {
            printf("cell %d height ERROR \n", i);
        }
        else {
            curModule->flg = StdCell;
            placementStdcellCNT++;
            total_std_area += curModule->area;
            continue;
        }
    }

    printf(
        "\n\n# std_cell = %d , #mov_mac = %d , # fixed_mac = %d\n\n\n",
        placementStdcellCNT, placementMacroCNT, terminalCNT);

    gmov_mac_cnt = placementMacroCNT;

    place_st_cnt = 0;

    place_st = (struct PLACE *)malloc(sizeof(struct PLACE) * row_cnt);

    int place_hei = 0;
    for(int i = 0; i < row_cnt; i++) {
        row = &row_st[i];
    
        if(i == 0) {
            place_hei = row->size.y;
        }
        else if(row->pmin.x != last_row->pmin.x ||
                row->pmax.x != last_row->pmax.x ||
                row->pmin.y != last_row->pmax.y ||
                row->pmin.z != last_row->pmin.z ||
                row->pmax.z != last_row->pmax.z) {
            place0 = &place_st[place_st_cnt++];

            place0->org.x = last_row->pmin.x;
            place0->org.y = last_row->pmax.y - place_hei;
            place0->org.z = last_row->pmin.z;

            place0->end.x = last_row->pmax.x;
            place0->end.y = last_row->pmax.y;
            place0->end.z = last_row->pmax.z;

            place0->stp.x = SITE_SPA;
            place0->stp.y = rowHeight;
            place0->stp.z = TIER_DEP;

            place0->cnt.x = place0->end.x - place0->org.x;
            place0->cnt.y = place0->end.y - place0->org.y;
            place0->cnt.z = place0->end.z - place0->org.z;

            place0->area = place0->cnt.x * place0->cnt.y * place0->cnt.z;

            place0->center.x = 0.5 * (place0->org.x + place0->end.x);
            place0->center.y = 0.5 * (place0->org.y + place0->end.y);
            place0->center.z = 0.5 * (place0->org.z + place0->end.z);

            place_hei = row->size.y;
        }
        else {
            place_hei += row->size.y;
        }
        last_row = row;
    }

    place0 = &place_st[place_st_cnt++];

    place0->org.x = last_row->pmin.x;
    place0->org.y = last_row->pmax.y - place_hei;
    place0->org.z = last_row->pmin.z;

    place0->end.x = last_row->pmax.x;
    place0->end.y = last_row->pmax.y;
    place0->end.z = last_row->pmax.z;

    place0->stp.x = SITE_SPA;
    place0->stp.y = rowHeight;
    place0->stp.z = TIER_DEP;

    place0->cnt.x = place0->end.x - place0->org.x;
    place0->cnt.y = place0->end.y - place0->org.y;
    place0->cnt.z = place0->end.z - place0->org.z;

    place0->area = place0->cnt.x * place0->cnt.y * place0->cnt.z;

    place0->center.x = 0.5 * (place0->org.x + place0->end.x);
    place0->center.y = 0.5 * (place0->org.y + place0->end.y);
    place0->center.z = 0.5 * (place0->org.z + place0->end.z);

    place_st = (PLACE *)realloc(place_st, sizeof(struct PLACE) * place_st_cnt);

    total_PL_area = 0;

    for(int i = 0; i < place_st_cnt; i++) {
        place0 = &place_st[i];
        total_PL_area += place0->area;
    }

    place_hei = 0;

    place.stp.x = SITE_SPA;
    place.stp.y = 1.0 // rowHeight //;
    place.stp.z = TIER_DEP;

    place.org = place_st[0].org;
    place.end = place_st[0].end;

    printf("# rows = %d , row_height = %d, tier_depth = %d\n", row_cnt,
           (int)(rowHeight + 0.5), (int)(TIER_DEP + 0.5));

    total_PL_area = 0;

    for(int i = 0; i < place_st_cnt; i++) {
        place0 = &place_st[i];


        //   printf("PL rect %d --> min = (%d , %d, %d) , \
        //   max = (%d , %d, %d) , dim = (%d , %d, %d)\n",
        //     i+1 ,
        //     (int) (place0->org.x + 0.5) ,
        //     (int) (place0->org.y + 0.5) ,
        //     (int) (place0->org.z + 0.5) ,
        //     (int) (place0->end.x + 0.5) ,
        //     (int) (place0->end.y + 0.5) ,
        //     (int) (place0->end.z + 0.5) ,
        //     (int) (place0->cnt.x + 0.5) ,
        //     (int) (place0->cnt.y + 0.5) ,
        //     (int) (place0->cnt.z + 0.5) );


        place.org.x = min(place.org.x, place0->org.x);
        place.org.y = min(place.org.y, place0->org.y);
        place.org.z = min(place.org.z, place0->org.z);

        place.end.x = max(place.end.x, place0->end.x);
        place.end.y = max(place.end.y, place0->end.y);
        place.end.z = max(place.end.z, place0->end.z);

        total_PL_area += place0->area;
    }

    ///////////////////////
    // exit(0);
    // place.end.z = 2;
    ///////////////////////

    //////////////////////////////////////////
    if(INPUT_FLG == IBM)
        place.end.y += rowHeight;
    //////////////////////////////////////////

    gmin = place.org;
    gmax = place.end;

    // gmin = terminalInstance[0].pmin ; //
    // gmax = terminalInstance[0].pmax ; //

    for(int i = 0; i < terminalCNT; i++) {
        curTerminal = &terminalInstance[i];

        gmin.x = min(gmin.x, curTerminal->pmin.x);
        gmin.y = min(gmin.y, curTerminal->pmin.y);
        gmin.z = min(gmin.z, curTerminal->pmin.z);

        gmax.x = max(gmax.x, curTerminal->pmax.x);
        gmax.y = max(gmax.y, curTerminal->pmax.y);
        gmax.z = max(gmax.z, curTerminal->pmax.z);
    }

    prec lx_mg = 0, ly_mg = 0;
    prec rx_mg = 0, ry_mg = 0;
    prec max_mg = 0;

    lx_mg = place.org.x - gmin.x;
    ly_mg = place.org.y - gmin.y;
    rx_mg = gmax.x - place.end.x;
    ry_mg = gmax.y - place.end.y;

    max_mg = max(max_mg, lx_mg);
    max_mg = max(max_mg, ly_mg);
    max_mg = max(max_mg, rx_mg);
    max_mg = max(max_mg, ry_mg);

    gmin.x = place.org.x - max_mg;
    gmin.y = place.org.y - max_mg;
    gmin.z = place.org.z;

    gmax.x = place.end.x + max_mg;
    gmax.y = place.end.y + max_mg;
    gmax.z = place.end.z;

    printf(" ==== %f %f %f %f %f %f\n", gmin.x, gmin.y, gmin.z, place.org.x,
           place.org.y, place.org.z);

    printf(" ==== %f %f %f %f %f %f\n", gmax.x, gmax.y, gmax.z, place.end.x,
           place.end.y, place.end.z);

    gwid.x = gmax.x - gmin.x;
    gwid.y = gmax.y - gmin.y;
    gwid.z = gmax.z - gmin.z;

    printf("place.org = (%d, %d, %d) , place.end = (%d, %d, %d)\n",
           (int)(place.org.x + 0.5), (int)(place.org.y + 0.5),
           (int)(place.org.z + 0.5), (int)(place.end.x + 0.5),
           (int)(place.end.y + 0.5), (int)(place.end.z + 0.5));

    place.cnt.x = place.end.x - place.org.x;
    place.cnt.y = place.end.y - place.org.y;
    place.cnt.z = place.end.z - place.org.z;

    place.center.x = 0.5 * (place.org.x + place.end.x);
    place.center.y = 0.5 * (place.org.y + place.end.y);
    place.center.z = 0.5 * (place.org.z + place.end.z);

    place.area = place.cnt.x * place.cnt.y * place.cnt.z;
}

int read_nodes (char *input) {
    FILE*fp=fopen(input,"r");
    char*token=NULL;
    char buf[BUF_SZ]// ,buf2[BUF_SZ] ;
    int buf_size=BUF_SZ;
    int i=0;

    MODULE *curModule=NULL;
    TERM *curTerminal=NULL;
    FPOS max_cell,min_cell,max_term,min_term;

    int idx=0;

    total_term_area = 0 ;

    max_cell.x = 0,max_cell.y = 0;
    min_cell.x = 100000000,min_cell.y=100000000;

    max_term.x = 0,max_term.y = 0;
    min_term.x = 100000000,min_term.y=100000000;

    avgCellSize.x = 0,avgCellSize.y = 0;
    avgTerminalSize.x = 0,avgTerminalSize.y = 0;


    do
    {
        fgets(buf,buf_size,fp);
        token=strtok(buf," \t\n") ;
    }
    while(!token || token[0]=='#' || !strcmp(token,"UCLA") );
    token = strtok(NULL," \t\n");
    token = strtok(NULL," \t\n");

    moduleCNT = atoi(token);

    do
    {
        fgets(buf,buf_size,fp);
        token=strtok(buf," \t\n") ;
    }
    while(!token || token[0]=='#' || !strcmp(token,"UCLA") );
    token = strtok(NULL," \t\n");
    token = strtok(NULL," \t\n");

    terminalCNT = atoi(token);

    moduleCNT -= terminalCNT;

    moduleInstance=(struct MODULE*)mkl_malloc(sizeof(struct
MODULE)*moduleCNT,64);
    terminalInstance=(struct TERM*)mkl_malloc(sizeof(struct
TERM)*terminalCNT,64);

    modu_map_cnt = 100000;
    term_map_cnt = 100000;

    modu_map = (int*)mkl_malloc(sizeof(int)*modu_map_cnt, 64);
    term_map = (int*)mkl_malloc(sizeof(int)*term_map_cnt, 64);

    for(i=0;i<modu_map_cnt;i++)
    {
        modu_map[i] = -1 ;
    }
    for(i=0;i<term_map_cnt;i++)
    {
        term_map[i] = -1 ;
    }


    for (i=0; i<moduleCNT; i++) {
        curModule = &moduleInstance[i];
        curModule->idx = i;

        do {
            fgets (buf, buf_size, fp);
            token = strtok(buf, " \t\n");
        } while (!token || token[0]=='#' || !strcmp(token, "UCLA"));

        strcpy (curModule->name, token);
        idx = atoi(token+1);

        if (idx > modu_map_cnt-1) {
            modu_map_cnt = max ( idx+1,modu_map_cnt*2 );
            modu_map = (int*)realloc(modu_map,sizeof(int)*modu_map_cnt);
        }

        modu_map [idx] = i ;

        token=strtok(NULL," \t\n");
        curModule->size.x = atoi(token) ;

        token=strtok(NULL," \t\n");
        curModule->size.y = atoi(token) ;

        curModule->area = curModule->size.x * curModule->size.y ;

        curModule->pof = NULL;
        curModule->pin = NULL;

        curModule->netCNTinObject=0;

        curModule->pinCNTinObject=0;

        if (max_cell.x < curModule->size.x) max_cell.x = curModule->size.x;
        if (min_cell.x > curModule->size.x) min_cell.x = curModule->size.x;
        if (max_cell.y < curModule->size.y) max_cell.y = curModule->size.y;
        if (min_cell.y > curModule->size.y) min_cell.y = curModule->size.y;
        avgCellSize.x += curModule->size.x;
        avgCellSize.y += curModule->size.y;
    }


    for(i=0;i<terminalCNT;i++)
    {
        curTerminal = & terminalInstance[i];

        curTerminal->idx = i ;
        do
        {
            fgets(buf,buf_size,fp);
            token=strtok(buf," \t\n") ;
        }
        while(!token || token[0]=='#' || !strcmp(token,"UCLA") );


        strcpy(curTerminal->name,token);

        idx = atoi(token+1);

        if(idx > term_map_cnt - 1)
        {
            term_map_cnt = max ( idx+1,term_map_cnt*2 );
            term_map = (int*)realloc(term_map,sizeof(int)*term_map_cnt);
        }

        term_map [idx] = i ;

        token=strtok(NULL," \t\n");
        curTerminal->size.x = atoi(token) ;

        token=strtok(NULL," \t\n");
        curTerminal->size.y = atoi(token) ;

        curTerminal->area = curTerminal->size.x * curTerminal->size.y ;

        curTerminal->pof = NULL ;
        curTerminal->pin = NULL ;

        curTerminal->netCNTinObject=0;
        curTerminal->pinCNTinObject=0;

        if(max_term.x < curTerminal->size.x)
            max_term.x = curTerminal->size.x ;
        if(min_term.x > curTerminal->size.x)
            min_term.x = curTerminal->size.x ;

        if(max_term.y < curTerminal->size.y)
            max_term.y = curTerminal->size.y ;
        if(min_term.y > curTerminal->size.y)
            min_term.y = curTerminal->size.y ;

        avgTerminalSize.x += curTerminal->size.x ;
        avgTerminalSize.y += curTerminal->size.y ;

        total_term_area += curTerminal->area ;
    }


    avgCellSize.x /= (prec)moduleCNT;
    avgCellSize.y /= (prec)moduleCNT;
    avgTerminalSize.x /= (prec)terminalCNT;
    avgTerminalSize.y /= (prec)terminalCNT;

    printf("INFO:  cell min (%.0f,%.0f) , max (%.0f,%.0f) avg (%.2f %.2f)\n",
            min_cell.x,min_cell.y,max_cell.x,max_cell.y,avgCellSize.x,avgCellSize.y);
    printf("curTerminal min (%.0f,%.0f) , max (%.0f,%.0f) avg (%.2f %.2f)\n",
            min_term.x,min_term.y,max_term.x,max_term.y,avgTerminalSize.x,avgTerminalSize.y);
    return 1;
}
*/

int read_shapes_3D(char *input) {

    FILE *fp = fopen(input, "r");
    if(fp == NULL) {
        cout << "** ERROR: Cannot open: " << input << " file (bookshelfIO)" << endl;
        exit(1);
    }

    char temp[BUFFERSIZE];
    char line[LINESIZE];

    bool flag = false;
//    bool isFirstNode = true;

    int curShapeIdx = 0;
    totalShapeCount = 0;

    shapeMap.set_empty_key("!@#!@#");

    char instName[255] = {0,};
    while(!feof(fp)) {
        *line = '\0';
        fgets(line, LINESIZE, fp);
        sscanf(line, "%s%*s", temp);

        if(strlen(line) < 5 || temp[0] == '#' || strcmp(temp, "shapes") == 0)
            continue;

        if(flag) {
            int curShapeCount;
            sscanf(line, "%*s %s", temp);

            if(strcmp(temp, ":") == 0) {
                sscanf(line, "%s", instName);
//                cout << instName << endl;

                sscanf(line, "%*s %*s %d", &curShapeCount);
                totalShapeCount += curShapeCount;
                curShapeIdx = 0;   
            }
            else {
                char name[255] = {0, };
                double llx, lly, width, height;

                sscanf(line, "%s %lf %lf %lf %lf", name, &llx, &lly, &width,
                       &height);
                
                auto it = shapeMap.find(instName);
                if( it == shapeMap.end() ) {
                    vector<int> tmpStor;
                    tmpStor.reserve(curShapeCount);
                    tmpStor.push_back(shapeStor.size());

                    shapeMap[instName] = tmpStor;
                }
                else{ 
                    shapeMap[instName].push_back( shapeStor.size() );
                }
                
                shapeStor.push_back( SHAPE(name, instName, curShapeIdx, 
                                (prec)llx, (prec)lly, (prec)width, (prec)height) );

                curShapeIdx++;
            }
        }

        if(strcmp(temp, "NumNonRectangularNodes") == 0) {
            sscanf(line, "%*s %*s %d", &numNonRectangularNodes);
            cout << "INFO:  #NonRectangularNodes=" << numNonRectangularNodes
                 << endl;
            flag = true;
        }
    }

/*
 * Check whether SHAPE's are correctly parsed 
    cout << "shapeMapSize: " << shapeMap.size() << endl;
    for(auto& curShape : shapeMap) {
        cout << "instName: " << curShape.first << endl;
        for(auto& curIdx : curShape.second) {
            shapeStor[curIdx].dump();
        }
        cout << endl;
    }
    */
    return 1;
}


int read_routes_3D(char *input) {
    char name[255];
    char *token = NULL;
    char temp[BUFFERSIZE];
    char line[LINESIZE];
    int numBlockedLayers = 0;
    bool blockageFlag = false;
    bool beolFlag = false;
    bool edgeFlag = false;
    int e1, e2, e3, e4, e5, e6, e7;
    double temp_gridLLx, temp_gridLLy;
    double temp_tileWidth, temp_tileHeight;
    double temp_blockagePorosity;
    vector< int > blockageLayers;

    FILE *fp;
    if((fp = fopen(input, "r")) == NULL) {
        sprintf(error_text, "bookshelf_IO: Cannot open: %s file", input);
        runtimeError(error_text);
    }

    while(!feof(fp)) {
        *line = '\0';
        fgets(line, LINESIZE, fp);
        sscanf(line, "%s%*s", temp);

        if(strlen(line) < 5 || temp[0] == '#' || strcmp(temp, "route") == 0)
            continue;
        if(strcmp(temp, "NumBlockageNodes") == 0) {
            blockageFlag = true;
            continue;
        }
        if(strcmp(temp, "NumEdgeCapacityAdjustments") == 0) {
            blockageFlag = false;
            edgeFlag = true;
            continue;
        }
        if(strcmp(temp, "Grid") == 0) {
            beolFlag = true;
        }
        if(strcmp(temp, "NumNiTerminals") == 0) {
            beolFlag = false;
            continue;
        }

        if(beolFlag) {
            sscanf(line, "%s :%*s", temp);
            if(strcmp(temp, "Grid") == 0) {
                sscanf(line, "%*s : %d %d %d", &nXgrids, &nYgrids, &nMetLayers);
            }
            else if(strcmp(temp, "VerticalCapacity") == 0) {
                verticalCapacity.clear();
                token = strtok(line, " \t\n");
                token = strtok(NULL, " \t\n");
                token = strtok(NULL, " \t\n");
                for(int i = 0; i < nMetLayers; i++) {
                    verticalCapacity.push_back(atoi(token));
                    token = strtok(NULL, " \t\n");
                }
            }
            else if(strcmp(temp, "HorizontalCapacity") == 0) {
                horizontalCapacity.clear();
                token = strtok(line, " \t\n");
                token = strtok(NULL, " \t\n");
                token = strtok(NULL, " \t\n");
                for(int i = 0; i < nMetLayers; i++) {
                    horizontalCapacity.push_back(atoi(token));
                    token = strtok(NULL, " \t\n");
                }
            }
            else if(strcmp(temp, "MinWireWidth") == 0) {
                minWireWidth.clear();
                token = strtok(line, " \t\n");
                token = strtok(NULL, " \t\n");
                token = strtok(NULL, " \t\n");
                for(int i = 0; i < nMetLayers; i++) {
                    minWireWidth.push_back(atof(token));
                    token = strtok(NULL, " \t\n");
                }
            }
            else if(strcmp(temp, "MinWireSpacing") == 0) {
                minWireSpacing.clear();
                token = strtok(line, " \t\n");
                token = strtok(NULL, " \t\n");
                token = strtok(NULL, " \t\n");
                for(int i = 0; i < nMetLayers; i++) {
                    minWireSpacing.push_back(atof(token));
                    token = strtok(NULL, " \t\n");
                }
            }
            else if(strcmp(temp, "ViaSpacing") == 0) {
                viaSpacing.clear();
                token = strtok(line, " \t\n");
                token = strtok(NULL, " \t\n");
                token = strtok(NULL, " \t\n");
                for(int i = 0; i < nMetLayers; i++) {
                    viaSpacing.push_back(atof(token));
                    token = strtok(NULL, " \t\n");
                }
            }
            else if(strcmp(temp, "GridOrigin") == 0) {
                sscanf(line, "%*s : %lf %lf", &temp_gridLLx, &temp_gridLLy);
                gridLLx = (prec)temp_gridLLx;
                gridLLy = (prec)temp_gridLLy;
            }
            else if(strcmp(temp, "TileSize") == 0) {
                sscanf(line, "%*s : %lf %lf", &temp_tileWidth,
                       &temp_tileHeight);
                tileWidth = (prec)temp_tileWidth;
                tileHeight = (prec)temp_tileHeight;
            }
            else if(strcmp(temp, "BlockagePorosity") == 0) {
                sscanf(line, "%*s : %lf", &temp_blockagePorosity);
                blockagePorosity = (prec)temp_blockagePorosity;
            }
        }

        if(blockageFlag) {
            blockageLayers.clear();
            sscanf(line, "%s %d%*s", name, &numBlockedLayers);
            token = strtok(line, " \t\n");
            token = strtok(NULL, " \t\n");
            token = strtok(NULL, " \t\n");
            for(int i = 0; i < numBlockedLayers; i++) {
                blockageLayers.push_back(atoi(token));
                token = strtok(NULL, " \t\n");
            }
            routeBlockageNodes[name] = blockageLayers;
        }
        if(edgeFlag) {
            sscanf(line, "%d %d %d %d %d %d %d ", &e1, &e2, &e3, &e4, &e5, &e6,
                   &e7);
            edgeCapAdj.push_back(std::make_tuple(e1, e2, e3, e4, e5, e6, e7));
        }
    }
    return 1;
}

// Input : .nodes
// output : moduleInstance, moduleCNT, module_size_max, 
// terminalInstance, terminalCNT, terminal_size_max,
// max_mac_dim ( same as max_cell )
//
// note that module == cell, area == size
// mostly focusing moduleInstance, teminalInstance
//

int read_nodes_3D(char *input) {
    nodesMap.set_empty_key("!@#!@#");

    FILE *fp = fopen(input, "r");
    char *token = NULL;
    char buf[255];

    MODULE *maxSizeCell = NULL;
    MODULE *minSizeCell = NULL;

    MODULE *curModule = NULL;
    TERM *curTerminal = NULL;

    
    int buf_size = 255;
    int idx = 0;
    int totalInstanceCNT = 0;
    // int isTerminal = 0;
    // int isTerminalNI = 0;
    int first_term = 1;
    int first_modu = 1;


    do {
        fgets(buf, buf_size, fp);
        token = strtok(buf, " \t\n");
    } while(!token || token[0] == '#' || !strcmp(token, "UCLA"));
    token = strtok(NULL, " \t\n");
    token = strtok(NULL, " \t\n");
    totalInstanceCNT = moduleCNT = atoi(token);

    do {
        fgets(buf, buf_size, fp);
        token = strtok(buf, " \t\n");
    } while(!token || token[0] == '#' || !strcmp(token, "UCLA"));
    token = strtok(NULL, " \t\n");
    token = strtok(NULL, " \t\n");
    terminalCNT = atoi(token);

    moduleCNT -= terminalCNT;
    moduleInstance =
        (struct MODULE *)mkl_malloc(sizeof(struct MODULE) * moduleCNT, 64);
    terminalInstance =
        (struct TERM *)mkl_malloc(sizeof(struct TERM) * terminalCNT, 64);

    // why this is required..
    // terminal_size_max = zeroFPoint;
    // module_size_max = zeroFPoint;

    // to find max size..
    prec maxSize = PREC_MIN;
    prec minSize = PREC_MAX;

    // for (int i=0; i<totalInstanceCNT; i++) {
    //

    // cursor for (module, terminal) Instance
    int moduleCur = 0, terminalCur = 0;
    int cnt = 0;
    
    // max & min initialize
    FPOS maxCell( PREC_MIN, PREC_MIN, TIER_DEP );
    FPOS minCell( PREC_MAX, PREC_MAX, TIER_DEP );
    FPOS maxTerm( PREC_MIN, PREC_MIN, TIER_DEP );
    FPOS minTerm( PREC_MAX, PREC_MAX, TIER_DEP );
    
    char nodeName[255], line[LINESIZE];
    
    FPOS avgCellSize, avgTerminalSize;
    
    while(!feof(fp)) {
        *line = '\0';
        fgets(line, LINESIZE, fp);
        sscanf(line, "%s\t%*s\n", nodeName);
        // cout << "current Line: " << line << endl;
        // cout << "current nodeName: " << nodeName << endl;

        if(strlen(line) < 5 || nodeName[0] == '#' ||
           strcmp(nodeName, "UCLA") == 0)
            continue;

        if(strcmp(nodeName, "NumNodes") == 0 ||
           strcmp(nodeName, "NumTerminals") == 0)
            continue;

        cnt++;
        char node_type[BUFFERSIZE];
        *node_type = '\0';

        prec x, y;
        prec z = TIER_DEP;

#if PREC_MODE == IS_FLOAT
        if(flg_3dic_io) {
            sscanf(line, "%s%f%f%f%s\n", nodeName, &x, &y, &z, node_type);
        }
        else {
            sscanf(line, "%s%f%f%s\n", nodeName, &x, &y, node_type);
        }
#elif PREC_MODE == IS_DOUBLE
        if(flg_3dic_io) {
            sscanf(line, "%s%lf%lf%lf%s\n", nodeName, &x, &y, &z, node_type);
        }
        else {
            sscanf(line, "%s%lf%lf%s\n", nodeName, &x, &y, node_type);
        }
#endif

        // cout << "current parsed: " << nodeName << ", " << x<< ", " <<
        // y << ", " << node_type << endl;
        // do {
        // fgets (buf, buf_size, fp);
        // token = strtok(buf, " \t\n");
        // }
        // while(!token || token[0]=='#' || !strcmp(token,"UCLA"));

        bool isTerminal = false, isTerminalNI = false;

        // get Index..
        if(strcmp(node_type, "") == 0) {
            // cout << "None!" << endl;
            idx = moduleCur++;
        }
        else if(strcmp(node_type, "terminal") == 0) {
            // cout << "terminal!!!" << endl;
            idx = terminalCur++;
            isTerminal = true;
        }
        else if(strcmp(node_type, "terminal_NI") == 0) {
            // cout << "terminalNI!!!" << endl;
            idx = terminalCur++;
            isTerminal = true;
            isTerminalNI = true;
        }

        // update NodesMap
        nodesMap[string(nodeName)] = NODES(idx, isTerminal, isTerminalNI);

        // moduleInstance update
        //
        // if normal basic nodes...( not terminal & not terminalNI )
        if(!isTerminal && !isTerminalNI) {
            curModule = &moduleInstance[idx];
            curModule->idx = idx;
            strcpy(curModule->name, nodeName);

            // width info
            // update size, half_size
            curModule->size.x = x;
            curModule->half_size.x = 0.5 * curModule->size.x;

            // height info
            // update size, half_size
            curModule->size.y = y;
            curModule->half_size.y = 0.5 * curModule->size.y;

            curModule->size.z = z;
            curModule->half_size.z = 0.5 * curModule->size.z;

            // update area
            curModule->area =
                curModule->size.x * curModule->size.y * curModule->size.z;

            // why this in here....
            curModule->pof = NULL;
            curModule->pin = NULL;
            curModule->netCNTinObject = 0;
            curModule->pinCNTinObject = 0;

            // max_min update for x.
            if(maxCell.x < curModule->size.x)
                maxCell.x = curModule->size.x;
            if(minCell.x > curModule->size.x)
                minCell.x = curModule->size.x;

            // max_min update for y
            if(maxCell.y < curModule->size.y)
                maxCell.y = curModule->size.y;
            if(minCell.y > curModule->size.y)
                minCell.y = curModule->size.y;

            // max update for area(size)
            if(maxSize < curModule->area) {
                maxSize = curModule->area;
                maxSizeCell = curModule;
            }

            // min update for area(size)
            if(minSize > curModule->area) {
                minSize = curModule->area;
                minSizeCell = curModule;
            }

            // update module_size_max (x,y)
            if(first_modu) {
                first_modu = 0;
                module_size_max = curModule->size;
            }
            else {
                module_size_max.x = max(module_size_max.x, curModule->size.x);
                module_size_max.y = max(module_size_max.y, curModule->size.y);
            }

            // for average calculation
            avgCellSize.x += curModule->size.x;
            avgCellSize.y += curModule->size.y;
        }
        // terminalInstance update
        // both of terminal & terminal_NI are saved in here!
        else {
            curTerminal = &terminalInstance[idx];
            curTerminal->idx = idx;
            strcpy(curTerminal->name, nodeName);

            curTerminal->size.x = x;
            curTerminal->size.y = y;

            // if Non-Image mode, ignore width & height
            if(isTerminalNI) {
                curTerminal->size.x = curTerminal->size.y = 0.0;
                curTerminal->isTerminalNI = true;
            }
            else {
                curTerminal->isTerminalNI = false;
            }

            if(flg_3dic) {
                curTerminal->size.z = TIER_DEP;
                curTerminal->area = curTerminal->size.x * curTerminal->size.y *
                                    curTerminal->size.z;
            }
            else {
                curTerminal->area = curTerminal->size.x * curTerminal->size.y;
            }

            // why this in here..
            curTerminal->pof = NULL;
            curTerminal->pin = NULL;
            curTerminal->netCNTinObject = 0;
            curTerminal->pinCNTinObject = 0;

            // max & min data update
            minTerm.x = (minTerm.x > curTerminal->size.x) ? curTerminal->size.x
                                                          : minTerm.x;
            minTerm.y = (minTerm.y > curTerminal->size.y) ? curTerminal->size.y
                                                          : minTerm.y;

            maxTerm.x = (maxTerm.x < curTerminal->size.x) ? curTerminal->size.x
                                                          : maxTerm.x;
            maxTerm.y = (maxTerm.y < curTerminal->size.y) ? curTerminal->size.y
                                                          : maxTerm.y;

            // update terminal_size_max (x,y)
            if(first_term) {
                first_term = 0;
                terminal_size_max = curTerminal->size;
            }
            else {
                terminal_size_max.x =
                    max(curTerminal->size.x, terminal_size_max.x);
                terminal_size_max.y =
                    max(curTerminal->size.y, terminal_size_max.y);
            }

            // average Terimanl Info update
            avgTerminalSize.x += curTerminal->size.x;
            avgTerminalSize.y += curTerminal->size.y;
        }
    }

    assert(cnt == totalInstanceCNT);

    // finally divide - module-cells
    avgCellSize.x /= (prec)moduleCNT;
    avgCellSize.y /= (prec)moduleCNT;
    avgCellSize.z = TIER_DEP;

    // finally divide - terminal
    avgTerminalSize.x /= (prec)terminalCNT;
    avgTerminalSize.y /= (prec)terminalCNT;
    avgTerminalSize.z = TIER_DEP;

    printf("INFO:  PLACEMENT INSTANCES INCLUDING STDCELL and MOVABLE MACRO\n");
    printf("INFO:    Instance   MinX=%.2lf, MaxX=%.2lf\n", minCell.x,
           maxCell.x);
    printf("INFO:               MinY=%.2lf, MaxY=%.2lf\n", minCell.y,
           maxCell.y);
    printf("INFO:               AvgX=%.2lf, AvgY=%.2lf\n", avgCellSize.x,
           avgCellSize.y);
    printf("INFO:      Smallest Instance  %10s  Size %.2f\n", minSizeCell->name,
           minSize);
    printf("INFO:      Largest  Instance  %10s  Size %.2f\n", maxSizeCell->name,
           maxSize);

    printf("INFO:  TERMINALS INCLUDING PAD and FIXED MACRO\n");
    printf("INFO:    Terminal   MinX=%.2lf, MaxX=%.2lf\n", minTerm.x,
           maxTerm.x);
    printf("INFO:               MinY=%.2lf, MaxY=%.2lf\n", minTerm.y,
           maxTerm.y);
    printf("INFO:               AvgX=%.2lf, AvgY=%.2lf\n", avgTerminalSize.x,
           avgTerminalSize.y);
    max_mac_dim.x = maxCell.x;
    max_mac_dim.y = maxCell.y;
    max_mac_dim.z = maxCell.z;
    return 1;
}

/*
int read_nets (char *input) {
  FILE*fp=fopen(input,"r");

  char*token=NULL;

  char buf[BUF_SZ]// ,buf2[BUF_SZ] ;
  char name[BUF_SZ];
  int buf_size=BUF_SZ;
  int pid=0;
  int moduleID=0;
  int lab=0;
  int io=0;

  int i=0,j=0;

  struct NET*tempNet=NULL;
  struct MODULE*curModule=NULL;
  struct TERM*curTerminal=NULL;
  struct PIN*pin=NULL;
  struct FPOS offset;

  int   max_net_deg = 0 ;

  do
    {
      fgets(buf,buf_size,fp);
      token=strtok(buf," \t\n") ;
    }
  while(!token || token[0]=='#' || !strcmp(token,"UCLA") );

  token = strtok(NULL," \t\n");
  token = strtok(NULL," \t\n");

  netCNT = atoi(token);
  netInstance = (struct NET*)mkl_malloc(sizeof(struct NET)*netCNT, 64);

  do
    {
      fgets(buf,buf_size,fp);
      token=strtok(buf," \t\n") ;
    }
  while(!token || token[0]=='#' || !strcmp(token,"UCLA") );
  token = strtok(NULL," \t\n");
  token = strtok(NULL," \t\n");

  pinCNT = atoi(token);

  pinInstance = (struct PIN*)mkl_malloc(sizeof(struct PIN)*pinCNT, 64);


  for(i=0;i<netCNT;i++)
    {
      tempNet = & netInstance[i];

      tempNet->idx = i ;

      do
        {
          fgets(buf,buf_size,fp);
          token=strtok(buf," \t\n") ;
        }
      while(!token || token[0]=='#' || !strcmp(token,"UCLA") );

      token=strtok(NULL," \t\n");
      token=strtok(NULL," \t\n");
      tempNet->pinCNTinObject = atoi(token);

      if(max_net_deg < tempNet->pinCNTinObject)
        max_net_deg = tempNet->pinCNTinObject ;

      token=strtok(NULL," \t\n");
      strcpy(tempNet->name,token);

      tempNet->pin = (struct PIN**)mkl_malloc(sizeof(struct
PIN*)*tempNet->pinCNTinObject, 64);

      for(j=0;j<tempNet->pinCNTinObject;j++)
        {
          fgets(buf,buf_size,fp);
          token=strtok(buf," \t\n");

          strcpy(name,token);
          name2id(name,&moduleID,&lab);

          if(moduleID < 0 || moduleID > moduleCNT-1)
            {
              g_rrr ++ ;
              printf("ERROR ERROR ERROR!\n");
              exit(1);
            }

          token=strtok(NULL," \t\n");
          if(token[0]=='I')
            io = 0 ;
          else if(token[0]=='O')
            io = 1 ;
          else // 'B'
            io = 2 ;

          token=strtok(NULL," \t\n"); // ":"
          token=strtok(NULL," \t\n"); // "offset_x"
          offset.x=atof(token);
          token=strtok(NULL," \t\n"); // "offset_y"
          offset.y=atof(token);

          offset.z = 0.0;

          pin = & pinInstance[pid];

          tempNet->pin[j] = pin ;

          if(lab==0)
            {
              curModule = & moduleInstance[moduleID];

              AddPinInfoForModuleAndTerminal ( & curModule->pin ,
                       & curModule->pof ,
                       curModule->pinCNTinObject++,
                       offset        ,
                       moduleID        ,
                       i          ,
                       j          ,
                       pid ++     ,
                       io         ,
                       lab        );
            }
          else
            {
              curTerminal = & terminalInstance[moduleID] ;
              AddPinInfoForModuleAndTerminal (
&curTerminal->pin,&curTerminal->pof,curTerminal->pinCNTinObject++,offset,moduleID,i,j,pid++,io,lab
);
            }

        }
    }

  printf("MAX_NET_DEG = %d\n",max_net_deg);
  printf("# NETS = %d\n",netCNT);

  return 1;
}
*/

int read_nets_3D(char *input) {
    FILE *fp = fopen(input, "r");
    char *token = NULL;
    char buf[255];
    char name[255];
    int buf_size = 255;
    int moduleID = 0;
    int io = 0;
    NET *curNet = NULL;
    MODULE *curModule = NULL;
    TERM *curTerminal = NULL;
    PIN *pin = NULL;
    FPOS offset;

    int isTerminal = 0;
    int max_net_deg = INT_MIN;

    do {
        fgets(buf, buf_size, fp);
        token = strtok(buf, " \t\n");
    } while(!token || token[0] == '#' || !strcmp(token, "UCLA"));
    token = strtok(NULL, " \t\n");
    token = strtok(NULL, " \t\n");
    netCNT = atoi(token);

    do {
        fgets(buf, buf_size, fp);
        token = strtok(buf, " \t\n");
    } while(!token || token[0] == '#' || !strcmp(token, "UCLA"));
    token = strtok(NULL, " \t\n");
    token = strtok(NULL, " \t\n");
    pinCNT = atoi(token);

    netInstance = (struct NET *)mkl_malloc(sizeof(struct NET) * netCNT, 64);
    pinInstance = (struct PIN *)mkl_malloc(sizeof(struct PIN) * pinCNT, 64);

    int pid = 0;
    for(int i = 0; i < netCNT; i++) {
        curNet = &netInstance[i];
        curNet->idx = i;

        do {
            fgets(buf, buf_size, fp);
            token = strtok(buf, " \t\n");
        } while(!token || strcmp(token, "NetDegree"));

        token = strtok(NULL, " \t\n");
        token = strtok(NULL, " \t\n");
        curNet->pinCNTinObject = atoi(token);

        token = strtok(NULL, " \t\n");

        if(token)
            strcpy(curNet->name, token);
        // lutong 05132016
        // if (strcmp(curNet->name, "CLK") && (strcmp(curNet->name, "GND"))){

        curNet->pin = (struct PIN **)mkl_malloc(
            sizeof(struct PIN *) * curNet->pinCNTinObject, 64);

        if(max_net_deg < curNet->pinCNTinObject) {
            max_net_deg = curNet->pinCNTinObject;
        }

        for(int j = 0; j < curNet->pinCNTinObject; j++) {
            fgets(buf, buf_size, fp);
            token = strtok(buf, " \t\n");
            strcpy(name, token);

            auto nodesMapPtr = nodesMap.find(string(name));
            if(nodesMapPtr == nodesMap.end()) {
                runtimeError(
                    string(string(".net files' nodename is mismatched! : ") +
                           string(name))
                        .c_str());
            }
            else {
                isTerminal = nodesMapPtr->second.isTerminal;
                moduleID = nodesMapPtr->second.index;
                // cout << "currentname: " << name << ", index: "
                // << nodesMapPtr->second.index << endl;
            }

            token = strtok(NULL, " \t\n");
            if(token[0] == 'I')
                io = 0;  // 'I'
            else if(token[0] == 'O')
                io = 1;  // 'O'
            else
                io = 2;                     // 'B'
            token = strtok(NULL, " \t\n");  // ":"

            if(token != 0) {
                token = strtok(NULL, " \t\n");  // "offset_x"
                offset.x = atof(token);

                token = strtok(NULL, " \t\n");  // "offset_y"
                offset.y = atof(token);
                offset.z = 0.0;
            }
            else {
                offset.x = 0.0;
                offset.y = 0.0;
                offset.z = 0.0;
            }

            pin = &pinInstance[pid];
            curNet->pin[j] = pin;

            // because we don't know, in the same nets,
            // what component goes into moduleInstance's pin or
            // terminalInstance's pin,
            //
            // It always reallocate all net's informations

            if(isTerminal == false) {
                curModule = &moduleInstance[moduleID];
                AddPinInfoForModuleAndTerminal(&curModule->pin, &curModule->pof,
                                               curModule->pinCNTinObject++,
                                               offset, moduleID, i, j, pid++,
                                               io, isTerminal);
            }
            else {
                curTerminal = &terminalInstance[moduleID];
                AddPinInfoForModuleAndTerminal(
                    &curTerminal->pin, &curTerminal->pof,
                    curTerminal->pinCNTinObject++, offset, moduleID, i, j,
                    pid++, io, isTerminal);
            }
        }
        //} else {
        // i      -= 1;
        // netCNT -= 1;
        // pinCNT -= curNet->pinCNTinObject;
        //}
    }

    printf("INFO:  #NET=%d\n", netCNT);
    printf("INFO:  #PIN=%d\n", pinCNT);
    printf("INFO:    Maximum Net Degree is %d\n", max_net_deg);
    return 1;
}

/*
int read_pl (char *input) {
    FILE *fp = fopen(input, "r");
    char *token = NULL;
    char buf[BUF_SZ];
    char name[BUF_SZ];
    int buf_size=BUF_SZ;
    int lab=0;
    int i=0;
    int moduleID=0;
    struct MODULE*curModule// ,*tmp=NULL ;
    struct TERM*curTerminal=NULL;

    for (i=0; i<moduleCNT; i++) {
        do {
            fgets (buf, buf_size, fp);
            token = strtok(buf, " \t\n");
        }
        while (!token || token[0]=='#' || !strcmp(token,"UCLA"));

        strcpy (name, token);
        name2id(name, &moduleID, &lab);

        if (lab!=0) {
            printf("read_pl: error 1\n");
            g_rrr++;
        }

        curModule = & moduleInstance[moduleID] ;

        token = strtok(NULL," \t\n");
        curModule->pmin.x = atof(token);

        token = strtok(NULL," \t\n");
        curModule->pmin.y = atof(token);

        curModule->pmin.z = 0.0;

        curModule->center.x = (prec) curModule->pmin.x + (prec)
curModule->size.x * 0.5 ;
        curModule->center.y = (prec) curModule->pmin.y + (prec)
curModule->size.y * 0.5 ;
        curModule->center.z = (prec) curModule->pmin.z + (prec)
curModule->size.z * 0.5 ;

        curModule->pmax.x = curModule->pmin.x + curModule->size.x ;
        curModule->pmax.y = curModule->pmin.y + curModule->size.y ;
        curModule->pmax.z = curModule->pmin.z + curModule->size.z ;
    }



    for(i=0;i<terminalCNT;i++)
    {
        do
        {
            fgets(buf,buf_size,fp);
            token=strtok(buf," \t\n") ;
        }
        while(!token || token[0]=='#' || !strcmp(token,"UCLA") );

        strcpy (name,token);
        name2id(name,&moduleID,&lab);
        curTerminal = & terminalInstance[moduleID] ;

        token = strtok(NULL," \t\n");
        curTerminal->pmin.x = (prec) atoi(token);

        token = strtok(NULL," \t\n");
        curTerminal->pmin.y = (prec) atoi(token);

        curTerminal->pmin.z = 0.0;

        curTerminal->center.x = (prec) curTerminal->pmin.x + (prec)
curTerminal->size.x * 0.5 ;
        curTerminal->center.y = (prec) curTerminal->pmin.y + (prec)
curTerminal->size.y * 0.5 ;
        curTerminal->center.z = (prec) curTerminal->pmin.z + (prec)
curTerminal->size.z * 0.5 ;

        curTerminal->pmax.x = curTerminal->pmin.x + curTerminal->size.x ;
        curTerminal->pmax.y = curTerminal->pmin.y + curTerminal->size.y ;
        curTerminal->pmax.z = curTerminal->pmin.z + curTerminal->size.z ;
    }

    return 1;
}

*/


// 
// update terminalInst & moduleInst's Info
//
int read_pl2(char *input) {
    FILE *fp = fopen(input, "r");
    char *token = NULL;
    char buf[255];
    int buf_size = 255;
    int moduleID = 0;
    char name[255];
    struct MODULE *curModule = NULL;
    struct TERM *curTerminal = NULL;

    bool isFirst = true;

    for(int i = 0; i < moduleCNT + terminalCNT; i++) {
        do {
            fgets(buf, buf_size, fp);
            token = strtok(buf, " \t\n");
        } while(!token || token[0] == '#' || !strcmp(token, "UCLA"));

        strcpy(name, token);
        // getCellIndex ( name , & moduleID , & isTerminal , &isTerminalNI)
        // ;

        bool isTerminal = false;
        auto nodesMapPtr = nodesMap.find(string(name));
        if(nodesMapPtr == nodesMap.end()) {
            runtimeError(
                string(string(".pl files' nodename is mismatched! : ") +
                       string(name))
                    .c_str());
        }
        else {
            isTerminal = nodesMapPtr->second.isTerminal;
            moduleID = nodesMapPtr->second.index;
            // cout << "currentname: " << name << ", index: " <<
            // nodesMapPtr->second.index << endl;
        }

        if(!isTerminal) {
            curModule = &moduleInstance[moduleID];
            token = strtok(NULL, " \t\n");
            curModule->pmin.x = atof(token);

            token = strtok(NULL, " \t\n");
            curModule->pmin.y = atof(token);

            if(flg_3dic) {
                if(flg_3dic_io) {
                    token = strtok(NULL, " \t\n");
                    curModule->pmin.z = atof(token);
                }
                else {
                    curModule->pmin.z = 0;
                }
            }

            curModule->center.x = curModule->pmin.x + curModule->half_size.x;
            curModule->center.y = curModule->pmin.y + curModule->half_size.y;

            if(flg_3dic) {
                curModule->center.z =
                    curModule->pmin.z + curModule->half_size.z;
            }

            curModule->pmax.x = curModule->pmin.x + curModule->size.x;
            curModule->pmax.y = curModule->pmin.y + curModule->size.y;

            if(flg_3dic) {
                curModule->pmax.z = curModule->pmin.z + curModule->size.z;
            }
        }
        else {
            curTerminal = &terminalInstance[moduleID];

            token = strtok(NULL, " \t\n");
            curTerminal->pmin.x = atof(token);

            token = strtok(NULL, " \t\n");
            curTerminal->pmin.y = atof(token);

            if(flg_3dic) {
                if(flg_3dic_io) {
                    token = strtok(NULL, " \t\n");
                    curTerminal->pmin.z = atof(token);
                }
                else {
                    curTerminal->pmin.z = 0;
                }
            }

            curTerminal->center.x =
                curTerminal->pmin.x + 0.5 * curTerminal->size.x;
            curTerminal->center.y =
                curTerminal->pmin.y + 0.5 * curTerminal->size.y;
            if(flg_3dic) {
                curTerminal->center.z =
                    curTerminal->pmin.z + 0.5 * curTerminal->size.z;
            }

            curTerminal->pmax.x = curTerminal->pmin.x + curTerminal->size.x;
            curTerminal->pmax.y = curTerminal->pmin.y + curTerminal->size.y;
            if(flg_3dic) {
                curTerminal->pmax.z = curTerminal->pmin.z + curTerminal->size.z;
            }

            if(isFirst) {
                terminal_pmin = curTerminal->pmin;
                terminal_pmax = curTerminal->pmax;

                isFirst = false;
            }
            else {
                terminal_pmin.x = min(terminal_pmin.x, curTerminal->pmin.x);
                terminal_pmin.y = min(terminal_pmin.y, curTerminal->pmin.y);

                terminal_pmax.x = max(terminal_pmax.x, curTerminal->pmax.x);
                terminal_pmax.y = max(terminal_pmax.y, curTerminal->pmax.y);
            }
        }
    }
    return 1;
}

// 
// build the row_st(ROW)
//
int read_scl(char *input) {
    char *token = NULL;
    char buf[BUF_SZ];
    int buf_size = BUF_SZ;

    FILE *fp = fopen(input, "r");
    do {
        fgets(buf, buf_size, fp);
        token = strtok(buf, " \t\n");
    } while(!token || token[0] == '#' || !strcmp(token, "UCLA"));

    token = strtok(NULL, " \t\n");
    token = strtok(NULL, " \t\n");

    row_cnt = atoi(token);
//    row_st = (struct ROW *)malloc(sizeof(struct ROW) * row_cnt);
    // call constructor
//    row_st = new ROW[row_cnt];
    row_st = (ROW*) mkl_malloc( sizeof(struct ROW)* row_cnt, 64 );
    fgets(buf, buf_size, fp);

    for(int i = 0; i < row_cnt; i++) {
        ROW* row = &row_st[i];
        new (row) ROW();

        fgets(buf, buf_size, fp);

        fgets(buf, buf_size, fp);
        token = strtok(buf, " \t\n");
        token = strtok(NULL, " \t\n");
        token = strtok(NULL, " \t\n");
        row->pmin.y = atoi(token);  // Coordinate

        fgets(buf, buf_size, fp);
        token = strtok(buf, " \t\n");
        token = strtok(NULL, " \t\n");
        token = strtok(NULL, " \t\n");
        row->size.y = atoi(token);  // Height
        row->pmax.y = row->pmin.y + row->size.y;

        fgets(buf, buf_size, fp);
        token = strtok(buf, " \t\n");
        token = strtok(NULL, " \t\n");
        token = strtok(NULL, " \t\n");
        row->site_wid = atoi(token);  // Sitewidth

        fgets(buf, buf_size, fp);
        token = strtok(buf, " \t\n");
        token = strtok(NULL, " \t\n");
        token = strtok(NULL, " \t\n");
        row->site_spa = atoi(token);  // Sitespacing

        fgets(buf, buf_size, fp);
        token = strtok(buf, " \t\n");
        token = strtok(NULL, " \t\n");
        token = strtok(NULL, " \t\n");
//        row->ori = atoi(token);  // Siteorient
        row->ori = string(token);  // Siteorient

        fgets(buf, buf_size, fp);
        token = strtok(buf, " \t\n");
        token = strtok(NULL, " \t\n");
        token = strtok(NULL, " \t\n");
//        row->sym = atoi(token);  // Sitesymmetry
//        row->sym = string(token);  // Sitesymmetry
        
        if( string(token) == "Y" ) {
            row->isYSymmetry = true;
        }

        fgets(buf, buf_size, fp);
        token = strtok(buf, " \t\n");
        token = strtok(NULL, " \t\n");
        token = strtok(NULL, " \t\n");
        row->pmin.x = atoi(token);  // SubrowOrigin

        token = strtok(NULL, " \t\n");
        token = strtok(NULL, " \t\n");
        token = strtok(NULL, " \t\n");
        row->x_cnt = atoi(token);  // NumSites

        row->size.x = row->x_cnt * row->site_spa;
        row->pmax.x = row->pmin.x + row->size.x;

        if(i == 0) {
            grow_pmin.Set(row->pmin);
            grow_pmax.Set(row->pmax);
        }
        else {
            grow_pmin.x = min(grow_pmin.x, (prec)row->pmin.x);
            grow_pmin.y = min(grow_pmin.y, (prec)row->pmin.y);
            grow_pmax.x = max(grow_pmax.x, (prec)row->pmax.x);
            grow_pmax.y = max(grow_pmax.y, (prec)row->pmax.y);
        }
        
        // For Z handling
        row->size.z = TIER_DEP;
        row->pmin.z = 0;
        row->pmax.z = TIER_DEP;

        fgets(buf, buf_size, fp);

        if(i == 0) {
            rowHeight = row->size.y;
            SITE_SPA = row->site_spa;
        }
        else if(rowHeight != row->size.y) {
            printf("Error: ROW HEIGHT INCONSISTENT!\n");
            exit(0);
        }
        else if(SITE_SPA != row->site_spa) {
            printf("Error: SITE SPACING INCONSISTENT!\n");
            exit(0);
        }
    }
    
    // check ROW instance
//    for(int i=0; i<row_cnt; i++) {
//        ROW* curRow = &row_st[i];
//        curRow->Dump(to_string(i));
//    }
//
    return 1;
}

/*

void name2id (char *name, int *idx,int *lab) {
    char    a[BUF_SZ];
    strcpy (a, name+1);
    int     b = atoi(a);
    if (b < moduleCNT) {
        *idx = modu_map[b];
        *lab = 0;
    } else {
        *idx = term_map[b];
        *lab = 1;
    }
}
*/

// int read_wts() {
// return  1;
//}

void output_mGP2D_pl(char *output) {
    FILE *fp = fopen(output, "w");
    struct MODULE *curModule = NULL;
    struct TERM *curTerminal = NULL;

    fputs("UCLA pl 1.0 \n", fp);
    fputs("# Created	:	Jan  6 2005\n", fp);
    fputs(
        "# User   	:	Gi-Joon Nam & Mehmet Yildiz at IBM Austin "
        "Research({gnam, mcan}@us.ibm.com)\n",
        fp);
    fputs("\n", fp);

    for(int i = 0; i < moduleCNT; i++) {
        curModule = &moduleInstance[i];
        if(flg_3dic) {
            fprintf(fp, "%s\t%.6lf\t%.6lf\t%.6lf\t: N\n", curModule->name,
                    curModule->pmin.x, curModule->pmin.y, curModule->pmin.z);
        }
        else {
            fprintf(fp, "%s\t%.6lf\t%.6lf\t: N\n", curModule->name,
                    curModule->pmin.x, curModule->pmin.y);
        }
    }

    for(int i = 0; i < terminalCNT; i++) {
        if(flg_3dic) {
            curTerminal = &terminalInstance[i];
            fprintf(fp, "%s %d\t%d\t%d\t: N /FIXED\n", curTerminal->name,
                    prec2int(curTerminal->pmin.x),
                    prec2int(curTerminal->pmin.y),
                    prec2int(curTerminal->pmin.z));
        }
        else {
            curTerminal = &terminalInstance[i];
            fprintf(fp, "%s %d\t%d\t: N /FIXED\n", curTerminal->name,
                    prec2int(curTerminal->pmin.x),
                    prec2int(curTerminal->pmin.y));
        }
    }
    fclose(fp);
}

// write pl files
void output_cGP2D_pl(char *output) {
    int i = 0;
    FILE *fp = fopen(output, "w");
    MODULE *curModule = NULL;
    TERM *curTerminal = NULL;

    fputs("UCLA pl 1.0 \n", fp);
    fputs("# Created	:	Jan  6 2005\n", fp);
    fputs(
        "# User   	:	Gi-Joon Nam & Mehmet Yildiz at IBM Austin "
        "Research({gnam, mcan}@us.ibm.com)\n",
        fp);
    fputs("\n", fp);

    for(i = 0; i < moduleCNT; i++) {
        curModule = &moduleInstance[i];

        if(flg_3dic) {
            if(curModule->flg == StdCell) {
                fprintf(fp, "%s\t%.6lf\t%.6lf\t%.6lf\t: N\n", curModule->name,
                        curModule->pmin.x, curModule->pmin.y,
                        curModule->pmin.z);
            }
            else {
                fprintf(fp, "%s\t%d\t%d\t%d\t: N /FIXED\n", curModule->name,
                        curModule->pmin_lg.x, curModule->pmin_lg.y,
                        curModule->pmin_lg.z);
            }
        }
        else {
            if(curModule->flg == StdCell) {
                fprintf(fp, "%s\t%.6lf\t%.6lf\t: N\n", curModule->name,
                        curModule->pmin.x, curModule->pmin.y);
            }
            else {
                fprintf(fp, "%s\t%d\t%d\t: N /FIXED\n", curModule->name,
                        curModule->pmin_lg.x, curModule->pmin_lg.y);
            }
        }
    }

    for(i = 0; i < terminalCNT; i++) {
        if(flg_3dic) {
            curTerminal = &terminalInstance[i];
            fprintf(fp, "%s %d\t%d\t%d\t: N /FIXED\n", curTerminal->name,
                    prec2int(curTerminal->pmin.x),
                    prec2int(curTerminal->pmin.y),
                    prec2int(curTerminal->pmin.z));
        }
        else {
            curTerminal = &terminalInstance[i];
            fprintf(fp, "%s %d\t%d\t: N /FIXED\n", curTerminal->name,
                    prec2int(curTerminal->pmin.x),
                    prec2int(curTerminal->pmin.y));
        }
    }

    fclose(fp);
}

void output_cGP3D_pl(char *output) {
    int i = 0;
    FILE *fp = fopen(output, "w");
    struct MODULE *curModule = NULL;
    struct TERM *curTerminal = NULL;

    fputs("UCLA pl 1.0 \n", fp);
    fputs("# Created	:	Jan  6 2005\n", fp);
    fputs(
        "# User   	:	Gi-Joon Nam & Mehmet Yildiz at IBM Austin "
        "Research({gnam, mcan}@us.ibm.com)\n",
        fp);
    fputs("\n", fp);

    for(i = 0; i < moduleCNT; i++) {
        curModule = &moduleInstance[i];

        if(flg_3dic) {
            if(curModule->flg == StdCell) {
                fprintf(fp, "%s\t%.6lf\t%.6lf\t%.6lf\t: N\n", curModule->name,
                        curModule->pmin.x, curModule->pmin.y,
                        curModule->pmin.z);
            }
            else {
                fprintf(fp, "%s\t%d\t%d\t%d\t: N /FIXED\n", curModule->name,
                        curModule->pmin_lg.x, curModule->pmin_lg.y,
                        curModule->pmin_lg.z);
            }
        }
        else {
            if(curModule->flg == StdCell) {
                fprintf(fp, "%s\t%.6lf\t%.6lf\t: N\n", curModule->name,
                        curModule->pmin.x, curModule->pmin.y);
            }
            else {
                fprintf(fp, "%s\t%d\t%d\t: N /FIXED\n", curModule->name,
                        curModule->pmin_lg.x, curModule->pmin_lg.y);
            }
        }
    }

    for(i = 0; i < terminalCNT; i++) {
        if(flg_3dic) {
            curTerminal = &terminalInstance[i];
            fprintf(fp, "%s %d\t%d\t%d\t: N /FIXED\n", curTerminal->name,
                    prec2int(curTerminal->pmin.x),
                    prec2int(curTerminal->pmin.y),
                    prec2int(curTerminal->pmin.z));
        }
        else {
            curTerminal = &terminalInstance[i];
            fprintf(fp, "%s %d\t%d\t: N /FIXED\n", curTerminal->name,
                    prec2int(curTerminal->pmin.x),
                    prec2int(curTerminal->pmin.y));
        }
    }

    fclose(fp);
}

/*
void gp_sol_trans_to_dim_one(void) {
    int i = 0, j = 0;
    MODULE *modu = NULL;
    TERM *curTerminal = NULL;

    for(i = 0; i < moduleCNT; i++) {
        modu = &moduleInstance[i];
        modu->center.x = (modu->center.x - place_backup.org.x) /
                         (place_backup.cnt.x * GP_SCAL);
        modu->center.y = (modu->center.y - place_backup.org.y) /
                         (place_backup.cnt.y * GP_SCAL);
        modu->center.z = (modu->center.z - place_backup.org.z) /
                         (place_backup.cnt.z * GP_SCAL);
        modu->size.x /= place_backup.cnt.x * GP_SCAL;
        modu->size.y /= place_backup.cnt.y * GP_SCAL;
        modu->size.z /= place_backup.cnt.z * GP_SCAL;
        modu->half_size.x = 0.5 * modu->size.x;
        modu->half_size.y = 0.5 * modu->size.y;
        modu->half_size.z = 0.5 * modu->size.z;
        modu->pmin.x = modu->center.x - modu->half_size.x;
        modu->pmin.y = modu->center.y - modu->half_size.y;
        modu->pmin.z = modu->center.z - modu->half_size.z;
        modu->pmax.x = modu->center.x + modu->half_size.x;
        modu->pmax.y = modu->center.y + modu->half_size.y;
        modu->pmax.z = modu->center.z + modu->half_size.z;
        modu->area = modu->size.x * modu->size.y * modu->size.z;

        for(j = 0; j < modu->pinCNTinObject; j++) {
            modu->pof[j].x /= place_backup.cnt.x * GP_SCAL;
            modu->pof[j].y /= place_backup.cnt.y * GP_SCAL;
            modu->pof[j].z /= place_backup.cnt.z * GP_SCAL;
        }
    }

    for(i = 0; i < terminalCNT; i++) {
        curTerminal = &terminalInstance[i];
        curTerminal->center.x = (curTerminal->center.x - place_backup.org.x) /
                                (place_backup.cnt.x * GP_SCAL);
        curTerminal->center.y = (curTerminal->center.y - place_backup.org.y) /
                                (place_backup.cnt.y * GP_SCAL);
        curTerminal->center.z = (curTerminal->center.z - place_backup.org.z) /
                                (place_backup.cnt.z * GP_SCAL);
        curTerminal->size.x /= place_backup.cnt.x * GP_SCAL;
        curTerminal->size.y /= place_backup.cnt.y * GP_SCAL;
        curTerminal->size.z /= place_backup.cnt.z * GP_SCAL;
        curTerminal->pmin.x = curTerminal->center.x - 0.5 * curTerminal->size.x;
        curTerminal->pmin.y = curTerminal->center.y - 0.5 * curTerminal->size.y;
        curTerminal->pmin.z = curTerminal->center.z - 0.5 * curTerminal->size.z;
        curTerminal->pmax.x = curTerminal->center.x + 0.5 * curTerminal->size.x;
        curTerminal->pmax.y = curTerminal->center.y + 0.5 * curTerminal->size.y;
        curTerminal->pmax.z = curTerminal->center.z + 0.5 * curTerminal->size.z;
        curTerminal->area =
            curTerminal->size.x * curTerminal->size.y * curTerminal->size.z;

        for(j = 0; j < curTerminal->pinCNTinObject; j++) {
            curTerminal->pof[j].x /= place_backup.cnt.x * GP_SCAL;
            curTerminal->pof[j].y /= place_backup.cnt.y * GP_SCAL;
            curTerminal->pof[j].z /= place_backup.cnt.z * GP_SCAL;
        }
    }
}

void gp_sol_trans_from_dim_one(void) {
    int i = 0, j = 0;
    MODULE *modu = NULL;
    TERM *curTerminal = NULL;

    for(i = 0; i < moduleCNT; i++) {
        modu = &moduleInstance[i];
        modu->center.x = modu->center.x * (place_backup.cnt.x * GP_SCAL) +
                         place_backup.org.x;
        modu->center.y = modu->center.y * (place_backup.cnt.y * GP_SCAL) +
                         place_backup.org.y;
        modu->center.z = modu->center.z * (place_backup.cnt.z * GP_SCAL) +
                         place_backup.org.z;
        modu->size.x *= place_backup.cnt.x * GP_SCAL;
        modu->size.y *= place_backup.cnt.y * GP_SCAL;
        modu->size.z *= place_backup.cnt.z * GP_SCAL;
        modu->half_size.x = 0.5 * modu->size.x;
        modu->half_size.y = 0.5 * modu->size.y;
        modu->half_size.z = 0.5 * modu->size.z;
        modu->pmin.x = modu->center.x - modu->half_size.x;
        modu->pmin.y = modu->center.y - modu->half_size.y;
        modu->pmin.z = modu->center.z - modu->half_size.z;
        modu->pmax.x = modu->center.x + modu->half_size.x;
        modu->pmax.y = modu->center.y + modu->half_size.y;
        modu->pmax.z = modu->center.z + modu->half_size.z;
        modu->area = modu->size.x * modu->size.y * modu->size.z;

        for(j = 0; j < modu->pinCNTinObject; j++) {
            modu->pof[j].x *= place_backup.cnt.x * GP_SCAL;
            modu->pof[j].y *= place_backup.cnt.y * GP_SCAL;
            modu->pof[j].z *= place_backup.cnt.z * GP_SCAL;
        }
    }

    for(i = 0; i < terminalCNT; i++) {
        curTerminal = &terminalInstance[i];
        curTerminal->center.x =
            curTerminal->center.x * (place_backup.cnt.x * GP_SCAL) +
            place_backup.org.x;
        curTerminal->center.y =
            curTerminal->center.y * (place_backup.cnt.y * GP_SCAL) +
            place_backup.org.y;
        curTerminal->center.z =
            curTerminal->center.z * (place_backup.cnt.z * GP_SCAL) +
            place_backup.org.z;
        curTerminal->size.x *= place_backup.cnt.x * GP_SCAL;
        curTerminal->size.y *= place_backup.cnt.y * GP_SCAL;
        curTerminal->size.z *= place_backup.cnt.z * GP_SCAL;
        curTerminal->pmin.x = curTerminal->center.x - 0.5 * curTerminal->size.x;
        curTerminal->pmin.y = curTerminal->center.y - 0.5 * curTerminal->size.y;
        curTerminal->pmin.z = curTerminal->center.z - 0.5 * curTerminal->size.z;
        curTerminal->pmax.x = curTerminal->center.x + 0.5 * curTerminal->size.x;
        curTerminal->pmax.y = curTerminal->center.y + 0.5 * curTerminal->size.y;
        curTerminal->pmax.z = curTerminal->center.z + 0.5 * curTerminal->size.z;
        curTerminal->area =
            curTerminal->size.x * curTerminal->size.y * curTerminal->size.z;

        for(j = 0; j < curTerminal->pinCNTinObject; j++) {
            curTerminal->pof[j].x *= place_backup.cnt.x * GP_SCAL;
            curTerminal->pof[j].y *= place_backup.cnt.y * GP_SCAL;
            curTerminal->pof[j].z *= place_backup.cnt.z * GP_SCAL;
        }
    }
}
*/

void output_pl(char *output) {
    FILE *fp = fopen(output, "w");
    struct MODULE *curModule = NULL;
    struct TERM *curTerminal = NULL;

    fputs("UCLA pl 1.0 \n", fp);
    fputs("# Created	:	Jan  6 2005\n", fp);
    fputs(
        "# User   	:	Gi-Joon Nam & Mehmet Yildiz at IBM Austin "
        "Research({gnam, mcan}@us.ibm.com)\n",
        fp);
    fputs("\n", fp);

    if(flg_3dic) {
        for(int i = 0; i < moduleCNT; i++) {
            curModule = &moduleInstance[i];
            fprintf(fp, "%s %.6lf\t%.6lf\t%.6lf\t: N\n", curModule->name,
                    curModule->pmin.x, curModule->pmin.y, curModule->pmin.z);
        }
        for(int i = 0; i < terminalCNT; i++) {
            curTerminal = &terminalInstance[i];
            // if (!curTerminal->isTerminalNI) {
            fprintf(fp, "%s %d\t%d\t%d\t: N /FIXED\n", curTerminal->name,
                    prec2int(curTerminal->pmin.x),
                    prec2int(curTerminal->pmin.y),
                    prec2int(curTerminal->pmin.z));
            //} else {
            // fprintf (fp, "%s %d\t%d\t%d\t: N /FIXED_NI\n",
            // curTerminal->name,
            // prec2int (curTerminal->pmin.x),
            // prec2int (curTerminal->pmin.y),
            // prec2int (curTerminal->pmin.z));
            //}
        }
    }
    else {
        for(int i = 0; i < moduleCNT; i++) {
            curModule = &moduleInstance[i];
            fprintf(fp, "%s %.6lf\t%.6lf\t: N\n", curModule->name,
                    curModule->pmin.x, curModule->pmin.y);
        }
        for(int i = 0; i < terminalCNT; i++) {
            curTerminal = &terminalInstance[i];
            // if (!curTerminal->isTerminalNI) {
            fprintf(fp, "%s %d\t%d\t: N /FIXED\n", curTerminal->name,
                    prec2int(curTerminal->pmin.x),
                    prec2int(curTerminal->pmin.y));
            //} else {
            // fprintf (fp, "%s %d\t%d\t: N /FIXED_NI\n",
            // curTerminal->name,
            // prec2int (curTerminal->pmin.x),
            // prec2int (curTerminal->pmin.y));
            //}
        }
    }
    fclose(fp);
}

// function is modified by mgwoo.
//
// this function will
// add pin's informations for moduleInstance & terminalInstance
// once it meets fixed or unfixed pins..
void AddPinInfoForModuleAndTerminal(PIN ***pin, FPOS **pof,
                                           int currentPinCount, FPOS offset,
                                           int curModuleIdx, int curNetIdx,
                                           int curPinIdxInNet, int curPinIdx,
                                           int curPinDirection,
                                           int isTerminal) {
    // pin : current Pin Object
    // pof : current Pin Offset
    //
    if(currentPinCount > 0) {
        currentPinCount++;
        *pin = (PIN **)realloc(*pin, sizeof(struct PIN *) * currentPinCount);
        *pof = (FPOS *)realloc(*pof, sizeof(struct FPOS) * currentPinCount);
        //*pin_tmp = (struct PIN**)mkl_malloc(sizeof(struct
        // PIN*)*currentPinCount, 64);
        //*pof_tmp = (struct FPOS*)mkl_malloc(sizeof(struct
        // FPOS)*currentPinCount, 64);
        // memcpy(*pin_tmp, *pin, currentPinCount*(sizeof(struct PIN*)));
        // memcpy(*pof_tmp, *pof, currentPinCount*(sizeof(struct FPOS)));
        // mkl_free(*pin);
        // mkl_free(*pof);
        //*pin = (struct PIN**)mkl_malloc(sizeof(struct PIN*)*currentPinCount,
        // 64);
        //*pof = (struct FPOS*)mkl_malloc(sizeof(struct FPOS)*currentPinCount,
        // 64);
        // memcpy(*pin, *pin_tmp, currentPinCount*(sizeof(struct PIN*)));
        // memcpy(*pof, *pof_tmp, currentPinCount*(sizeof(struct FPOS)));
        // mkl_free(*pin_tmp);
        // mkl_free(*pof_tmp);
    }
    else {
        currentPinCount++;
        *pin = (struct PIN **)malloc(sizeof(struct PIN *) * currentPinCount);
        *pof = (struct FPOS *)malloc(sizeof(struct FPOS) * currentPinCount);
    }

    (*pof)[currentPinCount - 1] = offset;

    (*pin)[currentPinCount - 1] = &pinInstance[curPinIdx];
    (*pin)[currentPinCount - 1]->term = isTerminal;
    (*pin)[currentPinCount - 1]->moduleID = curModuleIdx;
    (*pin)[currentPinCount - 1]->netID = curNetIdx;
    (*pin)[currentPinCount - 1]->pinIDinNet = curPinIdxInNet;
    (*pin)[currentPinCount - 1]->pinIDinModule = currentPinCount - 1;
    (*pin)[currentPinCount - 1]->gid = curPinIdx;
    (*pin)[currentPinCount - 1]->IO = curPinDirection;
}

void write_new_bench(void) {
    /////////////////////////////////

    int i = 0, j = 0, moduleID = 0, pinIDinModule = 0;

    sprintf(sing_fn_aux, "%s.aux", gbch);
    sprintf(sing_fn_nets, "%s.nets", gbch);
    sprintf(sing_fn_nodes, "%s.nodes", gbch);
    sprintf(sing_fn_pl, "%s.pl", gbch);
    sprintf(sing_fn_wts, "%s.wts", gbch);
    sprintf(sing_fn_scl, "%s.scl", gbch);

    sprintf(gTMP_bch_dir, "%s/%s", dir_bnd, gbch);

    // mkdir -p ../output/SB/experment0/tier/0
    char mkdir_cmd[BUF_SZ];
    sprintf(mkdir_cmd, "mkdir -p %s", gTMP_bch_dir);
    system(mkdir_cmd);

    sprintf(gTMP_bch_aux, "%s/%s.aux", gTMP_bch_dir, gbch);

    sprintf(gTMP_bch_nets, "%s/%s.nets", gTMP_bch_dir, gbch);

    sprintf(gTMP_bch_nodes, "%s/%s.nodes", gTMP_bch_dir, gbch);

    sprintf(gTMP_bch_pl, "%s/%s.pl", gTMP_bch_dir, gbch);

    sprintf(gTMP_bch_wts, "%s/%s.wts", gTMP_bch_dir, gbch);

    sprintf(gTMP_bch_scl, "%s/%s.scl", gTMP_bch_dir, gbch);

    printf("TMP BCH AUX %s NETS %s NODES %s PL %s WTS %s SCL %s\n", sing_fn_aux,
           sing_fn_nets, sing_fn_nodes, sing_fn_pl, sing_fn_wts, sing_fn_scl);


    PIN *pin = NULL;
    NET *net = NULL;
    MODULE *curModule = NULL;
    TERM *curTerminal = NULL;
    ROW *rwp = NULL;


    // write aux
    FILE *fp_aux = fopen(gTMP_bch_aux, "w");
    fprintf(fp_aux, "RowBasedPlacement :  %s  %s  %s  %s  %s\n", sing_fn_nodes,
            sing_fn_nets, sing_fn_wts, sing_fn_pl, sing_fn_scl);

    fclose(fp_aux);

    // write wts 
    FILE *fp_wts = fopen(gTMP_bch_wts, "w");
    fputs("UCLA wts 1.0\n", fp_wts);
    fputs("# Created\t:\tDec 27 2004\n", fp_wts);
    fputs(
        "# User   \t:\tGi-Joon Nam & Mehmet Yildiz at IBM Austin "
        "Research({gnam, mcan}@us.ibm.com)\n",
        fp_wts);

    fclose(fp_wts);

    // write nets
    FILE *fp_nets = fopen(gTMP_bch_nets, "w");
    fputs("UCLA nets 1.0 \n", fp_nets);
    fputs("# Created	:	Jan  6 2005\n", fp_nets);
    fputs(
        "# User   	:	Gi-Joon Nam & Mehmet Yildiz at IBM Austin "
        "Research({gnam, mcan}@us.ibm.com)\n",
        fp_nets);

    fputs("\n", fp_nets);

    fprintf(fp_nets, "NumNets : %d\n", netCNT);
    fprintf(fp_nets, "NumPins : %d\n", pinCNT);

    fputs("\n", fp_nets);

    char io[BUF_SZ];
    for(i = 0; i < netCNT; i++) {
        net = &netInstance[i];

        fprintf(fp_nets, "NetDegree : %d   %s\n", net->pinCNTinObject2,
                net->name);
        for(j = 0; j < net->pinCNTinObject2; j++) {
            pin = net->pin[j];

            moduleID = pin->moduleID;
            pinIDinModule = pin->pinIDinModule;

            if(pin->IO == 0)
                strcpy(io, "I");
            else if(pin->IO == 1)
                strcpy(io, "O");
            else  // B
                strcpy(io, "B");

            if(pin->term == 0) {
                curModule = &moduleInstance[moduleID];

                pinIDinModule = pin->pinIDinModule;

                fprintf(fp_nets, "\t%s\t%s : %.6lf\t%.6lf\n", curModule->name,
                        io, curModule->pof[pinIDinModule].x,
                        curModule->pof[pinIDinModule].y);
            }
            else {
                moduleID = pin->moduleID;
                curTerminal = &terminalInstance[moduleID];

                fprintf(fp_nets, "\t%s\t%s : %.6lf\t%.6lf\n", curTerminal->name,
                        io, curTerminal->pof[pinIDinModule].x,
                        curTerminal->pof[pinIDinModule].y);
            }
        }
    }

    fclose(fp_nets);

    // write nodes
    FILE *fp_nodes = fopen(gTMP_bch_nodes, "w");
    fputs("UCLA nodes 1.0\n", fp_nodes);
    fputs("# Created\t:\tJan  6 2005\n", fp_nodes);
    fputs(
        "# User   \t:\tGi-Joon Nam & Mehmet Yildiz at IBM Austin "
        "Research({gnam, mcan}@us.ibm.com)\n",
        fp_nodes);

    fputs("\n", fp_nodes);

    fprintf(fp_nodes, "NumNodes :  \t%d\n", moduleCNT + terminalCNT);
    fprintf(fp_nodes, "NumTerminals :  \t%d\n", terminalCNT + gmov_mac_cnt);

    for(i = 0; i < moduleCNT; i++) {
        curModule = &moduleInstance[i];

        if(curModule->flg == StdCell) {
            fprintf(fp_nodes, "\t%s\t%d\t%d\n", curModule->name,
                    prec2int(curModule->size.x), prec2int(curModule->size.y));
        }
        else {
            fprintf(fp_nodes, "\t%s\t%d\t%d\tterminal\n", curModule->name,
                    prec2int(curModule->size.x), prec2int(curModule->size.y));
        }
    }

    for(i = 0; i < terminalCNT; i++) {
        curTerminal = &terminalInstance[i];

        fprintf(fp_nodes, "\t%s\t%d\t%d\tterminal\n", curTerminal->name,
                prec2int(curTerminal->size.x), prec2int(curTerminal->size.y));
    }

    fclose(fp_nodes);

    // write pl files
    FILE *fp_pl = fopen(gTMP_bch_pl, "w");
    fputs("UCLA pl 1.0\n", fp_pl);
    fputs("# Created\t:\tJan  6 2005\n", fp_pl);
    fputs(
        "# User   \t:\tGi-Joon Nam & Mehmet Yildiz at IBM Austin "
        "Research({gnam, mcan}@us.ibm.com)\n",
        fp_pl);

    fputs("\n", fp_pl);

    for(i = 0; i < moduleCNT; i++) {
        curModule = &moduleInstance[i];

        if(curModule->flg == StdCell) {
            fprintf(fp_pl, "%s\t%.6lf\t%.6lf\t: N\n", curModule->name,
                    curModule->pmin.x, curModule->pmin.y);
        }
        else {
            fprintf(fp_pl, "%s\t%d\t%d\t: N /FIXED\n", curModule->name,
                    curModule->pmin_lg.x, curModule->pmin_lg.y);
        }
    }

    for(i = 0; i < terminalCNT; i++) {
        curTerminal = &terminalInstance[i];

        fprintf(fp_pl, "%s\t%d\t%d\t: N /FIXED\n", curTerminal->name,
                prec2int(curTerminal->pmin.x), prec2int(curTerminal->pmin.y));
    }

    fclose(fp_pl);

    // write scl files
    FILE *fp_scl = fopen(gTMP_bch_scl, "w");
    fputs("UCLA scl 1.0\n", fp_scl);
    fputs("# Created\t:\tJan  6 2005\n", fp_scl);
    fputs(
        "# User   \t:\tGi-Joon Nam & Mehmet Yildiz at IBM Austin "
        "Research({gnam, mcan}@us.ibm.com)\n",
        fp_scl);

    fputs("\n", fp_scl);
    fprintf(fp_scl, "NumRows :  \t%d\n", row_cnt);
    fputs("\n", fp_scl);

    for(i = 0; i < row_cnt; i++) {
        rwp = &row_st[i];

        fprintf(fp_scl, "CoreRow Horizontal\n");
        fprintf(fp_scl, "  Coordinate    :   %d\n", rwp->pmin.y);
        fprintf(fp_scl, "  Height        :   12\n");
        fprintf(fp_scl, "  Sitewidth     :    1\n");
        fprintf(fp_scl, "  Sitespacing   :    1\n");
        fprintf(fp_scl, "  Siteorient    :    1\n");
        fprintf(fp_scl, "  Sitesymmetry  :    1\n");
        fprintf(fp_scl, "  SubrowOrigin  :    %d\tNumSites  :  %d\n",
                rwp->pmin.x, rwp->x_cnt);
        fprintf(fp_scl, "End\n");
    }

    fclose(fp_scl);

    /////////////////////////////////
}

void runtimeError(string error_text) {
    fprintf(stderr, "ERROR: %s \n", error_text.c_str());
    fprintf(stderr, "Aborting !! \n");
    fflush(stdout);
    fflush(stderr);
    exit(1);
}

void read_routing_file(char *dir) {
    char route_file[BUFFERSIZE];
    char temp[BUFFERSIZE];
    char cmd[BUFFERSIZE];
    char netName[BUFFERSIZE];
    int netIdx;
    int fromX, fromY, fromL, toX, toY, toL, layer;
    char line[LINESIZE];
    bool flag = false;
    struct NET *net_tmp = NULL;
    struct FPOS from, to;

    sprintf(route_file, "%s/%s.est", dir, gbch);
    FILE *fp = fopen(route_file, "r");

    if((fp = fopen(route_file, "r")) == NULL) {
        sprintf(error_text, "ERROR: Cannot open: %s file", route_file);
        runtimeError(error_text);
    }

    while(!feof(fp)) {
        *line = '\0';
        fgets(line, LINESIZE, fp);
        sscanf(line, "%s%*s", temp);

        if(strlen(line) < 3)
            continue;

        if(temp[0] == 'n') {
            sscanf(line, "%s %d%*s", netName, &netIdx);
            flag = true;
            continue;
        }

        if(temp[0] == '!') {
            flag = false;
            continue;
        }

        if(flag) {
            if(temp[0] == '(') {
                sscanf(line, "(%d,%d,%d)-(%d,%d,%d)%*s", &fromX, &fromY, &fromL,
                       &toX, &toY, &toL);
                if(fromL != toL)
                    continue;
                from.x = fromX;
                from.y = fromY;
                from.z = fromL;
                to.x = toX;
                to.y = toY;
                to.z = toL;
                layer = fromL;
                net_tmp = &netInstance[netIdx];
                net_tmp->routing_tracks.push_back(
                    ROUTRACK(from, to, layer, netIdx));
            }
        }
    }
    sprintf(cmd, "rm -rf %s", route_file);
    system(cmd);
}

void delete_input_files_in(char *dir) {
    char cmd[BUF_SZ];
    sprintf(cmd, "rm -rf %s/%s.aux", dir, gbch);
    system(cmd);
    sprintf(cmd, "rm -rf %s/%s.nodes", dir, gbch);
    system(cmd);
    sprintf(cmd, "rm -rf %s/%s.nets", dir, gbch);
    system(cmd);
    sprintf(cmd, "rm -rf %s/%s.route", dir, gbch);
    system(cmd);
    sprintf(cmd, "rm -rf %s/%s.shapes", dir, gbch);
    system(cmd);
    sprintf(cmd, "rm -rf %s/%s.scl", dir, gbch);
    system(cmd);
    sprintf(cmd, "rm -rf %s/%s.wts", dir, gbch);
    system(cmd);
    sprintf(cmd, "rm -rf %s/%s_orig.pl", dir, gbch);
    system(cmd);
    sprintf(cmd, "rm -rf %s/%s.pl", dir, gbch);
    system(cmd);
}
void copy_original_SB_files_to_Dir(char *dir) {
    char cmd[BUF_SZ];
    char SBorigDirFile[BUF_SZ];

    sprintf(SBorigDirFile, "%s/SB_orig/%s/%s", currentDir, gbch, gbch);

    sprintf(cmd, "mkdir -p %s", dir);
    system(cmd);

    sprintf(cmd, "cp -rf %s.aux %s/.", SBorigDirFile, dir);
    system(cmd);
    sprintf(cmd, "cp -rf %s.nodes %s/.", SBorigDirFile, dir);
    system(cmd);
    sprintf(cmd, "cp -rf %s.nets %s/.", SBorigDirFile, dir);
    system(cmd);
    sprintf(cmd, "cp -rf %s.pl %s/%s_orig.pl", SBorigDirFile, dir, gbch);
    system(cmd);
    sprintf(cmd, "cp -rf %s.scl %s/.", SBorigDirFile, dir);
    system(cmd);
    sprintf(cmd, "cp -rf %s.wts %s/.", SBorigDirFile, dir);
    system(cmd);
    sprintf(cmd, "cp -rf %s.shapes %s/.", SBorigDirFile, dir);
    system(cmd);
    sprintf(cmd, "cp -rf %s.route %s/.", SBorigDirFile, dir);
    system(cmd);
}

void link_original_SB_files_to_Dir(char *dir) {
    char cmd[BUF_SZ];
    char SBorigDirFile[BUF_SZ];

    sprintf(SBorigDirFile, "%s/SB_orig/%s/%s", currentDir, gbch, gbch);

    sprintf(cmd, "mkdir -p %s", dir);
    system(cmd);

    sprintf(cmd, "ln -snf %s.aux %s/.", SBorigDirFile, dir);
    system(cmd);
    sprintf(cmd, "ln -snf %s.nodes %s/.", SBorigDirFile, dir);
    system(cmd);
    sprintf(cmd, "ln -snf %s.nets %s/.", SBorigDirFile, dir);
    system(cmd);
    sprintf(cmd, "cp -rf %s.pl %s/%s_orig.pl", SBorigDirFile, dir, gbch);
    system(cmd);
    sprintf(cmd, "ln -snf %s.scl %s/.", SBorigDirFile, dir);
    system(cmd);
    sprintf(cmd, "ln -snf %s.wts %s/.", SBorigDirFile, dir);
    system(cmd);
    sprintf(cmd, "ln -snf %s.shapes %s/.", SBorigDirFile, dir);
    system(cmd);
    sprintf(cmd, "ln -snf %s.route %s/.", SBorigDirFile, dir);
    system(cmd);
}


///////////////////////////////////////////////////
//  
//  bookshelf writing function.
//  (aux, nodes, shapes, net, pl, scl, wts)
//
///////////////////////////////////////////////////


///////////////////////////////////////////////////
// *.aux writing
void WriteAux( char* dir_tier, bool isShapeDrawing ) {
    char fn_aux[BUF_SZ] = {0, };
    sprintf(fn_aux, "%s/%s.aux", dir_tier, gbch);
    FILE* fp_aux = fopen(fn_aux, "w");

    if( isShapeDrawing ) {
        fprintf(fp_aux,
                "RowBasedPlacement : %s.nodes %s.nets %s.wts %s.pl %s.scl %s.shapes\n", gbch,
                gbch, gbch, gbch, gbch, gbch);
    } 
    else {
        fprintf(fp_aux,
                "RowBasedPlacement : %s.nodes %s.nets %s.wts %s.pl %s.scl\n", gbch,
                gbch, gbch, gbch, gbch);
    }

    fclose(fp_aux);
}

// *.shape writing - mgwoo
void WriteShapes( char* dir_tier, bool isShapeDrawing) {
    if( !isShapeDrawing ) {
        return;
    }

    char fn_shapes[BUF_SZ] = {0, };
    sprintf(fn_shapes, "%s/%s.shapes", dir_tier, gbch);
    FILE* fp_shapes = fopen(fn_shapes, "w");
    
    // file header..
    fputs("shapes 1.0 \n", fp_shapes);
    fputs("# Created	:	Jan  6 2005\n", fp_shapes);
    fputs(
        "# User   	:	Gi-Joon Nam & Mehmet Yildiz at IBM Austin "
        "Research({gnam, mcan}@us.ibm.com)\n",
        fp_shapes);
    fputs("\n", fp_shapes);


    fprintf(fp_shapes, "NumNonRectangularNodes  :  %d\n\n", numNonRectangularNodes);
    
    for(auto& curShapeNode : shapeMap) {
        fprintf(fp_shapes, "%s  :  %ld\n", curShapeNode.first.c_str(), 
                                           curShapeNode.second.size()); 
        for(auto& curIdx : curShapeNode.second) {
            fprintf(fp_shapes, "\t%s\t%d\t%d\t%d\t%d\n", 
            shapeStor[curIdx].name.c_str(), 
            (int)shapeStor[curIdx].llx, (int)shapeStor[curIdx].lly, 
            (int)shapeStor[curIdx].width, (int)shapeStor[curIdx].height ); 
        }
    }

    fclose(fp_shapes);
}

// *.weignt writing - nothing to write..
void WriteWts( char* dir_tier ) {
    
    char fn_wts[BUF_SZ] = {0, };
    sprintf(fn_wts, "%s/%s.wts", dir_tier, gbch);
    FILE* fp_wts = fopen(fn_wts, "w");
    
    fputs("UCLA wts 1.0\n", fp_wts);
    fputs("# Created\t:\tDec 27 2004\n", fp_wts);
    fputs(
        "# User   \t:\tGi-Joon Nam & Mehmet Yildiz at IBM Austin "
        "Research({gnam, mcan}@us.ibm.com)\n",
        fp_wts);

    fclose(fp_wts);
}

// *.nodes writing
void WriteNodes( char* dir_tier, int curLayer, int lab, int pin_term_cnt, bool isShapeDrawing ) {
    char fn_nodes[BUF_SZ] = {0, };
    sprintf(fn_nodes, "%s/%s.nodes", dir_tier, gbch);
    FILE* fp_nodes= fopen(fn_nodes, "w");

    fputs("UCLA pl 2.0 \n", fp_nodes);
    fputs("# Created	:	Jan  6 2005\n", fp_nodes);
    fputs(
        "# User   	:	Gi-Joon Nam & Mehmet Yildiz at IBM Austin "
        "Research({gnam, mcan}@us.ibm.com)\n",
        fp_nodes);
    fputs("\n", fp_nodes);

    TIER *tier = &tier_st[curLayer];
    tier->row_term_cnt = 0;
    
    int term_cnt = (lab)? tier->mac_cnt + terminalCNT + pin_term_cnt :
                          terminalCNT + pin_term_cnt;

    term_cnt = (!isShapeDrawing)? term_cnt + totalShapeCount - numNonRectangularNodes :
                                  term_cnt;

//    if(lab)
//        term_cnt = tier->mac_cnt + terminalCNT + pin_term_cnt;
//    else
//        term_cnt = terminalCNT + pin_term_cnt;

    int node_cnt = tier->modu_cnt + term_cnt;

    fprintf(fp_nodes, "NumNodes :  \t%d\n", node_cnt);
    fprintf(fp_nodes, "NumTerminals :  \t%d\n", term_cnt);

    for(int i = 0; i < tier->modu_cnt; i++) {
        MODULE* modu = tier->modu_st[i];
        fprintf(fp_nodes, "%s %d %d\n", modu->name, (int)(modu->size.x + 0.5),
                (int)(modu->size.y + 0.5));
    }

    for(int i = 0; i < terminalCNT; i++) {
        TERM* curTerminal = &terminalInstance[i];
        string termString = isShapeDrawing? ((curTerminal->isTerminalNI)? "terminal_NI" : "terminal") : "terminal";

        if( isShapeDrawing ) {
            fprintf(fp_nodes, "%s %d %d %s\n", curTerminal->name,
                    (int)(curTerminal->size.x + 0.5),
                    (int)(curTerminal->size.y + 0.5),
                    termString.c_str());
        }
        else {
            // considering shapeMaps - mgwoo 
            auto shapeMapIter = shapeMap.find(curTerminal->name);
            if(shapeMapIter == shapeMap.end()) {
                int sizeX = (curTerminal->isTerminalNI)? 0 : INT_CONVERT(curTerminal->size.x);
                int sizeY = (curTerminal->isTerminalNI)? 0 : INT_CONVERT(curTerminal->size.y);

                fprintf(fp_nodes, "%s %d %d %s\n", curTerminal->name,
                        sizeX, sizeY,
                        termString.c_str());
            }
            else {
                for(auto& curShapeIdx : shapeMap[curTerminal->name]) {
                    int sizeX = (curTerminal->isTerminalNI)? 0 : INT_CONVERT(shapeStor[curShapeIdx].width);
                    int sizeY = (curTerminal->isTerminalNI)? 0 : INT_CONVERT(shapeStor[curShapeIdx].height);

                    fprintf( fp_nodes, "%s %d %d %s\n", 
                             string( string(curTerminal->name) 
                                   + string("/") 
                                   + string(shapeStor[curShapeIdx].name) ).c_str(),
                             sizeX, sizeY,
                             termString.c_str() );

                }
            }
        }
    }

    if(lab)  // MMS-3D place
    {
        for(int i = 0; i < tier->mac_cnt; i++) {
            MODULE* mac = tier->mac_st[i];
            fprintf(fp_nodes, "%s %d %d terminal\n", mac->name,
                    (int)(mac->size.x + 0.5), (int)(mac->size.y + 0.5));
        }
    }

    for(int i = 0; i < pinCNT; i++) {
        PIN* pin = &pinInstance[i];
        if(pin->tier != curLayer && !pin->term) {
            fprintf(fp_nodes, "fakePin%d %d %d terminal\n", pin->gid, 0, 0);
        }
    }

    fclose(fp_nodes);
}

// *.nets writing
void WriteNet( char* dir_tier, int curLayer, int pin_cnt, int net_cnt, bool isShapeDrawing) {
    
    char fn_nets[BUF_SZ] = {0, };
    sprintf(fn_nets, "%s/%s.nets", dir_tier, gbch);
    FILE *fp_nets = fopen(fn_nets, "w");

    fputs("UCLA pl 2.0 \n", fp_nets);
    fputs("# Created	:	Jan  6 2005\n", fp_nets);
    fputs(
        "# User   	:	Gi-Joon Nam & Mehmet Yildiz at IBM Austin "
        "Research({gnam, mcan}@us.ibm.com)\n",
        fp_nets);
    fputs("\n", fp_nets);

    fprintf(fp_nets, "NumNets : %d\n", net_cnt);
    fprintf(fp_nets, "NumPins : %d\n", pin_cnt);

    fputs("\n", fp_nets);

    int total_pinCNTinObject2 = 0;

    for(int i = 0; i < netCNT; i++) {
        NET* net = &netInstance[i];

        total_pinCNTinObject2 += net->pinCNTinObject2;

        if(net->tier != curLayer)
            continue;

        fprintf(fp_nets, "NetDegree : %d   %s\n", net->pinCNTinObject_tier,
                net->name);

        for(int j = 0; j < net->pinCNTinObject2; j++) {
            PIN* pin = net->pin2[j];

            int moduleID = pin->moduleID;
            int pinIDinModule = pin->pinIDinModule;

            char io[BUF_SZ] = {0, };
            if(pin->IO == 0)
                strcpy(io, "I");
            else if(pin->IO == 1)
                strcpy(io, "O");
            else
                strcpy(io, "B");

            if(pin->term) {
                TERM* curTerminal = &terminalInstance[moduleID];
                if( isShapeDrawing ) {
                    fprintf(fp_nets, "\t%s\t%s : %.6lf\t%.6lf\n", curTerminal->name,
                            io, 
                            curTerminal->pof[pinIDinModule].x,
                            curTerminal->pof[pinIDinModule].y);
                }
                else {
                    if( shapeMap.find(curTerminal->name) == shapeMap.end() ) {
                        fprintf(fp_nets, "\t%s\t%s : %.6lf\t%.6lf\n", curTerminal->name,
                                io, 
                                curTerminal->pof[pinIDinModule].x,
                                curTerminal->pof[pinIDinModule].y);
                    }
                    else {
                        // convert into "o506100/Shape_0"
                        SHAPE* curShape = &shapeStor[ shapeMap[ curTerminal-> name ][0] ];
                        string concatedName = string( curTerminal->name ) +
                                              string("/") +
                                              string( curShape->name );

                        prec shapeCenterX = curShape->llx + curShape->width/2;
                        prec shapeCenterY = curShape->lly + curShape->height/2;
                                               
                        // covert as "o506100/shape_0"'s information
                        fprintf(fp_nets, "\t%s\t%s : %.6lf\t%.6lf\n", 
                                concatedName.c_str(),
                                io, 
                                curTerminal->pof[pinIDinModule].x 
                                    + (curTerminal->center.x - shapeCenterX),
                                curTerminal->pof[pinIDinModule].y
                                    + (curTerminal->center.y - shapeCenterY) );
                    }
                }
            }
            else if(pin->tier != curLayer) {
                fprintf(fp_nets, "\tfakePin%d\t%s : %.6lf\t%.6lf\n", pin->gid,
                        io, 0.0, 0.0);
            }
            else if(pin->tier == curLayer) {
                if(pin->term == 0) {
                    MODULE* modu = &moduleInstance[moduleID];

                    fprintf(fp_nets, "\t%s\t%s : %.6lf\t%.6lf\n", modu->name,
                            io, modu->pof[pinIDinModule].x,
                            modu->pof[pinIDinModule].y);
                }
                else {
                    TERM* curTerminal = &terminalInstance[moduleID];

                    fprintf(fp_nets, "\t%s\t%s : %.6lf\t%.6lf\n",
                            curTerminal->name, io,
                            curTerminal->pof[pinIDinModule].x,
                            curTerminal->pof[pinIDinModule].y);
                }
            }
        }
    }

    fclose(fp_nets);
}

// *.pl writing
void WritePl( char* dir_tier, int curLayer, int lab, bool isShapeDrawing ) {
    
    char fn_pl[BUF_SZ] = {0, };
    sprintf(fn_pl, "%s/%s.pl", dir_tier, gbch);
    FILE* fp_pl = fopen(fn_pl, "w");

    fputs("UCLA pl 1.0\n", fp_pl);
    fputs("# Created\t:\tJan  6 2005\n", fp_pl);
    fputs(
        "# User   \t:\tGi-Joon Nam & Mehmet Yildiz at IBM Austin "
        "Research({gnam, mcan}@us.ibm.com)\n",
        fp_pl);

    fputs("\n", fp_pl);
    
    TIER *tier = &tier_st[curLayer];

    for(int i = 0; i < tier->modu_cnt; i++) {
        MODULE* modu = tier->modu_st[i];
        fprintf(fp_pl, "%s\t%.6lf\t%.6lf : N\n", modu->name, modu->pmin.x,
                modu->pmin.y);
    }

    for(int i = 0; i < terminalCNT; i++) {
        TERM* curTerminal = &terminalInstance[i];
        string fixedStr = (curTerminal->isTerminalNI)? "FIXED_NI" : "FIXED";

        if( isShapeDrawing ) {
            fprintf(fp_pl, "%s\t%d\t%d\t: N /%s\n", curTerminal->name,
                    prec2int(curTerminal->pmin.x),
                    prec2int(curTerminal->pmin.y),
                    fixedStr.c_str());
        }
        else {
            if ( shapeMap.find( curTerminal->name ) == shapeMap.end() ) {
                fprintf(fp_pl, "%s\t%d\t%d\t: N /%s\n", curTerminal->name,
                        prec2int(curTerminal->pmin.x),
                        prec2int(curTerminal->pmin.y),
                        fixedStr.c_str());
            }
            else {
                for(auto& curShapeIdx : shapeMap[curTerminal->name]) {
                    fprintf(fp_pl, "%s\t%d\t%d\t: N /%s\n", 
                            string( string(curTerminal->name)
                                    + string("/")
                                    + string(shapeStor[curShapeIdx].name) ).c_str(),
                            prec2int(shapeStor[curShapeIdx].llx),
                            prec2int(shapeStor[curShapeIdx].lly),
                                fixedStr.c_str());
                }
            }
        }
    }

    if(lab)  // MMS-3D place
    {
        for(int i = 0; i < tier->mac_cnt; i++) {
            MODULE* mac = tier->mac_st[i];
            fprintf(fp_pl, "%s\t%d\t%d\t: N /FIXED\n", mac->name,
                    mac->pmin_lg.x, mac->pmin_lg.y);
        }
    }

    for(int i = 0; i < pinCNT; i++) {
        PIN* pin = &pinInstance[i];
        if(!pin->term && pin->tier != curLayer)
            fprintf(fp_pl, "fakePin%d\t%.6lf\t%.6lf\t: N /FIXED\n", pin->gid,
                    pin->fp.x, pin->fp.y);
    }

    fclose(fp_pl);

}

// *.scl writing
void WriteScl( char* dir_tier, int curLayer ) {
    char fn_scl[BUF_SZ];
    sprintf(fn_scl, "%s/%s.scl", dir_tier, gbch);
    FILE* fp_scl = fopen(fn_scl, "w");

    fputs("UCLA scl 1.0\n", fp_scl);
    fputs("# Created\t:\tJan  6 2005\n", fp_scl);
    fputs(
        "# User   \t:\tGi-Joon Nam & Mehmet Yildiz at IBM Austin "
        "Research({gnam, mcan}@us.ibm.com)\n",
        fp_scl);

    TIER *tier = &tier_st[curLayer];

    fputs("\n", fp_scl);
    fprintf(fp_scl, "NumRows :  \t%d\n", tier->row_cnt + tier->row_term_cnt);
    fputs("\n", fp_scl);

    for(int i = 0; i < tier->row_cnt; i++) {
        ROW* row = &(tier->row_st[i]);

        fprintf(fp_scl, "CoreRow Horizontal\n");
        fprintf(fp_scl, "  Coordinate    :   %d\n", row->pmin.y);
        fprintf(fp_scl, "  Height        :   %d\n", (int)(rowHeight+0.5f));
        fprintf(fp_scl, "  Sitewidth     :    %d\n", row->site_wid);
        fprintf(fp_scl, "  Sitespacing   :    %d\n", row->site_spa);
        fprintf(fp_scl, "  Siteorient    :    %s\n", (row->isYSymmetry)? "Y" : "1");
        fprintf(fp_scl, "  Sitesymmetry  :    %s\n", row->ori.c_str());
        fprintf(fp_scl, "  SubrowOrigin  :    %d\tNumSites  :  %d\n",
                row->pmin.x, row->x_cnt);
        fprintf(fp_scl, "End\n");
    }

    fclose(fp_scl);
}

// write bookshelf's output
// z : current tier's index
void WriteBookshelfWithTier(int z, int lab, bool isShapeDrawing) {

    PIN *pin = NULL;
    NET *net = NULL;

    // create directory for tier's result
    char dir_tier[BUF_SZ] = {0, };
    sprintf(dir_tier, "%s/tiers/%d", dir_bnd, z);

    char cmd[BUF_SZ] = {0, };
    sprintf(cmd, "mkdir -p %s", dir_tier);
    system(cmd);
    

    /////////////////////////////////////////////////////////////////////////////
    // count for pin_term_cnt, net_cnt 
    // (Required for *.nodes)
    int pin_term_cnt = 0;
    int net_cnt = 0;

    // to calculate net_cnt
    for(int i=0; i<netCNT; i++) {
        netInstance[i].tier = -1;
    }

    for(int i = 0; i < pinCNT; i++) {
        pin = &pinInstance[i];
        if(pin->term) {
        }
        else if(pin->tier != z) {
            pin_term_cnt++;
        }
        else if(pin->tier == z) {
            net = &netInstance[pin->netID];
            if(net->tier < z) {
                net->tier = z;
                net_cnt++;
            }
        }
    }
    
    /////////////////////////////////////////////////////////////////////////////
    // count for pin_cnt
    // (Required for *.net)
    
    int pin_cnt = 0;
    for(int i = 0; i < netCNT; i++) {
        net = &netInstance[i];
        net->pinCNTinObject_tier = 0;

        if(net->tier != z)
            continue;

        for(int j = 0; j < net->pinCNTinObject2; j++) {
            pin = net->pin[j];

            /* if(pin->tier <= z) */
            /*   { */
            net->pinCNTinObject_tier++;
            pin_cnt++;
            /* } */
        }
    }

    WriteAux( dir_tier, isShapeDrawing );
//    cout << "aux finish" << endl;
    WriteShapes( dir_tier, isShapeDrawing );
//    cout << "shape finish" << endl;
    WriteNodes( dir_tier, z, lab, pin_term_cnt, isShapeDrawing );
//    cout << "nodes finish" << endl;
    WriteNet( dir_tier, z, pin_cnt, net_cnt, isShapeDrawing );
//    cout << "net finish" << endl;
    WritePl( dir_tier, z, lab, isShapeDrawing );
//    cout << "pl finish" << endl;
    
    // shapeSupport doesn't affect to below function
    WriteScl( dir_tier, z ); 
    WriteWts( dir_tier ); 
}


// 
//
// read *.pl back to the replace to do some post-processing -- nodesMap
//
// See also ReadPlLefDef(const char* fileName) in lefdefIO.cpp 
//
void ReadPlBookshelf(const char *fileName) {
    FPOS pof;
    FILE *fp = fopen(fileName, "r");
    if( !fp ) {
        runtimeError( fileName + string( " is not Exist!" ) );
    }

    char buf[BUF_SZ], name[BUF_SZ];
    int moduleID = 0;
    while(fgets(buf, BUF_SZ, fp)) {
        char* token = strtok(buf, " \t\n");
        if(!token || token[0] == '#' || !strcmp(token, "UCLA"))
            continue;

        strcpy(name, token);

        if(name[0] == 'f' && name[1] == 'a' && name[2] == 'k' && name[3] == 'e')
            continue;

        // if(!getCellIndex ( name , & moduleID , & isTerminal,
        // &isTerminalNI))
        // continue;

        bool isTerminal = false;
        auto nodesMapPtr = nodesMap.find(string(name));
        if(nodesMapPtr == nodesMap.end()) {
            continue;
        }
        else {
            isTerminal = nodesMapPtr->second.isTerminal;
            moduleID = nodesMapPtr->second.index;
        }

        if(!isTerminal) {
            MODULE* curModule = &moduleInstance[moduleID];

            if(curModule->flg == Macro)
                continue;

            token = strtok(NULL, " \t\n");
            curModule->pmin.x = atof(token);

            token = strtok(NULL, " \t\n");
            curModule->pmin.y = atof(token);

            curModule->pmax.x = curModule->pmin.x + curModule->size.x;
            curModule->pmax.y = curModule->pmin.y + curModule->size.y;

            curModule->center.x = 0.5 * (curModule->pmin.x + curModule->pmax.x);
            curModule->center.y = 0.5 * (curModule->pmin.y + curModule->pmax.y);

            for(int i = 0; i < curModule->pinCNTinObject; i++) {
                PIN* pin = curModule->pin[i];
                pof = curModule->pof[i];
                pin->fp.x = curModule->center.x + pof.x;
                pin->fp.y = curModule->center.y + pof.y;
            }
        }
        else {
            continue;
        }
    }

    fclose(fp);
}

void output_tier_pl_global_router(char *dir, int z, int lab) {
    char fn_pl[BUF_SZ];
    // char            dir_router[BUF_SZ];
    // char            cmd[BUF_SZ];

    FILE *fp_pl = NULL;

    struct TIER *tier = &tier_st[z];
    struct MODULE *modu = NULL;
    struct MODULE *mac = NULL;
    struct TERM *curTerminal = NULL;
    struct PIN *pin = NULL;

    // sprintf (dir_router, "%s/router/tier%d/inflarion_iter%d", dir_bnd, z,
    // infl_iter);

    // sprintf (cmd, "mkdir -p %s", dir_router);
    // system (cmd);

    cout << "INFO:  "
         << "Temp. Working Dir= " << endl
         << "       " << dir << endl;
    sprintf(fn_pl, "%s/%s.pl", dir, gbch);

    fp_pl = fopen(fn_pl, "w");

    fputs("UCLA pl 1.0\n", fp_pl);
    fputs("# Created\t:\tJan  6 2005\n", fp_pl);
    fputs(
        "# User   \t:\tGi-Joon Nam & Mehmet Yildiz at IBM Austin "
        "Research({gnam, mcan}@us.ibm.com)\n",
        fp_pl);
    fputs("\n", fp_pl);

    for(int i = 0; i < tier->modu_cnt; i++) {
        modu = tier->modu_st[i];
        fprintf(fp_pl, "%s\t%.6lf\t%.6lf : N\n", modu->name, modu->pmin.x,
                modu->pmin.y);
    }

    for(int i = 0; i < terminalCNT; i++) {
        curTerminal = &terminalInstance[i];
        if(!curTerminal->isTerminalNI) {
            fprintf(fp_pl, "%s\t%d\t%d\t: N /FIXED\n", curTerminal->name,
                    prec2int(curTerminal->pmin.x),
                    prec2int(curTerminal->pmin.y));
        }
        else {
            fprintf(fp_pl, "%s\t%d\t%d\t: N /FIXED_NI\n", curTerminal->name,
                    prec2int(curTerminal->pmin.x),
                    prec2int(curTerminal->pmin.y));
        }
    }

    if(lab) {  // MMS-3D place
        for(int i = 0; i < tier->mac_cnt; i++) {
            mac = tier->mac_st[i];
            fprintf(fp_pl, "%s\t%d\t%d\t: N /FIXED\n", mac->name,
                    mac->pmin_lg.x, mac->pmin_lg.y);
        }
    }

    for(int i = 0; i < pinCNT; i++) {
        pin = &pinInstance[i];
        if(!pin->term && pin->tier != z)
            fprintf(fp_pl, "fakePin%d\t%.6lf\t%.6lf\t: N /FIXED\n", pin->gid,
                    pin->fp.x, pin->fp.y);
    }
    fclose(fp_pl);
    cout << "INFO:  Temp. placement solution for Global Router has been "
            "written at "
         << endl
         << "       " << fn_pl << endl
         << endl;
}

void get_3d_dimension(struct POS *tier_min, struct POS *tier_max,
                      int *tier_row_cnt) {
    int xlen = 0;
    int ylen = 0;

    xlen = grow_pmax.x - grow_pmin.x;
    ylen = grow_pmax.y - grow_pmin.y;
    prec aspect_ratio = (prec)ylen / (prec)xlen;

    printf("INFO:  ASPECT RATIO=%.6lf\n", aspect_ratio);
    printf("INFO:  ROW = [(MinX=%f, MinY=%f), (MaxX=%f, MaxY=%f)]\n",
           grow_pmin.x, grow_pmin.y, grow_pmax.x, grow_pmax.y);
    printf("INFO:  TERM: (%lf, %lf) - (%lf, %lf)\n", terminal_pmin.x,
           terminal_pmin.y, terminal_pmax.x, terminal_pmax.y);

    shrunk_ratio = zeroFPoint;
    struct FPOS shrunk_len = zeroFPoint;

    // shrunk_len.x = pow ((1.0 + ExtraWSfor3D)*xlen*ylen
    //// / target_cell_den
    // / (aspect_ratio), 0.5);
    // shrunk_len.y = pow ((1.0 + ExtraWSfor3D)*xlen*ylen
    //// / target_cell_den
    // * aspect_ratio, 0.5);

    if(strcmp(bmFlagCMD.c_str(), "etc") == 0) {
        shrunk_len.x = pow(xlen * ylen / aspect_ratio, 0.5);
        shrunk_len.y = pow(xlen * ylen * aspect_ratio, 0.5);
    }
    else {
        shrunk_len.x = pow((1.0 + ExtraWSfor3D) * xlen * ylen /
                               target_cell_den / (numLayer * aspect_ratio),
                           0.5);
        shrunk_len.y =
            pow((1.0 + ExtraWSfor3D) * xlen * ylen / target_cell_den *
                    aspect_ratio / ((prec)numLayer),
                0.5);
    }

    // if ((strcmp(gbch, "newblue3")==0
    // ||strcmp(gbch, "newblue2")==0
    // ||strcmp(gbch, "bigblue3")==0) && numLayer >= 4) {
    // printf("bigblue3 OR newblue2 OR newblue3: accomodate all the macros by
    // inserting %2.0f%% whitespace to each tier!\n",
    // ExtraWSfor3D
    // /target_cell_den
    // *100);

    // shrunk_len.x = pow ((1.0 + MaxExtraWSfor3D)*xlen*ylen
    // / target_cell_den
    // / ((prec)numLayer*aspect_ratio), 0.5);
    // shrunk_len.y = pow ((1.0 + MaxExtraWSfor3D)*xlen*ylen
    // / target_cell_den
    // * aspect_ratio/((prec)numLayer), 0.5);

    if(shrunk_len.x < max_mac_dim.x || shrunk_len.y < max_mac_dim.y) {
        printf(
            "Cannot accomodate the largest macro by inserting %2.0f%% "
            "whitespace to each tier!\n",
            ExtraWSfor3D
                // /target_cell_den
                * 100);
        shrunk_len.x = pow((1.0 + MaxExtraWSfor3D) * xlen * ylen
                               // / target_cell_den
                               / ((prec)numLayer * aspect_ratio),
                           0.5);
        shrunk_len.y = pow((1.0 + MaxExtraWSfor3D) * xlen * ylen
                               // / target_cell_den
                               * aspect_ratio / ((prec)numLayer),
                           0.5);

        if(shrunk_len.x < max_mac_dim.x || shrunk_len.y < max_mac_dim.y) {
            printf(
                "INFO:  Cannot accomodate the largest macro by inserting "
                "%2.0f%% whitespace to each tier!\n",
                MaxExtraWSfor3D
                    // /target_cell_den
                    * 100);
            printf(" === Cannot conduct 3D-IC placement!\n");
            printf(" === Placement exits!\n");
            exit(0);
        }
        else {
            printf("INFO:  NOTE:  %2.0f%% Whitespace is Added to Each Tier.\n",
                   MaxExtraWSfor3D
                       // /target_cell_den
                       * 100);
            printf("INFO:    LARGEST MACRO SIZE (X,Y) = (%d, %d)\n",
                   max_mac_dim.x, max_mac_dim.y);
            printf("INFO:    TIER DIMENSION (X,Y) = (%.0lf, %.0lf)\n",
                   shrunk_len.x, shrunk_len.y);
        }
    }
    else {
        printf("INFO:  NOTE:  %2.0f%% Whitespace is Added to Each Tier.\n",
               ExtraWSfor3D
                   // /target_cell_den
                   * 100);
        printf("INFO:    LARGEST MACRO SIZE (X,Y) = (%d, %d)\n", max_mac_dim.x,
               max_mac_dim.y);
        printf("INFO:    TIER DIMENSION (X,Y) = (%.0lf, %.0lf)\n", shrunk_len.x,
               shrunk_len.y);
    }

    shrunk_ratio.x = shrunk_len.x / xlen;
    shrunk_ratio.y = shrunk_len.y / ylen;

    *tier_row_cnt = (int)((prec)row_cnt * shrunk_ratio.y + 0.5);  //+ 1;  igkang

    struct POS shrunk_len_int = zeroPoint;

    shrunk_len_int.x = (int)(shrunk_len.x + 0.5);
    shrunk_len_int.y = (*tier_row_cnt) * rowHeight;

    shrunk_ratio.x = (prec)shrunk_len_int.x / xlen;
    shrunk_ratio.y = (prec)shrunk_len_int.y / ylen;

    tier_max->x = grow_pmin.x + shrunk_len_int.x;
    tier_max->y = grow_pmin.y + shrunk_len_int.y;

    tier_min->x = grow_pmin.x;
    tier_min->y = grow_pmin.y;

    printf("shrunk_(xlen, ylen) = (%d, %d)\n", shrunk_len_int.x,
           shrunk_len_int.y);
}

// 
// update tier_min, tier_max, and tier_row_cnt
//
void get_mms_3d_dim(POS *tier_min, POS *tier_max,
                    int *tier_row_cnt) {
    int xlen = 0;
    int ylen = 0;

    xlen = grow_pmax.x - grow_pmin.x;
    ylen = grow_pmax.y - grow_pmin.y;
    prec aspect_ratio = (prec)ylen / (prec)xlen;

    printf("INFO:  ASPECT RATIO=%.6lf\n", aspect_ratio);
    printf("INFO:  ROW = (MinX=%f, MinY=%f)\n", grow_pmin.x, grow_pmin.y);
    printf("INFO:  ROW = (MaxX=%f, MaxY=%f)\n", grow_pmax.x, grow_pmax.y);
    printf("INFO:  TERM: (%lf, %lf) - (%lf, %lf)\n", terminal_pmin.x,
           terminal_pmin.y, terminal_pmax.x, terminal_pmax.y);

    shrunk_ratio = zeroFPoint;
    struct FPOS shrunk_len = zeroFPoint;

    // shrunk_len.x = pow ((1.0 + ExtraWSfor3D)*xlen*ylen
    //// / target_cell_den
    // / (aspect_ratio), 0.5);
    // shrunk_len.y = pow ((1.0 + ExtraWSfor3D)*xlen*ylen
    //// / target_cell_den
    // * aspect_ratio, 0.5);
    if(numLayer == 1) {
        // lutong
        // if (strcmp(bmFlagCMD.c_str(), "etc") == 0) {
        // shrunk_len.x = pow (xlen*ylen/aspect_ratio, 0.5);
        // shrunk_len.y = pow (xlen*ylen*aspect_ratio, 0.5);
        //} else {
        shrunk_len.x = xlen;
        shrunk_len.y = ylen;
        //}
    }
    else {
        if(strcmp(bmFlagCMD.c_str(), "etc") == 0) {
            shrunk_len.x = pow(xlen * ylen / aspect_ratio, 0.5);
            shrunk_len.y = pow(xlen * ylen * aspect_ratio, 0.5);
        }
        else {
            shrunk_len.x = pow((1.0 + ExtraWSfor3D) * xlen * ylen /
                                   target_cell_den / (numLayer * aspect_ratio),
                               0.5);
            shrunk_len.y =
                pow((1.0 + ExtraWSfor3D) * xlen * ylen / target_cell_den *
                        aspect_ratio / ((prec)numLayer),
                    0.5);
        }
    }

    // if ((strcmp(gbch, "newblue3")==0
    // ||strcmp(gbch, "newblue2")==0
    // ||strcmp(gbch, "bigblue3")==0) && numLayer >= 4) {
    // printf("bigblue3 OR newblue2 OR newblue3: accomodate all the macros by
    // inserting %2.0f%% whitespace to each tier!\n",
    // ExtraWSfor3D
    // /target_cell_den
    // *100);

    // shrunk_len.x = pow ((1.0 + MaxExtraWSfor3D)*xlen*ylen
    // / target_cell_den
    // / ((prec)numLayer*aspect_ratio), 0.5);
    // shrunk_len.y = pow ((1.0 + MaxExtraWSfor3D)*xlen*ylen
    // / target_cell_den
    // * aspect_ratio/((prec)numLayer), 0.5);
    if(shrunk_len.x < max_mac_dim.x || shrunk_len.y < max_mac_dim.y) {
        printf(
            "Cannot accomodate the largest macro by inserting %2.0f%% "
            "whitespace to each tier!\n",
            ExtraWSfor3D
                // /target_cell_den
                * 100);
        shrunk_len.x = pow((1.0 + MaxExtraWSfor3D) * xlen * ylen
                               // / target_cell_den
                               / ((prec)numLayer * aspect_ratio),
                           0.5);
        shrunk_len.y = pow((1.0 + MaxExtraWSfor3D) * xlen * ylen
                               // / target_cell_den
                               * aspect_ratio / ((prec)numLayer),
                           0.5);

        if(shrunk_len.x < max_mac_dim.x || shrunk_len.y < max_mac_dim.y) {
            printf(
                "INFO:  Cannot accomodate the largest macro by inserting "
                "%2.0f%% whitespace to each tier!\n",
                MaxExtraWSfor3D
                    // /target_cell_den
                    * 100);
            printf(" === Cannot conduct 3D-IC placement!\n");
            printf(" === Placement exits!\n");
            exit(0);
        }
        else {
            printf("INFO:  NOTE:  %2.0f%% Whitespace is Added to Each Tier.\n",
                   MaxExtraWSfor3D
                       // /target_cell_den
                       * 100);
            printf("INFO:    LARGEST MACRO SIZE (X,Y) = (%d, %d)\n",
                   max_mac_dim.x, max_mac_dim.y);
            printf("INFO:    TIER DIMENSION (X,Y) = (%.0lf, %.0lf)\n",
                   shrunk_len.x, shrunk_len.y);
        }
    }
    else {
        printf("INFO:  NOTE:  %2.0f%% Whitespace is Added to Each Tier.\n",
               ExtraWSfor3D
                   // /target_cell_den
                   * 100);
        printf("INFO:    LARGEST MACRO SIZE (X,Y) = (%d, %d)\n", max_mac_dim.x,
               max_mac_dim.y);
        printf("INFO:    TIER DIMENSION (X,Y) = (%.0lf, %.0lf)\n", shrunk_len.x,
               shrunk_len.y);
    }

    shrunk_ratio.x = shrunk_len.x / xlen;
    shrunk_ratio.y = shrunk_len.y / ylen;

    //*tier_row_cnt = (int)((prec)row_cnt * shrunk_ratio.y + 0.5);    //+ 1;
    // igkang
    *tier_row_cnt = row_cnt;  // lutong

    struct POS shrunk_len_int = zeroPoint;

    // shrunk_len_int.x = (int)(shrunk_len.x + 0.5);
    shrunk_len_int.x = xlen;  // lutong
    shrunk_len_int.y = (*tier_row_cnt) * rowHeight;

    shrunk_ratio.x = (prec)shrunk_len_int.x / xlen;
    shrunk_ratio.y = (prec)shrunk_len_int.y / ylen;

    tier_max->x = grow_pmin.x + shrunk_len_int.x;
    tier_max->y = grow_pmin.y + shrunk_len_int.y;

    tier_min->x = grow_pmin.x;
    tier_min->y = grow_pmin.y;

    if(numLayer > 1) {
        printf("shrunk_(xlen, ylen) = (%d, %d)\n", shrunk_len_int.x,
               shrunk_len_int.y);
    }
}

void get_ibm_3d_dim(struct POS *tier_min, struct POS *tier_max,
                    int *tier_row_cnt) {
    if(INPUT_FLG != IBM) {
        printf("Error function call of get_3d_dim, benchmark suite not IBM.\n");
        exit(1);
    }

    tier_min->x = grow_pmin.x;
    tier_min->y = grow_pmin.y;

    // directly use 3D-IC core region from NTUplace3

    if(!strcmp(gbch, "ibm01")) {
        tier_max->x = 2337;
        tier_max->y = 2304;
        *tier_row_cnt = 36;
    }
    else if(!strcmp(gbch, "ibm03")) {
        tier_max->x = 2859;
        tier_max->y = 2816;
        *tier_row_cnt = 44;
    }
    else if(!strcmp(gbch, "ibm04")) {
        tier_max->x = 3330;
        tier_max->y = 3328;
        *tier_row_cnt = 52;
    }
    else if(!strcmp(gbch, "ibm06")) {
        tier_max->x = 3266;
        tier_max->y = 3264;
        *tier_row_cnt = 51;
    }
    else if(!strcmp(gbch, "ibm07")) {
        tier_max->x = 4231;
        tier_max->y = 4224;
        *tier_row_cnt = 66;
    }
    else if(!strcmp(gbch, "ibm08")) {
        tier_max->x = 4414;
        tier_max->y = 4352;
        *tier_row_cnt = 68;
    }
    else if(!strcmp(gbch, "ibm09")) {
        tier_max->x = 4484;
        tier_max->y = 4480;
        *tier_row_cnt = 70;
    }
    else if(!strcmp(gbch, "ibm13")) {
        tier_max->x = 5442;
        tier_max->y = 5440;
        *tier_row_cnt = 85;
    }
    else if(!strcmp(gbch, "ibm15")) {
        tier_max->x = 7602;
        tier_max->y = 7552;
        *tier_row_cnt = 118;
    }
    else if(!strcmp(gbch, "ibm18")) {
        tier_max->x = 9490;
        tier_max->y = 9472;
        *tier_row_cnt = 148;
    }
    else {
        printf("Error: %s not found in IBM 3D-IC benchmarks.\n", gbch);
        exit(1);
    }

    tier_max->x += tier_min->x;
    tier_max->y += tier_min->y;

    shrunk_ratio.x = (tier_max->x - tier_min->x) / (grow_pmax.x - grow_pmin.x);
    // shrunk_ratio.y = (tier_max->y - tier_min->y) / (grow_pmax.y -
    // grow_pmin.y);
    shrunk_ratio.y = (prec)(*tier_row_cnt) / (prec)row_cnt;
}

void output_gp_net_hpwl(char *fn) {
    int i = 0;
    FILE *fp = fopen(fn, "w");
    struct FPOS net_hpwl = zeroFPoint;
    struct NET *net = NULL;

    for(i = 0; i < netCNT; i++) {
        net = &netInstance[i];
        net_hpwl.x = net->max_x - net->min_x;
        net_hpwl.y = net->max_y - net->min_y;
        net_hpwl.z = net->max_z - net->min_z;

        if(GP_DIM_ONE) {
            net_hpwl.x *= place_backup.cnt.x * GP_SCAL;
            net_hpwl.y *= place_backup.cnt.y * GP_SCAL;
            net_hpwl.z *= place_backup.cnt.z * GP_SCAL;
        }

        fprintf(fp, "%s %.6lf %.6lf %.6lf\n", net->name, net_hpwl.x, net_hpwl.y,
                net_hpwl.z);
    }
}

/*
void output_net_hpwl(char *fn) {
    int i = 0;
    FILE *fp = fopen(fn, "w");
    struct FPOS net_hpwl = zeroFPoint;
    struct NET *net = NULL;

    for(i = 0; i < netCNT; i++) {
        net = &netInstance[i];
        net_hpwl.x = net->max_x - net->min_x;
        net_hpwl.y = net->max_y - net->min_y;
        net_hpwl.z = net->max_z - net->min_z;

        fprintf(fp, "%s %.6lf %.6lf %.6lf\n", net->name, net_hpwl.x, net_hpwl.y,
                net_hpwl.z);
    }
}
*/

void output_final_pl(char *output) {
    int i = 0;
    FILE *fp = fopen(output, "w");
    struct MODULE *curModule = NULL;
    struct TERM *curTerminal = NULL;

    fputs("UCLA pl 1.0 \n", fp);
    fputs("# Created	:	Jan  6 2005\n", fp);
    fputs(
        "# User   	:	Gi-Joon Nam & Mehmet Yildiz at IBM Austin "
        "Research({gnam, mcan}@us.ibm.com)\n",
        fp);
    fputs("\n", fp);

    for(i = 0; i < moduleCNT; i++) {
        curModule = &moduleInstance[i];

        if(flg_3dic) {
            if(curModule->flg == StdCell) {
                fprintf(fp, "%s\t%d\t%d\t%d\t: N\n", curModule->name,
                        (int)(curModule->pmin.x + 0.5),
                        (int)(curModule->pmin.y + 0.5),
                        (int)(curModule->pmin.z + 0.5));
            }
            else {
                fprintf(fp, "%s\t%d\t%d\t%d\t: N\n", curModule->name,
                        curModule->pmin_lg.x, curModule->pmin_lg.y,
                        curModule->pmin_lg.z);
            }
        }
        else {
            if(curModule->flg == StdCell) {
                fprintf(fp, "%s\t%d\t%d\t: N\n", curModule->name,
                        (int)(curModule->pmin.x + 0.5),
                        (int)(curModule->pmin.y + 0.5));
            }
            else {
                fprintf(fp, "%s\t%d\t%d\t: N\n", curModule->name,
                        curModule->pmin_lg.x, curModule->pmin_lg.y);
            }
        }
    }

    for(i = 0; i < terminalCNT; i++) {
        if(flg_3dic) {
            curTerminal = &terminalInstance[i];
            if(!curTerminal->isTerminalNI) {
                fprintf(fp, "%s %d\t%d\t%d\t: N /FIXED\n", curTerminal->name,
                        prec2int(curTerminal->pmin.x),
                        prec2int(curTerminal->pmin.y),
                        prec2int(curTerminal->pmin.z));
            }
            else {
                fprintf(fp, "%s %d\t%d\t%d\t: N /FIXED_NI\n", curTerminal->name,
                        prec2int(curTerminal->pmin.x),
                        prec2int(curTerminal->pmin.y),
                        prec2int(curTerminal->pmin.z));
            }
        }
        else {
            curTerminal = &terminalInstance[i];
            if(!curTerminal->isTerminalNI) {
                fprintf(fp, "%s %d\t%d\t: N /FIXED\n", curTerminal->name,
                        prec2int(curTerminal->pmin.x),
                        prec2int(curTerminal->pmin.y));
            }
            else {
                fprintf(fp, "%s %d\t%d\t: N /FIXED_NI\n", curTerminal->name,
                        prec2int(curTerminal->pmin.x),
                        prec2int(curTerminal->pmin.y));
            }
        }
    }
    fclose(fp);
}

/*
void output_hgraph_nofix(char *fn) {
    int i = 0, j = 0;
    int moduleID = 0, netCNTinObject = 0;
    FILE *fp = fopen(fn, "w");
    struct NET *net = NULL;
    struct MODULE *modu = NULL;

    for(i = 0; i < netCNT; i++) {
        net = &netInstance[i];
        if(net->pinCNTinObject > 1)
            netCNTinObject++;
    }

    fprintf(fp, "%d %d 10\n", netCNTinObject, moduleCNT);

    for(i = 0; i < netCNT; i++) {
        net = &netInstance[i];

        if(net->pinCNTinObject <= 1)
            continue;

        for(j = 0; j < net->pinCNTinObject; j++) {
            moduleID = net->pin[j]->moduleID;

            fprintf(fp, "%d ", moduleID);
        }
        fprintf(fp, "\n");
    }

    for(i = 0; i < moduleCNT; i++) {
        modu = &moduleInstance[i];
        fprintf(fp, "%d\n", (int)(modu->area + 0.5));
    }

    fclose(fp);
}
*/

void write_3d_bench(void) {
    int i = 0, j = 0, moduleID = 0, pinIDinModule = 0;

    char mkdir_cmd[BUF_SZ], bch_dir[BUF_SZ];
    sprintf(bch_dir, "%s/3d", dir_bnd);
    sprintf(mkdir_cmd, "mkdir -p %s", bch_dir);

    system(mkdir_cmd);

    char fn_aux[BUF_SZ], fn_nets[BUF_SZ], fn_nodes[BUF_SZ], fn_pl[BUF_SZ],
        fn_wts[BUF_SZ], fn_scl[BUF_SZ];
    sprintf(fn_aux, "%s/%s.aux", bch_dir, gbch);
    sprintf(fn_nets, "%s/%s.nets", bch_dir, gbch);
    sprintf(fn_nodes, "%s/%s.nodes", bch_dir, gbch);
    sprintf(fn_pl, "%s/%s.pl", bch_dir, gbch);
    sprintf(fn_wts, "%s/%s.wts", bch_dir, gbch);
    sprintf(fn_scl, "%s/%s.scl", bch_dir, gbch);

    printf(
        "3D-BCH DIR %s AUX %s.aux NETS %s.nets \
NODES %s.nodes PL %s.pl WTS %s.pl SCL %s.scl\n",
        bch_dir, gbch, gbch, gbch, gbch, gbch, gbch);

    FILE *fp_aux = fopen(fn_aux, "w");
    FILE *fp_nets = fopen(fn_nets, "w");
    FILE *fp_nodes = fopen(fn_nodes, "w");
    FILE *fp_pl = fopen(fn_pl, "w");
    FILE *fp_wts = fopen(fn_wts, "w");
    FILE *fp_scl = fopen(fn_scl, "w");

    struct PIN *pin = NULL;
    struct NET *net = NULL;
    struct MODULE *curModule = NULL;
    struct TERM *curTerminal = NULL;
    struct ROW *rwp = NULL;

    char io[BUF_SZ];

    fprintf(fp_aux,
            "RowBasedPlacement :  %s.nodes  %s.nets  %s.wts  %s.pl  %s.scl\n",
            gbch, gbch, gbch, gbch, gbch);

    fclose(fp_aux);

    fputs("UCLA wts 1.0\n", fp_wts);
    fputs("# Created\t:\tDec 27 2004\n", fp_wts);
    fputs(
        "# User   \t:\tGi-Joon Nam & Mehmet Yildiz at IBM Austin "
        "Research({gnam, mcan}@us.ibm.com)\n",
        fp_wts);

    fclose(fp_wts);

    fputs("UCLA nets 1.0 \n", fp_nets);
    fputs("# Created	:	Jan  6 2005\n", fp_nets);
    fputs(
        "# User   	:	Gi-Joon Nam & Mehmet Yildiz at IBM Austin "
        "Research({gnam, mcan}@us.ibm.com)\n",
        fp_nets);

    fputs("\n", fp_nets);

    fprintf(fp_nets, "NumNets : %d\n", netCNT);
    fprintf(fp_nets, "NumPins : %d\n", pinCNT);

    fputs("\n", fp_nets);

    for(i = 0; i < netCNT; i++) {
        net = &netInstance[i];

        fprintf(fp_nets, "NetDegree : %d   %s\n", net->pinCNTinObject,
                net->name);
        for(j = 0; j < net->pinCNTinObject; j++) {
            pin = net->pin[j];

            moduleID = pin->moduleID;
            pinIDinModule = pin->pinIDinModule;

            if(pin->IO == 0)
                strcpy(io, "I");
            else if(pin->IO == 1)
                strcpy(io, "O");
            else  // B
                strcpy(io, "B");

            if(pin->term == 0) {
                curModule = &moduleInstance[moduleID];

                pinIDinModule = pin->pinIDinModule;

                fprintf(fp_nets, "\t%s\t%s : %.6lf\t%.6lf\t%.6lf\n",
                        curModule->name, io, curModule->pof[pinIDinModule].x,
                        curModule->pof[pinIDinModule].y,
                        curModule->pof[pinIDinModule].z);
            }
            else {
                moduleID = pin->moduleID;
                curTerminal = &terminalInstance[moduleID];

                fprintf(fp_nets, "\t%s\t%s : %.6lf\t%.6lf\t%.6lf\n",
                        curTerminal->name, io,
                        curTerminal->pof[pinIDinModule].x,
                        curTerminal->pof[pinIDinModule].y,
                        curTerminal->pof[pinIDinModule].z);
            }
        }
    }

    fclose(fp_nets);

    fputs("UCLA nodes 1.0\n", fp_nodes);
    fputs("# Created\t:\tJan  6 2005\n", fp_nodes);
    fputs(
        "# User   \t:\tGi-Joon Nam & Mehmet Yildiz at IBM Austin "
        "Research({gnam, mcan}@us.ibm.com)\n",
        fp_nodes);

    fputs("\n", fp_nodes);

    fprintf(fp_nodes, "NumNodes :  \t%d\n", moduleCNT + terminalCNT);
    fprintf(fp_nodes, "NumTerminals :  \t%d\n", terminalCNT);

    for(i = 0; i < moduleCNT; i++) {
        curModule = &moduleInstance[i];

        fprintf(fp_nodes, "\t%s\t%d\t%d\t%d\n", curModule->name,
                (int)(curModule->size.x + 0.5), (int)(curModule->size.y + 0.5),
                (int)(curModule->size.z + 0.5));
    }

    for(i = 0; i < terminalCNT; i++) {
        curTerminal = &terminalInstance[i];

        fprintf(fp_nodes, "\t%s\t%d\t%d\t%d\tterminal\n", curTerminal->name,
                (int)(curTerminal->size.x + 0.5),
                (int)(curTerminal->size.y + 0.5),
                (int)(curTerminal->size.z + 0.5));
    }

    fclose(fp_nodes);

    fputs("UCLA pl 1.0\n", fp_pl);
    fputs("# Created\t:\tJan  6 2005\n", fp_pl);
    fputs(
        "# User   \t:\tGi-Joon Nam & Mehmet Yildiz at IBM Austin "
        "Research({gnam, mcan}@us.ibm.com)\n",
        fp_pl);

    fputs("\n", fp_pl);

    for(i = 0; i < moduleCNT; i++) {
        curModule = &moduleInstance[i];

        fprintf(fp_pl, "%s\t0\t0\t0\t: N\n", curModule->name);
    }

    for(i = 0; i < terminalCNT; i++) {
        curTerminal = &terminalInstance[i];

        fprintf(fp_pl, "%s\t%d\t%d\t%d\t: N /FIXED\n", curTerminal->name,
                (int)(curTerminal->pmin.x + 0.5),
                (int)(curTerminal->pmin.y + 0.5),
                (int)(curTerminal->pmin.z + 0.5));
    }

    fclose(fp_pl);

    fputs("UCLA scl 1.0\n", fp_scl);
    fputs("# Created\t:\tJan  6 2005\n", fp_scl);
    fputs(
        "# User   \t:\tGi-Joon Nam & Mehmet Yildiz at IBM Austin "
        "Research({gnam, mcan}@us.ibm.com)\n",
        fp_scl);

    fputs("\n", fp_scl);

    int z = 0;
    int tot_row_cnt = 0;
    struct TIER *tier = NULL;

    for(z = 0; z < numLayer; z++) {
        tot_row_cnt += tier_st[z].row_cnt;
    }

    fprintf(fp_scl, "NumRows :  \t%d\n", tot_row_cnt);
    fprintf(fp_scl, "NumTiers :  \t%d\n", numLayer);
    fputs("\n", fp_scl);

    for(z = 0; z < numLayer; z++) {
        tier = &tier_st[z];
        for(i = 0; i < tier->row_cnt; i++) {
            rwp = &tier->row_st[i];

            fprintf(fp_scl, "CoreRow Horizontal\n");
            fprintf(fp_scl, "  Coordinate    :   %d\n", rwp->pmin.y);
            fprintf(fp_scl, "  Height        :   %d\n", rwp->size.y);
            fprintf(fp_scl, "  CoordinateZ   :   %d\n", rwp->pmin.z);
            fprintf(fp_scl, "  HeightZ       :   %d\n", rwp->size.z);
            fprintf(fp_scl, "  Sitewidth     :    1\n");
            fprintf(fp_scl, "  Sitespacing   :    1\n");
            fprintf(fp_scl, "  Siteorient    :    1\n");
            fprintf(fp_scl, "  Sitesymmetry  :    1\n");
            fprintf(fp_scl, "  SubrowOrigin  :    %d\tNumSites  :  %d\n",
                    rwp->pmin.x, rwp->x_cnt);
            fprintf(fp_scl, "End\n");
        }
    }
    fclose(fp_scl);
}

void extract_dir(char *f, char *d) {
    int i = 0, idx = 0;
    strcpy(d, f);

    int len = strlen(f);

    for(i = 0; i < len; i++) {
        if(f[i] == '/') {
            idx = i + 1;
            continue;
        }
    }
    d[idx] = '\0';
}

void extract_sfx(char *f, char *s) {
    int i = 0, last_dot = 0;
    int len = strlen(f);

    for(i = 0; i < len; i++) {
        if(f[i] == '.') {
            last_dot = i;
            continue;
        }
    }

    strcpy(s, &f[last_dot + 1]);
}

