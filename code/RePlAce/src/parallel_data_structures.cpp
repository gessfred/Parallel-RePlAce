///////////////////////////////////////////////////////////////////////////////
// Author: Frederic Gessler
//          (supervisor: Dr. Mirjana Stojilovic),
//          based on Dr. Jingwei Lu with ePlace and ePlace-MS
//          and code by Ilgweon Kang and Lutong Wang
//
//          Many subsequent improvements were made by Mingyu Woo
//          leading up to the initial release.
//
// BSD 3-Clause License
//
// Copyright (c) 2019, The Regents of EPFL
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

#include "parallel_data_structures.hpp"

vector<int> refIo(cell_phy_t* cells, size_t numberOfCells) {
    vector<int> cellsRefs;
    for(size_t i = 0; i < numberOfCells; ++i) 
        cellsRefs.push_back(cells[i].pinCNT);
    return cellsRefs;
}

vector<int> refNets(net_t* nets, size_t numberOfNets) {
    vector<int> cellsRefs;
    for(size_t i = 0; i < numberOfNets; ++i) 
        cellsRefs.push_back(nets[i].pinCNT);
    return cellsRefs;
}

vector<int> refCells(cell_den_t* cells, size_t numberOfCells) {
    vector<int> cellsRefs;
    for(size_t i = 0; i < numberOfCells; ++i) 
        cellsRefs.push_back((cells[i].binEnd.y-cells[i].binStart.y)*(cells[i].binEnd.x-cells[i].binStart.x));
    return cellsRefs;
}

