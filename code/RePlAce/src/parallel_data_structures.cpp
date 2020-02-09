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

