#include "parallel_data_structures.hpp"

vector<int> refIo(io_t* cells, size_t numberOfCells) {
    vector<int> cellsRefs;
    for(size_t i = 0; i < numberOfCells; ++i) 
        cellsRefs.push_back(cells[i].pinCNT);
    return cellsRefs;
}

vector<int> refNets(field_t* nets, size_t numberOfNets) {
    vector<int> cellsRefs;
    for(size_t i = 0; i < numberOfNets; ++i) 
        cellsRefs.push_back(nets[i].pinCNT);
    return cellsRefs;
}

vector<int> refCells(Cell_t* cells) {
    vector<int> cellsRefs;
    for(size_t i = 0; i < cells->size; ++i) 
        cellsRefs.push_back((cells->b1_x[i]-cells->b0_x[i])*(cells->b1_y[i]-cells->b0_y[i]));
    return cellsRefs;
}
