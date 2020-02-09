#include <bitset>
#include <vector>
#include <map>

using namespace std;

#define __getMeta__()   meta[0] = pin->term; \
                        meta[1] = pin->X_MIN; \
                        meta[2] = pin->X_MAX; \
                        meta[3] = pin->Y_MIN; \
                        meta[4] = pin->Y_MAX; \
                        meta[5] = pin->flg1.x;\
                        meta[6] = pin->flg1.y;\
                        meta[7] = pin->flg2.x;\
                        meta[8] = pin->flg2.y;

/**
Original data structures
From "data_structures.h"
*/
struct FPOS;
struct POS;
struct CELLx;
struct PIN;
struct NET;
struct BIN;
struct MODULE;

/**
pos2_t represents a 2d integer point
*/
struct pos2_t {
    public:
    int x;
    int y;
    inline void set(int x_, int y_) {
        x = x_;
        y = y_;
    }
    inline void Legalize(pos2_t lowerBound, pos2_t upperBound) {
        if(x < lowerBound.x) x = 0;
        if(x > upperBound.x - 1) x = upperBound.x - 1;
        if(y < lowerBound.y) y = 0;
        if(y > upperBound.y - 1) y = upperBound.y - 1;
    }
};

/**
fpos2_t represents a 2d real point
*/
struct fpos2_t {
    public:
    float x;
    float y;
    inline void set(float x_, float y_) {
        x = x_;
        y = y_;
    }
    void from(FPOS);
} ;


/**
pin_t represents a pin as a charge
*/
struct pin_t {
    fpos2_t expL;
    fpos2_t expR;
    fpos2_t coord;
    int pinID;
    int moduleID;
    int netID;
    bitset<9> meta;// flg2(x, y) flg1(x, y) term
    
    void from(PIN*);
    void to(PIN*);
};

/**
cell_den_t represents a cell as a rectangle
*/
struct cell_den_t {
    pos2_t binStart;
    pos2_t binEnd;
    fpos2_t min;
    fpos2_t max;
    fpos2_t size; // half_size
    float scale;
    char type;//()
    void from(CELLx*);
    void to(CELLx*);    
};



struct net_t {
    pin_t* pinArray;
    fpos2_t sumNumL;
    fpos2_t sumNumR;
    fpos2_t sumDenomL;
    fpos2_t sumDenomR;
    fpos2_t min;
    fpos2_t max;
    int pinCNT;

    inline void from(NET* net); 
    inline void copy(NET* net); 
};

//cell as a collection of pins
struct cell_phy_t {
    pin_t** pinArrayPtr;
    int pinCNT;
    char type;
    
    inline void from(CELLx* cell, net_t* nets);    
};

//Represent all the cells in of itself
struct Cell_t {
    float* center_x;
    float* center_y;
    float* den_pmin_x;
    float* den_pmin_y;
    float* den_pmax_x;
    float* den_pmax_y;
    float* pmin_x;
    float* pmin_y;
    float* pmax_x;
    float* pmax_y;
    float* half_size_x;
    float* half_size_y;
    float* scale;
    char* flg;
    float* b0_x;
    float* b1_x;
    float* b0_y;
    float* b1_y;
    size_t size;

    inline void build(size_t N) {
        center_x = (float*)calloc(sizeof(float), N);
        center_y = (float*)calloc(sizeof(float), N);
        den_pmin_x = (float*)calloc(sizeof(float), N);
        den_pmin_y = (float*)calloc(sizeof(float), N);
        den_pmax_x = (float*)calloc(sizeof(float), N);
        den_pmax_y = (float*)calloc(sizeof(float), N);
        pmin_x = (float*)calloc(sizeof(float), N);
        pmin_y = (float*)calloc(sizeof(float), N);
        pmax_x = (float*)calloc(sizeof(float), N);
        pmax_y = (float*)calloc(sizeof(float), N);
        half_size_x = (float*)calloc(sizeof(float), N);
        half_size_y = (float*)calloc(sizeof(float), N);
	    scale = (float*)calloc(sizeof(float), N);
        flg = (char*)calloc(sizeof(char), N);
        b0_x = (float*)calloc(sizeof(float), N);
        b0_y = (float*)calloc(sizeof(float), N);
        b1_x = (float*)calloc(sizeof(float), N);
        b1_y = (float*)calloc(sizeof(float), N);
        size = N;
    }

    inline void copy(CELLx* origin); 
    inline void copyback(CELLx* destination);
    //inline void Copy(cell_t* origin); 
    //inline void Copyback(cell_t* destination);
    inline void destroy() {
        free(center_x);
        free(center_y);
        free(den_pmin_x);
        free(den_pmin_y);
        free(den_pmax_x);
        free(den_pmax_y);
        free(pmin_x);
        free(pmin_y);
        free(pmax_x);
        free(pmax_y);
        free(half_size_x);
        free(half_size_y);
        free(flg);
        free(b0_x);
        free(b0_y);
        free(b1_x);
        free(b1_y);
        size = 0;
    }
};

struct timing_t {
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

struct bin_t {
    fpos2_t max;
    fpos2_t min;
    fpos2_t field;
    float cellArea;
    float fillerArea;
    float potential;
    inline void from(BIN* bin);
    inline void to(BIN* bin);
};
//bins as a collection of area data
struct area_t {
    float virtArea;
    long terminArea;
    float binDensity;
    float fillerDensity;
    pos2_t coord;
    inline void from(BIN* bin);
    inline void to(BIN* bin);
};
struct circuit_t {
    net_t* nets;
    cell_den_t* rects; // cell rectangles
    cell_phy_t* cells;
    bin_t* bins;
    area_t* areas;
    fpos2_t** pinOffsets;

    size_t numberOfBins;
    size_t numberOfCells;
    size_t numberOfNets;
    size_t numberOfThreads;

    map<string, float> timing;

    float** cellAreas;
    float** fillerAreas;
    //schedules
    size_t* constPinsPerCell;
    size_t* constPinsPerNet;
    //bin matrix

    public:
    circuit_t () {
    }

    circuit_t(size_t numberOfThreads) {
        this->numberOfThreads = numberOfThreads;
    }

    inline circuit_t* withCells(CELLx* cells, size_t numberOfCells);
    inline circuit_t* withModules(MODULE* modules, size_t numberOfModules);
    inline circuit_t* withNets(NET* nets, size_t numberOfNets);
    inline circuit_t* withBins(BIN* bins, size_t numberOfBins);
    inline void destroy(CELLx* cells, NET* nets, BIN* bins);
};

vector<int> refIo(cell_phy_t* cells, size_t numberOfCells);
vector<int> refNets(net_t* nets, size_t numberOfNets);
vector<int> refCells(cell_den_t* cells, size_t numberOfCells);

inline size_t* schedule(vector<int> data, unsigned int numberOfThreads) {
	size_t tid=0;
	float work=0, workUpTo= 0;
	for(size_t i = 0; i < data.size(); ++i) 
		work+=data[i];
	work /= numberOfThreads;
	size_t* schedule = (size_t*)calloc(numberOfThreads, sizeof(size_t));
	for(size_t i = 0; i < data.size(); ++i) {
		workUpTo+=data[i];
		if(workUpTo>work) {
			schedule[tid++] = i;
			workUpTo=0;
		}
	}
	return schedule;	
}
