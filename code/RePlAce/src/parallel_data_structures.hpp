#include <bitset>
#include <vector>

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

//From "data_structures.h"
struct FPOS;
struct POS;
struct CELLx;
struct PIN;
struct NET;
struct pos_t {
    public:
    int x;
    int y;
    inline void set(int x_, int y_) {
        x = x_;
        y = y_;
    }
    inline void Legalize(pos_t lowerBound, pos_t upperBound) {
        if(x < lowerBound.x) x = 0;
        if(x > upperBound.x - 1) x = upperBound.x - 1;
        if(y < lowerBound.y) y = 0;
        if(y > upperBound.y - 1) y = upperBound.y - 1;
    }
};


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


//pins as a charge 
struct charge_t {
    fpos2_t e1;
    fpos2_t e2;
    fpos2_t fp;
    int idx;
    int moduleID;
    int netID;
    bitset<9> meta;// flg2(x, y) flg1(x, y) term
    
    void from(PIN*);
    void to(PIN*);
};

//for density area_share...
struct cell_t {
    pos_t b0;
    pos_t b1;
    fpos2_t den_pmin;
    fpos2_t den_pmax;
    fpos2_t pmin;
    fpos2_t pmax;
    fpos2_t half_size;
    float scale;
    char flg;//()
    void from(CELLx*);
    void to(CELLx*);    
};

struct field_t{
    charge_t* pin;
    fpos2_t sum_num1;
    fpos2_t sum_num2;
    fpos2_t sum_denom1;
    fpos2_t sum_denom2;
    fpos2_t min;
    //fpos2_t terMin;
    //fpos2_t terMax;
    fpos2_t max;
    int pinCNT;
 inline void from(NET* net); 
                        
 inline void copy(NET* net); 
};

//cell as a collection of pins
struct io_t{
    charge_t** pins;
    int pinCNT;
    char flg;
inline void from(CELLx* cell, field_t* nets);    
};

//pin as a position
struct pin_t {
    fpos2_t fp;
    int moduleID;
    bitset<5> meta;// Y_MAX | X_MAX | Y_MIN | X_MIN | term
 inline void from(PIN* pin);
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
inline void Copy(cell_t* origin); 
    inline void Copyback(cell_t* destination);
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

vector<int> refIo(io_t* cells, size_t numberOfCells);
vector<int> refNets(field_t* nets, size_t numberOfNets);
vector<int> refCells(Cell_t* cells);
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