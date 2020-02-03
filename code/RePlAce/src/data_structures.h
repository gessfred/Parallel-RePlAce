#include <string>
using std::string;
using std::min;
using std::max;
using std::cout;
using std::endl;
using std::fixed;
using std::setprecision;
using std::vector;
typedef float prec;
struct POS;
struct fpos2_t;
struct pos_t;
struct net_t;
struct net_t;
struct pin_t;
struct FPOS {
    prec x;
    prec y;
    prec z;

    FPOS() {
        SetZero();
    };
    
    FPOS(prec _x, prec _y) : x(_x), y(_y), z(0) {};
    FPOS(prec _x, prec _y, prec _z) : x(_x), y(_y), z(_z) {};
    
    inline void Set(prec a) {
        x = y = z = a;
    }
    inline void Set(FPOS a) {
        x = a.x;
        y = a.y;
        z = a.z;
    }
    inline void Set(prec _x, prec _y, prec _z) {
        x = _x; y = _y; z = _z;
    }

    inline void Set(POS a);

    inline void SetZero() {
        x = y = z = 0.0f;
    }
    inline void Add(FPOS a) {
        x += a.x;
        y += a.y;
        z += a.z;
    }
    inline void SetAdd(FPOS a, FPOS b) {
        x = a.x + b.x;
        y = a.y + b.y;
        z = a.z + b.z;
    }
    inline void Min(FPOS a) {
        x = min(x, a.x);
        y = min(y, a.y);
        z = min(z, a.z);
    }
    inline void SetMin(FPOS a, FPOS b) {
        x = min(a.x, b.x);
        y = min(a.y, b.y);
        z = min(a.z, b.z);
    }
    inline void Max(FPOS a) {
        x = max(x, a.x);
        y = max(y, a.y);
        z = max(z, a.z);
    }
    inline void SetMax(FPOS a, FPOS b) {
        x = max(a.x, b.x);
        y = max(a.y, b.y);
        z = max(a.z, b.z);
    }

    inline prec GetProduct() {
        return x*y*z;
    }
    inline void Dump() {
        cout << "(" << x << " " << y << " " << z << ")" << endl;
    }  
    inline void Dump(string a) {
        cout << a << ": (" << x << " " << y << " " << z << ")" << endl;   
    }

    inline void from(fpos2_t fp);
};

struct POS {
    int x;
    int y;
    int z;

    POS() {
        SetZero();
    };

    POS(int _x, int _y) : x(_x), y(_y), z(0) {};
    POS(int _x, int _y, int _z) : x(_x), y(_y), z(_z) {};

    inline void Set(int a) {
        x = y = z = a;
    }
    inline void Set(POS a) {
        x = a.x;
        y = a.y;
        z = a.z;
    }
    
    inline void Set(int _x, int _y, int _z) {
        x = _x; y = _y; z = _z;
    }

    inline void Set(FPOS fp);

    inline void SetZero() {
        x = y = z = 0.0f;
    }
    inline void Add(POS a) {
        x += a.x;
        y += a.y;
        z += a.z;
    }

    inline void SetAdd(POS a, POS b) {
        x = a.x + b.x;
        y = a.y + b.y;
        z = a.z + b.z;
    }
    inline void Min(POS a) {
        x = min(x, a.x);
        y = min(y, a.y);
        z = min(z, a.z);
    }
    inline void SetMin(POS a, POS b) {
        x = min(a.x, b.x);
        y = min(a.y, b.y);
        z = min(a.z, b.z);
    }
    inline void Max(POS a) {
        x = max(x, a.x);
        y = max(y, a.y);
        z = max(z, a.z);
    }
    inline void SetMax(POS a, POS b) {
        x = max(a.x, b.x);
        y = max(a.y, b.y);
        z = max(a.z, b.z);
    }
    inline int GetProduct() {
        return x*y*z;
    }
    inline void SetXProjection(int a, int b) {
        x = (x<a)? a : (x>b)? b : x;
    }
    inline void SetYProjection(int a, int b) {
        y = (y<a)? a : (y>b)? b : y;
    }
    inline void SetZProjection(int a, int b) {
        z = (z<a)? a : (z>b)? b : z;
    }
    inline void SetProjection(POS a, POS b) {
        SetXProjection(a.x, b.x);
        SetYProjection(a.y, b.y);
        SetZProjection(a.z, b.z);
    }
    inline void SetXYProjection(POS a, POS b) {
        SetXProjection(a.x, b.x);
        SetYProjection(a.y, b.y);
    }
    inline void Dump() {
        cout << "(" << x << " " << y << " " << z << ")" << endl;
    }   
    inline void Dump(string a) {
        cout << a << ": (" << x << " " << y << " " << z << ")" << endl;   
    }

    inline void from(pos_t p);
};

struct RECT {
    FPOS pmin;
    FPOS pmax;
    RECT() { pmin.SetZero(); pmax.SetZero(); };
    void Dump() {
        pmin.Dump("RECT: pmin");
        pmax.Dump("RECT: pmax");
        cout << endl;
    }
};

// Pin Instance
struct PIN {
    FPOS fp;
    FPOS e1;
    FPOS e2;
    POS flg1;
    POS flg2;
    int moduleID;
    int pinIDinModule;
    int netID;
    int pinIDinNet;
    int gid;   // current Pin's idx
    int IO;    // I -> 0; O -> 1
    int term;  // term -> 1, move -> 0
    int tier;
    int X_MIN;
    int Y_MIN;
    int X_MAX;
    int Y_MAX;
    int netCNT;
    inline void from(pin_t pin);
};

// *.nodes -> not isTerminal
// Module Instance
struct MODULE {
    FPOS pmin;
    FPOS pmax;
    FPOS size;
    FPOS half_size;
    FPOS center;
    FPOS *pof;
    PIN **pin;
    prec area;
    char name[255];
    int idx;
    int netCNTinObject;
    int pinCNTinObject;
    int flg;
    int tier;
    int mac_idx;
    int ovlp_flg;
    POS pmin_lg;
    POS pmax_lg;

    MODULE() : pof(0), pin(0), area(0.0f), name(""), 
    idx(0), netCNTinObject(0), pinCNTinObject(0),
    flg(0), tier(0), mac_idx(0), ovlp_flg(0) {
        pmin.SetZero();
        pmax.SetZero();
        size.SetZero();
        half_size.SetZero();
        center.SetZero();
        pmin_lg.SetZero();
        pmax_lg.SetZero();
    }
    void Dump(string a) {
        cout << fixed << setprecision(0) << a << endl;
        cout << idx << ", name: " << name << endl;
        cout << "tier: " << tier << endl;
        cout << "mac_idx: " << mac_idx << endl;
        cout << "ovlp_flg: " << ovlp_flg << endl;
        pmin.Dump("pmin");
        pmax.Dump("pmax");
        size.Dump("size");
        half_size.Dump("half_size");
        center.Dump("center");
        cout << "area: " << area << endl;
        cout << "netCNTinObject: " << netCNTinObject << endl;
        cout << "pinCNTinObject: " << pinCNTinObject << endl;
        pmin_lg.Dump("pmin_lg");
        pmax_lg.Dump("pmax_lg");
        cout << endl;
    }
};

// *.nodes -> isTerminal // isTerminalNI
// Terminal Instance
struct TERM {
    // struct  POS p;   // left_x,  low_y
    // struct  POS mp;  // right_x, up_y
    FPOS pmin;
    FPOS pmax;
    prec area;
    FPOS size;
    FPOS center;
    FPOS *pof;
    PIN **pin;
    char name[255];
    int idx;
    int netCNTinObject;
    int pinCNTinObject;
    int IO;    // I -> 0; O -> 1
    int tier;
    bool isTerminalNI;
    
    prec PL_area;
    
    TERM() : area(0.0f), pof(0), pin(0), name(""), 
    idx(0), netCNTinObject(0), pinCNTinObject(0),
    IO(0), tier(0), isTerminalNI(0), PL_area(0.0f)  {
        pmin.SetZero();
        pmax.SetZero();
        size.SetZero();
        center.SetZero();
    }

    void Dump() {
        printf("terminal[%d]: name: %s \n",
               idx, name);
        fflush(stdout);
        cout << fixed << setprecision(0) ;
        cout << "isTerminalNI: " << (isTerminalNI? "YES" : "NO") << endl;
        cout << "IO: " << ((IO==0)? "Input" : "Output") << endl;
        cout << "tier: " << tier << endl;
        pmin.Dump("pmin"); 
        pmax.Dump("pmax"); 
        cout << "area: " << area << endl;
        size.Dump("size");
        center.Dump("center");
        cout << "netCNTinObject: " << netCNTinObject << endl;
        cout << "pinCNTinObject: " << pinCNTinObject << endl;
        cout << "PL_area: " << PL_area << endl;
        cout << endl;
    }
};

struct CELLx {
    FPOS pmin;
    FPOS pmax;
    FPOS den_pmin;
    FPOS den_pmax;
    FPOS center;
    prec area;
    FPOS *pof;
    PIN **pin;
    prec x_grad;
    prec y_grad;
    int tier;
    FPOS size;
    prec den_scal;
    FPOS half_size;
    FPOS half_den_size;
    char name[255];
    int idx;
    int pinCNTinObject;
    int netCNTinObject;
    int flg;
    prec area_org_befo_bloating;
    prec inflatedNewArea;
    prec inflatedNewAreaDelta;
    prec den_scal_org_befo_bloating;
    FPOS size_org_befo_bloating;
    FPOS half_size_org_befo_bloating;
    FPOS half_den_size_org_befo_bloating;
    FPOS *pof_tmp;
    PIN **pin_tmp;
};

class UFPin {
   public:
    // int id;
    int parent;
    int rank;
    int modu;

    UFPin() {
    }

    UFPin(int moduleID) {
        // id = 0;
        parent = moduleID;
        rank = 0;
        modu = moduleID;
    }

    ~UFPin() {
    }
};

class TwoPinNets {
   public:
    bool selected;
    int start_modu;
    int end_modu;
    int rect_dist;
    int i;
    int j;

    TwoPinNets() {
    }

    TwoPinNets(int start_modu, int end_modu, prec rect_dist, int i, int j) {
        selected = false;
        this->start_modu = start_modu;
        this->end_modu = end_modu;
        this->rect_dist = rect_dist;
        this->i = i;
        this->j = j;
    }

    ~TwoPinNets() {
    }
};

bool TwoPinNets_comp(TwoPinNets x, TwoPinNets y);

int UFFind(struct NET *net, int moduleID);
void UFUnion(struct NET *net, int idx, int idy);

class ROUTRACK {
   public:
    struct FPOS from;
    struct FPOS to;
    int layer;  // 1:M_Layer1, 2:M_Layer2, ..., etc.
    int netIdx;
    ROUTRACK(struct FPOS f, struct FPOS t, int lay, int idx) {
        from.x = f.x;
        from.y = f.y;
        from.z = f.z;
        to.x = t.x;
        to.y = t.y;
        to.z = t.z;
        layer = lay;
        netIdx = idx;
    };
};


struct NET {
    char name[255];
    std::map< int, UFPin > mUFPin;
    vector< TwoPinNets > two_pin_nets;
    vector< ROUTRACK > routing_tracks;

    prec min_x;
    prec min_y;
    prec min_z;
    prec max_x;
    prec max_y;
    prec max_z;
    FPOS sum_num1;
    FPOS sum_num2;
    FPOS sum_denom1;
    FPOS sum_denom2;
    PIN **pin;
    PIN **pin2;
    FPOS terminalMin;
    FPOS terminalMax;
    prec hpwl_x;
    prec hpwl_y;
    prec hpwl_z;
    prec hpwl;
    int pinCNTinObject;
    int pinCNTinObject2;
    int pinCNTinObject_tier; // used for writing bookshelf
    int HPWL;
    int idx;
    int mod_idx;
    int flg;
    int tier;
    prec stn_cof;  // lutong
    prec wl_rsmt;  // lutong
    //inline void from(net_t net);
    inline void copy(net_t* net);
};
