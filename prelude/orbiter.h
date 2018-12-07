typedef union {
    double data[3];               ///< array data interface
    struct { double x, y, z; };   ///< named data interface
} VECTOR3;

inline VECTOR3 _V(double x, double y, double z) {     VECTOR3 vec = {x,y,z}; return vec; }
