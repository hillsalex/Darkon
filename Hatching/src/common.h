/**
  Vector constructs and utilities.

  author - psastras

**/

#pragma once

#include <iostream>
#include <math.h>

#define GL_GLEXT_LEGACY
#define GL_GLEXT_PROTOTYPES

#include <glm/glm.h>
#include <qgl.h>
#include "glext.h"

#ifdef WIN32
#undef near
#undef far
#endif

#define SAFE_DELETE(x) if((x)) { delete (x); (x) = NULL; }
#define MAX(x, y) (x) > (y) ? (x) : (y)
#define MIN(x, y) (x) < (y) ? (x) : (y)
#define PI 3.14159265358979323846


struct float3 {


    float3(float v0 = 0, float v1 = 0, float v2 = 0) : x(v0), y(v1), z(v2) { }
    float3(float *data) { x = data[0]; y = data[1]; z = data[2]; }

    static inline float3 zero() { return float3(0,0,0); }

    #define VECOP_PCW(op) { x op rhs.x; y op rhs.y; z op rhs.z; return *this; }
    #define VECOP_SCA(op) { x op rhs;   y op rhs  ; z op rhs  ; return *this; }

    inline float3& operator  = (const float3& rhs) VECOP_PCW( =) /// equality assignment
    inline float3& operator += (const float3& rhs) VECOP_PCW(+=) /// piecewise addition operator
    inline float3& operator -= (const float3& rhs) VECOP_PCW(-=) /// piecewise subtraction operator


    inline float3  operator  + (const float3& rhs) const { return float3(*this) += rhs; } /// piecewise addition
    inline float3  operator  - (const float3& rhs) const { return float3(*this) -= rhs; } /// piecewise subtraction

    inline float3& operator += (const float  rhs)  VECOP_SCA(+=) /// scalar addition operator
    inline float3& operator -= (const float  rhs)  VECOP_SCA(-=) /// scalar subtraction operator
    inline float3& operator *= (const float  rhs)  VECOP_SCA(*=) /// scalar multiplication operator
    inline float3& operator /= (const float  rhs)  VECOP_SCA(/=) /// scalar division operator

    inline float3  operator  + (const float  rhs) const { return float3(*this) += rhs; } /// piecewise addition
    inline float3  operator  - (const float  rhs) const { return float3(*this) -= rhs; } /// piecewise subtraction
    inline float3  operator  * (const float  rhs) const { return float3(*this) *= rhs; } /// piecewise multiplication
    inline float3  operator  / (const float  rhs) const { return float3(*this) /= rhs; } /// piecewise multiplication

    #undef VECOP_PCW
    #undef VECOP_SCA

    inline float dot(const float3 &rhs) const {
            return x * rhs.x + y * rhs.y + z * rhs.z;
    }

    inline float normalize() {
            float m = getMagnitude();
            x /= m, y /= m, z /= m;
            return m;
    }

    inline float3 getNormalized() {
            float m = getMagnitude();
            return float3(x / m, y / m, z / m);
    }

    inline float getMagnitude() const {
            return sqrt(getMagnitude2());
    }

    inline float getMagnitude2() const {
            return x * x + y * y + z * z;
    }

    inline float getDistance(const float3 &rhs) const {
            return sqrt(getDistance2(rhs));
    }

    inline float getDistance2(const float3 &rhs) const {
            return (rhs.x - x) * (rhs.x - x) + (rhs.y - y) * (rhs.y - y) +
                    (rhs.z - z) * (rhs.z - z);
    }

    inline float3 cross(const float3& rhs) const {
            return float3(data[1] * rhs.data[2] - data[2] * rhs.data[1],
                           data[2] * rhs.data[0] - data[0] * rhs.data[2],
                           data[0] * rhs.data[1] - data[1] * rhs.data[0]);
    }

    inline bool operator==(const float3 &rhs) {
            return (x == rhs.x && y == rhs.y && z == rhs.z);
    }

    inline bool operator!=(const float3 &rhs) {
            return (x != rhs.x || y != rhs.y || z != rhs.z);
    }
    /**
       Assuming *this is incident to the surface and the result is pointing
       away from the surface.
    **/
    inline float3 reflectVector(const float3 &normal) {
          return (*this) - (normal * ((*this).dot(normal))) * 2.0;
    }

    union {
        struct {
            float x, y, z;
        };
        struct {
            float r, g, b;
        };
        float data[3];
    };
};

inline static float dot(const float3 &v1, const float3 &v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

inline float3 operator*(const float scale, const float3 &rhs) {
    return float3(rhs.x * scale, rhs.y * scale, rhs.z * scale);
}
inline float3 operator-(const float3 &rhs) {
    return float3(-rhs.x, -rhs.y, -rhs.z);
}

inline std::ostream& operator<<(std::ostream& os, const float3& f) {
        os <<"[";
        for (unsigned i = 0; i < 3; ++i) {
            os << f.data[i] << ",";
        }
        os << "]";
        return os;
}

struct float2 {
    union {
        struct {
            float x, y;
        };
        float data[2];
    };
};

struct Model {
    GLMmodel *model;
    GLuint idx;
};

struct Camera {
    float3 eye, center, up;
    float fovy, near, far;
};

//DOUBLE2 STUFFS PLEASE WORK.
struct double2 {


    double2(double v0 = 0, double v1 = 0) : x(v0), y(v1) { }
    double2(double *data) { x = data[0]; y = data[1]; }

    static inline double2 zero() { return double2(0,0); }

    #define VECOP_PCW(op) { x op rhs.x; y op rhs.y; return *this; }
    #define VECOP_SCA(op) { x op rhs;   y op rhs  ; return *this; }

    inline double2& operator  = (const double2& rhs) VECOP_PCW( =) /// equality assignment
    inline double2& operator += (const double2& rhs) VECOP_PCW(+=) /// piecewise addition operator
    inline double2& operator -= (const double2& rhs) VECOP_PCW(-=) /// piecewise subtraction operator


    inline double2  operator  + (const double2& rhs) const { return double2(*this) += rhs; } /// piecewise addition
    inline double2  operator  - (const double2& rhs) const { return double2(*this) -= rhs; } /// piecewise subtraction

    inline double2& operator += (const double  rhs)  VECOP_SCA(+=) /// scalar addition operator
    inline double2& operator -= (const double  rhs)  VECOP_SCA(-=) /// scalar subtraction operator
    inline double2& operator *= (const double  rhs)  VECOP_SCA(*=) /// scalar multiplication operator
    inline double2& operator /= (const double  rhs)  VECOP_SCA(/=) /// scalar division operator

    inline double2  operator  + (const double  rhs) const { return double2(*this) += rhs; } /// piecewise addition
    inline double2  operator  - (const double  rhs) const { return double2(*this) -= rhs; } /// piecewise subtraction
    inline double2  operator  * (const double  rhs) const { return double2(*this) *= rhs; } /// piecewise multiplication
    inline double2  operator  / (const double  rhs) const { return double2(*this) /= rhs; } /// piecewise multiplication

    #undef VECOP_PCW
    #undef VECOP_SCA


    inline double magnitude2() { return x*x + y*y; }
    inline double magnitude() { return sqrtf(magnitude2()); }

    inline bool operator==(const double2 &rhs) {
            return (x == rhs.x && y == rhs.y);
    }

    inline bool operator!=(const double2 &rhs) {
            return (x != rhs.x || y != rhs.y);
    }

    union {
        struct {
            double x, y;
        };
        struct {
            double s, t;
        };
        double data[2];
    };
};


inline double2 operator*(const double scale, const double2 &rhs) {
    return double2(rhs.x * scale, rhs.y * scale);
}
inline double2 operator-(const double2 &rhs) {
    return double2(-rhs.x, -rhs.y);
}

inline std::ostream& operator<<(std::ostream& os, const double2& f) {
        os <<"[";
        for (unsigned i = 0; i < 2; ++i) {
            os << f.data[i] << ",";
        }
        os << "]";
        return os;
}


// Dense matrix ops
// Naming convention first two letters are types
// ex: mm is matrix matrix, mv is matrix vector
// third letter after _ is type.  Final letters are op.
// memory must be alocated to the output pointer by the caller

inline void vv_dmult(double *a, double *b, double &out, int m) {
    out = 0.f;
    for(int i=0; i<m; i++) {
        out +=a[i]*b[i];
    }
}
#include <assert.h>
 inline void mm_dmult(double *A, double *B, double *out, int rowsA, int colsA, int rowsB, int colsB) {
    assert(colsA == rowsB);
    for(int i = 0; i < rowsA; i++)
    {
        for(int j = 0; j < colsB; j++)
        {
            out[i*colsB+j] = 0;
            for(int k = 0; k < colsA; k++)
            {
                out[i*colsB+j] += A[i*colsA + k] * B[k*colsB + j];
            }
        }
    }
    /*for(int i=0, s=0; i<m; i++) {
        for(int j=0; j<n; j++, s++) {
            out[s] = 0.f;
            for(int k=0; k<n; k++) {
                out[s] += A[(i*n)+k] * B[(k*n)+j];
            }
        }
    }*/
}


inline void mv_dmult(double *A, double *b, double *out, int m, int n) {
    for(int i=0; i<m; i++) {
        out[i] = 0.0;
        for(int j=0; j<n; j++){
            out[i] += A[i*n+j] * b[j];
        }
    }
}

inline void m_dtrans(double *A, double *At, int m, int n) {
    for(int i=0, s=0; i<n; i++) {
        for(int j=0; j<m; j++,s++) {
            At[s] = A[(j*n)+i];
        }
    }
}
