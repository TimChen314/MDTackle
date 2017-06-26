#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <math.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/time.h>
#include <vector>
using namespace std;
#ifndef __DATA_STRUCT__
#define __DATA_STRUCT__
//each of "bond" types, "angle" types and "vec" types has two different kind

template <typename T>
inline vector<T>& VecResize(vector<T>& vec, int np)
{
    if (vec.size() != 0 && np >= 0)
        vec.resize(np);
    return vec;
}

struct Bond {
    Bond(){};
    Bond(string bond_type, unsigned int tag_a, unsigned int tag_b)
        : type(bond_type)
        , a(tag_a)
        , b(tag_b)
    {
    }
    string type;
    unsigned int a;
    unsigned int b;

    friend bool operator==(const Bond& lift, const Bond& right);
    friend bool operator<(const Bond& lift, const Bond& right);

private:
    friend ostream& operator<<(ostream& os, const Bond& bond);
};
inline bool operator==(const Bond& lift, const Bond& right)
{
    if (lift.type == right.type && lift.a == right.a && lift.b == right.b)
        return true;
    return false;
}
inline bool operator<(const Bond& lift, const Bond& right)
{
    if (lift.a < right.a)
        return true;
    else if (lift.a == right.a && lift.b < right.b)
        return true;
    return false;
}

struct Angle {
    Angle(){};
    Angle(string angle_type, unsigned int tag_a, unsigned int tag_b, unsigned int tag_c)
        : type(angle_type)
        , a(tag_a)
        , b(tag_b)
        , c(tag_c)
    {
    }
    string type;    //!< The type index of the angle
    unsigned int a; //!< The tag of the first particle in the angle
    unsigned int b; //!< The tag of the second particle in the angle
    unsigned int c; //!< The tag of the third particle in the angle
    inline friend bool operator==(const Angle& lift, const Angle& right)
    {
        if (lift.type == right.type && lift.a == right.a && lift.b == right.b && lift.c == right.c)
            return true;
        return false;
    }
    inline friend bool operator<(const Angle& lift, const Angle& right)
    {
        if (lift.a < right.a)
            return true;
        else if (lift.a == right.a && lift.b < right.b)
            return true;
        else if (lift.a == right.a && lift.b == right.b && lift.c < right.c)
            return true;
        return false;
    }

private:
    friend ostream& operator<<(ostream& os, const Angle& angle);
};

struct Dihedral {

    Dihedral(unsigned int dihedral_type, unsigned int tag_a, unsigned int tag_b, unsigned int tag_c, unsigned int tag_d)
        : type(dihedral_type)
        , a(tag_a)
        , b(tag_b)
        , c(tag_c)
        , d(tag_d)
    {
    }
    unsigned int type; //!< The type index of the dihedral
    unsigned int a;    //!< The tag of the first particle in the dihedral
    unsigned int b;    //!< The tag of the second particle in the dihedral
    unsigned int c;    //!< The tag of the third particle in the dihedral
    unsigned int d;    //!< The tag of the forth particle in the dihedral
};

struct IntegratorVariables {
    std::string type;             //!<The type of integrator (NVT, NPT, etc.)
    std::vector<double> variable; //!<Variables that define the integration state
};

struct vec4D {

    vec4D()
        : x(0.0)
        , y(0.0)
        , z(0.0)
        , w(0.0)
    {
    }

    vec4D(double xp, double yp, double zp, double wp)
        : x(xp)
        , y(yp)
        , z(zp)
        , w(wp)
    {
    }
    double x; //!< x-component
    double y; //!< y-component
    double z; //!< z-component
    double w;
};

struct vec4Dfloat {

    vec4Dfloat()
        : x(0.0)
        , y(0.0)
        , z(0.0)
        , w(0.0)
    {
    }

    vec4Dfloat(float xp, float yp, float zp, float wp)
        : x(xp)
        , y(yp)
        , z(zp)
        , w(wp)
    {
    }
    float x; //!< x-component
    float y; //!< y-component
    float z; //!< z-component
    float w;
};

struct vec {

    vec()
        : x(0.0)
        , y(0.0)
        , z(0.0)
    {
    }

    vec(double xp, double yp, double zp)
        : x(xp)
        , y(yp)
        , z(zp)
    {
    }
    double x; //!< x-component
    double y; //!< y-component
    double z; //!< z-component
    inline vec& operator+=(const vec&);
    inline vec& operator-=(const vec&);
    inline vec& operator*=(const vec&);
    inline vec& operator/=(const vec&);
    inline vec& operator=(const vec&);

    // multiplication of vector and normal substraction
    inline vec& operator/=(const double&);
    inline vec& operator-=(const double&);

private:
    friend ostream& operator<<(ostream& os, const vec& v);
};
inline void pbc(const vec&, const vec&, vec&);

//! simple integer vec for storing particle data
struct vec_int {
    //! Default construtor
    vec_int()
        : x(0)
        , y(0)
        , z(0)
    {
    }
    //! Constructs a vec with given components
    /*! \param xp x-component
      \param yp y-component
      \param zp z-component
      */
    vec_int(int xp, int yp, int zp)
        : x(xp)
        , y(yp)
        , z(zp)
    {
    }
    int x; //!< x-component
    int y; //!< y-component
    int z; //!< z-component
};

const unsigned int NO_BODY = 0xffffffff;

// Function
// Calculate
void FoldBackToBox(vector<vec>& pos, vec& box);
void QuestRealCoorCMAndSetToOri(vector<vec>& pos, vector<vec>& pos_cm, const unsigned nm_tot, const vector<unsigned int>& MolNA);
// Convert variant type
void Vector_vec_To_Vector_float(const vector<vec>& vv, vector<float>& vx, vector<float>& vy, vector<float>& vz);
// Print Process Info
void ProcessBar(unsigned percent);
long unsigned getAutoCorrelationProcessPercent(long unsigned ifram, long unsigned nfram_remain);
void printOnePercentTime(timeval& start);
// Dump
void DumpAs3202Version(const string fname, const vector<vec>& pos, const vector<vec>& vel, const vec& box, const unsigned timestep, const bool m_output_velocity);

inline void pbc(const vec& box, const vec& boxINV, vec& p)
{
    p.x -= box.x * rintf(p.x * boxINV.x);
    p.y -= box.y * rintf(p.y * boxINV.y);
    p.z -= box.z * rintf(p.z * boxINV.z);
}
// unitary operator of vec
inline vec& vec::operator+=(const vec& v1)
{
    x += v1.x;
    y += v1.y;
    z += v1.z;
    return *this;
}
inline vec& vec::operator-=(const vec& dv1)
{
    x -= dv1.x;
    y -= dv1.y;
    z -= dv1.z;
    return *this;
}
inline vec& vec::operator*=(const vec& dv1)
{
    x *= dv1.x;
    y *= dv1.y;
    z *= dv1.z;
    return *this;
}
inline vec& vec::operator/=(const vec& dv1)
{
    x /= dv1.x;
    y /= dv1.y;
    z /= dv1.z;
    return *this;
}
inline vec& vec::operator=(const vec& dv1)
{
    x = dv1.x;
    y = dv1.y;
    z = dv1.z;
    return *this;
}

inline vec& vec::operator/=(const double& d1)
{
    x /= d1;
    y /= d1;
    z /= d1;
    return *this;
}
inline vec& vec::operator-=(const double& d1)
{
    x -= d1;
    y -= d1;
    z -= d1;
    return *this;
}

// binary operator of vec
inline vec operator+(const vec& dv1, const vec& dv2)
{
    vec dv(dv1.x + dv2.x, dv1.y + dv2.y, dv1.z + dv2.z);
    return dv;
}
inline vec operator-(const vec& dv1, const vec& dv2)
{
    vec dv(dv1.x - dv2.x, dv1.y - dv2.y, dv1.z - dv2.z);
    return dv;
}
inline vec operator*(const vec& dv1, const vec& dv2)
{
    vec dv(dv1.x * dv2.x, dv1.y * dv2.y, dv1.z * dv2.z);
    return dv;
}
inline vec operator*(const vec& dv1, const vec_int& vi2)
{
    vec dv(dv1.x * vi2.x, dv1.y * vi2.y, dv1.z * vi2.z);
    return dv;
}
inline vec operator/(const vec& dv1, const vec& dv2)
{
    vec dv(dv1.x / dv2.x, dv1.y / dv2.y, dv1.z / dv2.z);
    return dv;
}
inline double ScalarProduct(const vec& dv1, const vec& dv2)
{
    return dv1.x * dv2.x + dv1.y * dv2.y + dv1.z * dv2.z;
}

// reload  <<
inline ostream& operator<<(ostream& os, const Bond& bond)
{
    os << bond.type << " " << bond.a << " " << bond.b;
    return os;
}

inline ostream& operator<<(ostream& os, const Angle& angle)
{
    os << angle.type << " " << angle.a << " " << angle.b << " " << angle.c;
    return os;
}

inline ostream& operator<<(ostream& os, const vec& v)
{
    os << v.x << " " << v.y << " " << v.z;
    return os;
}

inline ostream& operator<<(ostream& os, const vec_int& v)
{
    os << v.x << " " << v.y << " " << v.z;
    return os;
}

#endif
