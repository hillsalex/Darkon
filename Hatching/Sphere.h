#ifndef SPHERE_H
#define SPHERE_H
#include "math/CS123Algebra.h"
#include <qgl.h>
class Sphere
{
public:
    Sphere(int stacks, int slices);
    void makeGeometry(int stacks,int slices);
    Vector4* verts;
    Vector4* norms;
    Vector4* tans;
    vec2<REAL>* texs;
    int numVerts;
    void draw();
};

#endif // SPHERE_H
