#ifndef LAPPEDORIENT_H
#define LAPPEDORIENT_H

#include <shot.h>
#include "Sphere.h"
#include "LappedUtils.h"
#include <math/CS123Algebra.h>
#include <src/SparseMatrix.h>
class TAMstroke;



class LappedOrient
{
public:
    LappedOrient();

    ~LappedOrient();

    void orientTexture(LappedPatch patch);
    double2 baryCoord(double2 A, double2 B, double2 C, double2 T);

    QVector<GLMtriangle*> *m_patch;
    QVector<double2*> *m_verts, *m_tvecs, *m_svecs;
};

#endif // LAPPEDORIENT_H
