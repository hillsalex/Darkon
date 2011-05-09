#ifndef BYSANDBOX_H
#define BYSANDBOX_H

#include <shot.h>
#include "Sphere.h"
#include <math/CS123Algebra.h>
#include <src/SparseMatrix.h>
class TAMstroke;



class BYSandbox
{
public:
    BYSandbox(DrawEngine* parent);

    ~BYSandbox();

    void orientTexture();
    double2 baryCoord(double2 A, double2 B, double2 C, double2 T);
    double2 relativeCoord(double2 A, double2 B, double2 C);
    //void printMatrix(SparseMatrix* mat) {


    QVector<GLMtriangle*> *m_patch;
    QVector<double2*> *m_verts, *m_tvecs, *m_svecs;
};

#endif // BYSANDBOX_H
