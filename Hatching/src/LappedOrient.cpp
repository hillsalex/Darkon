#include "LappedOrient.h"
#include "LappedUtils.h"
#include "drawengine.h"
#include <QVector>
#include "glm.h"
#include "SparseMatrix.h"
#include "common.h"
#include "src/alglib/solvers.h"
#include "src/alglib/ap.h"
#include "src/alglib/alglibmisc.h"
#include "src/alglib/alglibinternal.h"

using namespace alglib;
using namespace std;

LappedOrient::LappedOrient()
{
//    QImage *original2D = new QImage(1000, 1000, QImage::Format_ARGB32);
//    original2D->fill(0xffffffff);
//    QPainter painter_black(original2D);
//    painter_black.setPen(QPen(QColor(0,0,0,255)));

//    for (int i = 0; i < m_patch->size(); i++) {
//        GLMtriangle* currTri = m_patch->at(i);

//        QPoint vert1(m_verts.at(currTri->vindices[0])->x*100+500, -m_verts.at(currTri->vindices[0])->y*100+500);
//        QPoint vert2(m_verts.at(currTri->vindices[1])->x*100+500, -m_verts.at(currTri->vindices[1])->y*100+500);
//        QPoint vert3(m_verts.at(currTri->vindices[2])->x*100+500, -m_verts.at(currTri->vindices[2])->y*100+500);

//        painter_black.drawLine(vert1, vert2);
//        painter_black.drawLine(vert1, vert3);
//        painter_black.drawLine(vert2, vert3);

//        double avgx = (m_verts.at(currTri->vindices[0])->x + m_verts.at(currTri->vindices[1])->x
//                       + m_verts.at(currTri->vindices[2])->x)/3;
//        double avgy = (m_verts.at(currTri->vindices[0])->y + m_verts.at(currTri->vindices[1])->y
//                       + m_verts.at(currTri->vindices[2])->y)/3;
//        double2 avg = double2(avgx, avgy);
//        QPoint zero(avg.x*100+500, -avg.y*100+500);

//        QPoint t((avg.x+m_tvecs->at(i)->x)*100+500, -(avg.y+m_tvecs->at(i)->y)*100+500);
//        QPoint s((avg.x+m_svecs->at(i)->x)*100+500, -(avg.y+m_svecs->at(i)->y)*100+500);
//        painter_black.drawLine(zero, t);
//        painter_black.drawLine(zero, s);
//    }

//    original2D->save("before.png");
//    orientTexture();
}

LappedOrient::~LappedOrient() {}

void LappedOrient::orientTexture(LappedPatch* patch) {
    //insert each vertex number into a list of vertices.
    QVector<PatchVert*> verts;
    for (int i = 0; i < patch->tris->size(); i++) {
        PatchTri* currtri = patch->tris->at(i);
        if (!verts.contains(currtri->v0)) verts.append(currtri->v0);
        if (!verts.contains(currtri->v1)) verts.append(currtri->v1);
        if (!verts.contains(currtri->v2)) verts.append(currtri->v2);
    }

    real_2d_array H;
    H.setlength(2*verts.size(), 2*verts.size());
    for (int i = 0; i < 2*verts.size(); i++) {
        for (int j = 0; j < 2*verts.size(); j++) {
            H[i][j] = 0;
        }
    }

    double at, bt, as, bs;

    //calculating f0.
    double f0[verts.size()*2];
    memset(f0, 0, sizeof(f0));

    //entering in values for H00 for all vertices.
    for (int tri = 0; tri < patch->tris->size(); tri++) {
        PatchTri* ctri = patch->tris->at(tri);

        //idx locations in H
        int idx_ax = verts.indexOf(ctri->v0)*2;
        int idx_ay = verts.indexOf(ctri->v0)*2+1;
        int idx_bx = verts.indexOf(ctri->v1)*2;
        int idx_by = verts.indexOf(ctri->v1)*2+1;
        int idx_cx = verts.indexOf(ctri->v2)*2;
        int idx_cy = verts.indexOf(ctri->v2)*2+1;

        double2 v0pos = double2(patch->uvs->value(ctri->v0).x, patch->uvs->value(ctri->v0).y);
        double2 v1pos = double2(patch->uvs->value(ctri->v1).x, patch->uvs->value(ctri->v1).y);
        double2 v2pos = double2(patch->uvs->value(ctri->v2).x, patch->uvs->value(ctri->v2).y);
        double2 t = double2(ctri->tangent.x,  ctri->tangent.y);
        double2 s = double2(ctri->tangent.y, -ctri->tangent.x);

        //CONVERT METHOD TO 3D to get alpha and beta
        double2 coeff_t = baryCoord(v0pos, v1pos, v2pos, t); // NEED TO FIX THISSSS
        double2 coeff_s = baryCoord(v0pos, v1pos, v2pos, s); // NEED TO FIX THISSSS
        at = coeff_t.x;
        as = coeff_s.x;
        bt = coeff_t.y;
        bs = coeff_s.y;

        //ax ax = at^2 + as^2
        H[idx_ax][idx_ax] += at*at + as*as;
        //ax bx = at*bt + as*bs
        H[idx_ax][idx_bx] += at*bt + as*bs;
        //ax cx = (- at^2 - at*bt) + (- as^2 - as_t*bs)
        H[idx_ax][idx_cx] += -at*at - at*bt - as*as - as*bs;
        //ay ay = at^2
        H[idx_ay][idx_ay] += at*at + as*as;
        //ay by = at*beta
        H[idx_ay][idx_by] += at*bt + as*bs;
        //ay cy = - at^2 - at*beta
        H[idx_ay][idx_cy] += -at*at - at*bt - as*as - as*bs;
        //bx ax = at*beta
        H[idx_bx][idx_ax] += at*bt + as*bs;
        //bx bx = 0.707/sqrt(2*0beta^2
        H[idx_bx][idx_bx] += bt*bt + bs*bs;
        //bx cx = -at*beta - beta^2
        H[idx_bx][idx_cx] += -bt*bt - at*bt - bs*bs - as*bs;
        //by ay = at*beta
        H[idx_by][idx_ay] += at*bt + as*bs;
        //by by = beta^2
        H[idx_by][idx_by] += bt*bt + bs*bs;
        //by cy = -at*beta - beta^2
        H[idx_by][idx_cy] += - bt*bt - at*bt - bs*bs - as*bs;
        //cx ax = - at^2 - at*beta
        H[idx_cx][idx_ax] += - at*at - at*bt - as*as - as*bs;
        //cx bx = -at*beta - beta^2
        H[idx_cx][idx_bx] += -bt*bt - at*bt - bs*bs - as*bs;
        //cx cx = (at+beta)^2
        H[idx_cx][idx_cx] += at*at + 2*(at*bt) + bt*bt + as*as + 2*(as*bs) + bs*bs;
        //cy ay = - at^2 - at*beta
        H[idx_cy][idx_ay] += - at*at - at*bt - as*as - as*bs;
        //cy by = -at*beta - beta^2
        H[idx_cy][idx_by] += - bt*bt - at*bt - bs*bs - as*bs;
        //cy cy = (at+beta)^2
        H[idx_cy][idx_cy] += at*at + 2*at*bt + bt*bt + as*as + 2*as*bs + bs*bs;

        //filling in f0.
        f0[idx_ax] += -2*as; //-2at*s_x
        f0[idx_ay] += -2*at; //-2at*t_y
        f0[idx_bx] += -2*bs;  //-2beta*s_x
        f0[idx_by] += -2*bt;  //-2beta*t_y
        f0[idx_cx] += 2*(as + bs); //2*at*s_x+2*beta*s_x
        f0[idx_cy] += 2*(at + bt); //2*at*t_y+2*beta*t_y
    }

    //making H01, H10
    SparseMatrix* H10 = new SparseMatrix(2,verts.size()*2);
    SparseMatrix* H01 = new SparseMatrix(verts.size()*2,2);

    //adding seed barycenter = 1/3(phi(A) + phi(B) + phi(C))
    H10->addValue(0, verts.indexOf(patch->seed->v0)*2,   1.0/6.0);
    H10->addValue(1, verts.indexOf(patch->seed->v0)*2+1, 1.0/6.0);
    H10->addValue(0, verts.indexOf(patch->seed->v1)*2,   1.0/6.0);
    H10->addValue(1, verts.indexOf(patch->seed->v1)*2+1, 1.0/6.0);
    H10->addValue(0, verts.indexOf(patch->seed->v2)*2,   1.0/6.0);
    H10->addValue(1, verts.indexOf(patch->seed->v2)*2+1, 1.0/6.0);

    H01->addValue(verts.indexOf(patch->seed->v0)*2,   0, 1.0/6.0);
    H01->addValue(verts.indexOf(patch->seed->v0)*2+1, 1, 1.0/6.0);
    H01->addValue(verts.indexOf(patch->seed->v1)*2,   0, 1.0/6.0);
    H01->addValue(verts.indexOf(patch->seed->v1)*2+1, 1, 1.0/6.0);
    H01->addValue(verts.indexOf(patch->seed->v2)*2,   0, 1.0/6.0);
    H01->addValue(verts.indexOf(patch->seed->v2)*2+1, 1, 1.0/6.0);

    //D
    SparseMatrix *D = new SparseMatrix(H01->operator +(H10->getTranspose()));

    //calculating D*q
    double DQ[verts.size()*2];
    memset(DQ, 0, sizeof(DQ));
    double q[2];
    memset(q, 0, sizeof(q));
    q[0] = (patch->uvs->value(patch->seed->v0).x + patch->uvs->value(patch->seed->v1).x + patch->uvs->value(patch->seed->v2).x)/3.0;
    q[1] = (patch->uvs->value(patch->seed->v0).y + patch->uvs->value(patch->seed->v1).y + patch->uvs->value(patch->seed->v2).y)/3.0;
    D->multiply(DQ, q, 2, 1);

    //calculate -D*q-f0
    real_1d_array Dqf;
    Dqf.setlength(verts.size()*2);
    for (int row = 0; row < D->getM(); row++) {
        Dqf[row] = -DQ[row] - f0[row];
    }

    //GET THE RMATRIXSVD breakdown so that you can do the singular value decomposition
    //reverse the sigma matrix to get the pseudoinverse and then
    //multiply by -dqf to get X.

    real_2d_array u, vt, sigma, Hpseudo, Hpseudo_temp;
    u.setlength(2*verts.size(), 2*verts.size());
    vt.setlength(2*verts.size(), 2*verts.size());
    sigma.setlength(2*verts.size(), 2*verts.size());
    Hpseudo.setlength(verts.size()*2, verts.size()*2);
    Hpseudo_temp.setlength(verts.size()*2, verts.size()*2);

    real_1d_array w, unconstrained;
    w.setlength(2*verts.size());
    unconstrained.setlength(2*verts.size());

    for (int i = 0; i < 2*verts.size(); i++) {
        for (int j = 0; j < 2*verts.size(); j++) {
            u[i][j] = 0;
            vt[i][j] = 0;
            sigma[i][j] = 0;
            Hpseudo_temp[i][j] = 0;
            Hpseudo[i][j] = 0;
        }
        w[i] = 0;
        unconstrained[i] = 0;
    }

    rmatrixsvd(H, 2*verts.size(), 2*verts.size(), 2, 2, 0, w, u, vt);

    //fill in sigma
    for (int i = 0; i < 2*verts.size(); i++) {
        if (w[i] > 0.000000001) sigma[i][i] = 1.0/w[i];
    }

    //Hpseudo = u*sigma*vt;
    rmatrixgemm(verts.size()*2, verts.size()*2, verts.size()*2, 1, u, 0, 0, 0, sigma, 0, 0, 0, 0, Hpseudo_temp, 0, 0);
    rmatrixgemm(verts.size()*2, verts.size()*2, verts.size()*2, 1, Hpseudo_temp, 0, 0, 0, vt, 0, 0, 0, 0, Hpseudo, 0, 0);

    //unconstrained = Hpseudo*Dqf
    rmatrixmv(verts.size()*2, verts.size()*2, Hpseudo, 0, 0, 0, Dqf, 0, unconstrained, 0);

    double scale_oldx = patch->uvs->value(patch->seed->v0).x - patch->uvs->value(patch->seed->v1).x;
    double scale_oldy = patch->uvs->value(patch->seed->v0).y - patch->uvs->value(patch->seed->v1).y;
    double scale_old = sqrt(scale_oldx*scale_oldx + scale_oldy*scale_oldy);

    double scale_newx = unconstrained[verts.indexOf(patch->seed->v0)*2] - unconstrained[verts.indexOf(patch->seed->v1)*2];
    double scale_newy = unconstrained[verts.indexOf(patch->seed->v0)*2+1] - unconstrained[verts.indexOf(patch->seed->v1)*2+1];
    double scale_new = sqrt(scale_newx*scale_newx + scale_newy*scale_newy);

    double scale = scale_old/scale_new;

    for (int i = 0; i < verts.size(); i++) {
        unconstrained[i*2]   *= scale;
        unconstrained[i*2+1] *= scale;
//        cout << "vert " << i << " is (" << unconstrained[i*2] << ", " << unconstrained[i*2+1] << ")" << endl;
    }

    for (int i = 0; i < patch->tris->size(); i++) {
        PatchTri* currtri = patch->tris->at(i);
        patch->uvs->remove(currtri->v0);
        patch->uvs->remove(currtri->v1);
        patch->uvs->remove(currtri->v2);
        float v0x = unconstrained[verts.indexOf(currtri->v0)*2];
        float v0y = unconstrained[verts.indexOf(currtri->v0)*2+1];
        float v1x = unconstrained[verts.indexOf(currtri->v1)*2];
        float v1y = unconstrained[verts.indexOf(currtri->v1)*2+1];
        float v2x = unconstrained[verts.indexOf(currtri->v2)*2];
        float v2y = unconstrained[verts.indexOf(currtri->v2)*2+1];
        patch->uvs->insert(currtri->v0, vec2<float>(v0x, v0y));
        patch->uvs->insert(currtri->v1, vec2<float>(v1x, v1y));
        patch->uvs->insert(currtri->v2, vec2<float>(v2x, v2y));
    }
//    QImage *new2D = new QImage(1000, 1000, QImage::Format_ARGB32);
//    new2D->fill(0xffffffff);
//    QPainter painter_black(new2D);
//    painter_black.setPen(QPen(QColor(0,0,0,255)));

//    for (int i = 0; i < m_patch->size(); i++) {
//        GLMtriangle* currTri = m_patch->at(i);

//        QPoint vert1(unconstrained[currTri->vindices[0]*2]*100+500, -unconstrained[currTri->vindices[0]*2+1]*100+500);
//        QPoint vert2(unconstrained[currTri->vindices[1]*2]*100+500, -unconstrained[currTri->vindices[1]*2+1]*100+500);
//        QPoint vert3(unconstrained[currTri->vindices[2]*2]*100+500, -unconstrained[currTri->vindices[2]*2+1]*100+500);

//        painter_black.drawLine(vert1, vert2);
//        painter_black.drawLine(vert1, vert3);
//        painter_black.drawLine(vert2, vert3);
//    }

//    new2D->save("after.png");
}

/* baryCoord: Barycentric Coordinates
 * To express phi(T) as a linear combination of points phi(A), phi(B), phi(C),
 * phi(T) = at * phi(A) + beta * phi(B) + gamma * phi(C)
 * at + beta + gamma = 0;
 * assuming T is coplanar with A, B, C
 * 05.07 : pretty sure this is correct now.
 */
double2 LappedOrient::baryCoord(double2 A, double2 B, double2 C, double2 T) {
    double det   = 1/((A.x-C.x)*(B.y-C.y)-(A.y-C.y)*(B.x-C.x));
    double alpha = ((B.y-C.y)*T.x + (-B.x+C.x)*T.y)*det;
    double beta  = ((-A.y+C.y)*T.x + (A.x-C.x)*T.y)*det;
    return double2(alpha, beta);
}
