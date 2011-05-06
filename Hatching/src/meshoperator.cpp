#include "meshoperator.h"
#include "src/alglib/solvers.h"
#include "src/alglib/ap.h"
#include "src/alglib/alglibmisc.h"
#include "src/alglib/alglibinternal.h"
#include "qhash.h"
#include "qset.h"
#include "SparseMatrix.h"
#include "math/CS123Algebra.h"
#include <iostream>
#include "math.h"

using namespace alglib;
using namespace std;

MeshOperator::MeshOperator()
{
}

void printVec3(vec3<float>* v)
{

     cout << '[' << v->x << " , " << v->y << " , " << v->z <<']' << endl;
}

void printMatrix(SparseMatrix* s)
{
    for (int i=0;i<s->getM();i++)
    {
        cout << '[';
        for (int j=0;j<s->getN();j++)
        {
            cout << s->getValue(i,j) << ',';
        }
        cout <<']' << endl;
    }
}


vec3<float>* mattoVec3(SparseMatrix* m)
{
    vec3<float>* v = new vec3<float>();
    v->x= m->getValue(0,0);
    v->y= m->getValue(1,0);
    v->z= m->getValue(2,0);
}


SparseMatrix* vec3toMat(vec3<float> v)
{
    SparseMatrix* m = new SparseMatrix(3,1);
    m->setValue(0,0,v.x);
    m->setValue(1,0,v.y);
    m->setValue(2,0,v.z);
}
void MeshOperator::calculateCurvatures(GLMmodel* model)
{
    QHash<int,QSet<int>* > adjacencyMatrix;
    QHash<int, int> vertexNormalComps;
    for(int i=0;i<model->numtriangles;i++) //Create adjacency matrix and vertex-normal mapping

    {
        GLMtriangle tri = model->triangles[i];
        int v0 = tri.vindices[0];
        int v1 = tri.vindices[1];
        int v2 = tri.vindices[2];

        if (adjacencyMatrix.contains(v0))
        {
            adjacencyMatrix.value(v0)->insert(v1);
            adjacencyMatrix.value(v0)->insert(v2);
        }
        else
        {
            QSet<int>* set = new QSet<int>();
            set->insert(v1);
            set->insert(v2);
            adjacencyMatrix.insert(v0,set);
        }

        if (adjacencyMatrix.contains(v1))
        {
            adjacencyMatrix.value(v1)->insert(v0);
            adjacencyMatrix.value(v1)->insert(v2);
        }
        else
        {
            QSet<int>* set = new QSet<int>;
            set->insert(v0);
            set->insert(v2);
            adjacencyMatrix.insert(v1,set);
        }

        if (adjacencyMatrix.contains(v2))
        {
            adjacencyMatrix.value(v2)->insert(v1);
            adjacencyMatrix.value(v2)->insert(v0);
        }
        else
        {
            QSet<int>* set = new QSet<int>;
            set->insert(v1);
            set->insert(v0);
            adjacencyMatrix.insert(v2,set);
        }
        if (!vertexNormalComps.contains(v0)) vertexNormalComps.insert(v0,tri.nindices[0]);
        if (!vertexNormalComps.contains(v1)) vertexNormalComps.insert(v1,tri.nindices[1]);
        if (!vertexNormalComps.contains(v2)) vertexNormalComps.insert(v2,tri.nindices[2]);
    }
    QList<int> keys = adjacencyMatrix.keys();
    for (int i=0;i<keys.size();i++)
    {

        int vertexIndex = keys[i]; //Index of current vertex
        int normalIndex = vertexNormalComps.value(vertexIndex); //Normal index of current vertex

        vec3<float> normal; //Normal of current index
        vec3<float> vert; //Position of current vertex
        normal.x = model->normals[normalIndex*3];
        normal.y = model->normals[normalIndex*3+1];
        normal.z = model->normals[normalIndex*3+2];
        vert.x = model->vertices[vertexIndex*3];
        vert.y = model->vertices[vertexIndex*3+1];
        vert.z = model->vertices[vertexIndex*3+2];

        QList<int> adjacentVertices = adjacencyMatrix.value(vertexIndex)->values(); //Adjacent vertices
        int v0index = adjacentVertices[0]; //first adjacent vertex, for calculating basis

        vec3<float> v0;
        v0.x = model->vertices[v0index*3];
        v0.y = model->vertices[v0index*3+1];
        v0.z = model->vertices[v0index*3+2];
        //Translate v0 to vert as origin, project to uv plane
        normal.normalize();
        vec3<float> u = v0-vert;
        u.normalize();
        u = u - normal * u.dot(normal);
        u.normalize();

        vec3<float> v = normal.cross(u);
        v.normalize();

        SparseMatrix M = SparseMatrix(3,3);
        M.setValue(0,0,u.x);
        M.setValue(0,1,u.y);
        M.setValue(0,2,u.z);

        M.setValue(1,0,v.x);
        M.setValue(1,1,v.y);
        M.setValue(1,2,v.z);

        M.setValue(2,0,normal.x);
        M.setValue(2,1,normal.y);
        M.setValue(2,2,normal.z);
        real_2d_array Q;
        Q.setlength(adjacentVertices.size(),3);
        real_1d_array B;
        B.setlength(adjacentVertices.size());
        //SparseMatrix *Q = new SparseMatrix(adjacentVertices.size(),3);
        //double *B = new double[adjacentVertices.size()];
        for (int adjIndex=0;adjIndex<adjacentVertices.size();adjIndex++)
        {
            int adjVertIndex = adjacentVertices[adjIndex];
            vec3<float> adjVert;
            adjVert.x = model->vertices[adjVertIndex*3];
            adjVert.y = model->vertices[adjVertIndex*3+1];
            adjVert.z = model->vertices[adjVertIndex*3+2];
            SparseMatrix* wi = vec3toMat(adjVert);
            wi = &(M * (*wi));
            float x = wi->getValue(0,0);
            float y = wi->getValue(1,0);
            float z = wi->getValue(2,0);
            Q[adjIndex][0]=x*x;
            Q[adjIndex][1]=2*x*y;
            Q[adjIndex][2]=y*y;
            B[adjIndex]=z;
        }

        real_1d_array X;
        X.setlength(3);
        densesolverlsreport myRep;
        ae_int_t  myAEint = 1;
        rmatrixsolvels(Q,adjacentVertices.size(),3,B,0.0,myAEint,myRep,X);
        real_2d_array eigMat = real_2d_array();
        eigMat.setlength(2,2);
        eigMat[0][0]=X[0];
        eigMat[0][1]=X[1];
        eigMat[1][0]=X[1];
        eigMat[1][1]=X[2];
        real_1d_array evals;
        real_2d_array evecs;
        smatrixevd(eigMat,2,1,false,evals,evecs);
        vec3<float> curvature;
        int which=0;
        if (fabs(evals[0])<fabs(evals[1]))
            which=1;
        curvature.x = evecs[0][which];
        curvature.y = evecs[1][which];
        curvature.x = 0;
        SparseMatrix* curvatureMat = vec3toMat(curvature);
        SparseMatrix Mtrans = M.getTranspose();
        curvatureMat = &(Mtrans * (*curvatureMat));
        curvature = *mattoVec3(curvatureMat);
        curvature.normalize();
        model->curvatures[vertexIndex*3]=curvature.x;
        model->curvatures[vertexIndex*3+1]=curvature.y;
        model->curvatures[vertexIndex*3+2]=curvature.z;
    }
}


