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

        if (fabs(normal.x)>0.01)
        {
            vert = vec3<float>(-(normal.y+normal.z)/normal.x,normal.x,normal.x);
            vert=vert/normal.x;
        }
        else
        {
            if (fabs(normal.y)>0.01)
            {
                vert= vec3<float>(normal.y,-(normal.x+normal.z),normal.y);
                vert=vert/normal.y;
            }

            else
                if (fabs(normal.z) > 0.01)
                {
                vert= vec3<float>(normal.z,normal.z,-(normal.y+normal.x));
                vert=vert/normal.z;
            }
        }
        /*
        vert.x = model->vertices[vertexIndex*3];
        vert.y = model->vertices[vertexIndex*3+1];
        vert.z = model->vertices[vertexIndex*3+2];
*/
        QList<int> adjacentVertices = adjacencyMatrix.value(vertexIndex)->values(); //Adjacent vertices
        //int v0index = adjacentVertices[0]; //first adjacent vertex, for calculating basis

        //vec3<float> v0;
        //v0.x = model->vertices[v0index*3];
        //v0.y = model->vertices[v0index*3+1];
        //v0.z = model->vertices[v0index*3+2];
        //Translate v0 to vert as origin, project to uv plane
        normal.normalize();
        //vec3<float> u = v0-vert;
        vec3<float> u = vert;
        u.normalize();
        u = u - normal * u.dot(normal);
        u.normalize();
        vert.normalize();
        vec3<float> v = normal.cross(u);
        v.normalize();

        //cout << "U/V/N/CURVE" << endl;
        //printVec3(&u);
        //printVec3(&v);
        //printVec3(&normal);
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
        vec3<float> curvature1,curvature2,curvature,minCurvature;
        curvature1.x = evecs[0][0];
        curvature1.y = evecs[1][0];
        curvature1.z = 0;

        curvature2.x = evecs[0][1];
        curvature2.y = evecs[1][1];
        curvature2.z = 0;

        SparseMatrix* curvatureMat = vec3toMat(curvature1);
        SparseMatrix Mtrans = M.getTranspose();
        curvatureMat = &(Mtrans * (*curvatureMat));
        curvature1 = *mattoVec3(curvatureMat);
        curvature1.normalize();

        curvatureMat = vec3toMat(curvature2);
        curvatureMat = &(Mtrans * (*curvatureMat));
        curvature2 = *mattoVec3(curvatureMat);
        curvature2.normalize();
        //cout << evals[0] << "," << evals[1] << ',' << (fabs(evals[0])-fabs(evals[1])) << endl;
        if (fabs(evals[0])<fabs(evals[1]))

        {
            curvature=curvature2;
            minCurvature=curvature1;
        }

        else
        {
            if (fabs(curvature2.dot(curvature1))>0.9)
            {
                if (curvature1.x > curvature2.x)
                {

                    curvature = curvature1;
                    minCurvature=curvature2;
                }
                else
                {
                    if (curvature1.y > curvature2.y)
                    {

                        curvature = curvature1;
                        minCurvature=curvature2;
                    }
                    else
                    {
                        if (curvature1.z > curvature2.z)
                        {

                            curvature = curvature1;
                            minCurvature=curvature2;
                        }
                    }
                }
            }
        }
        curvature.normalize();
        curvature = curvature - normal * curvature.dot(normal);
        curvature.normalize();

        minCurvature.normalize();
        minCurvature = minCurvature - normal * minCurvature.dot(normal);
        minCurvature.normalize();

        //printVec3(&curvature);

        model->vertMaxCurvatures[vertexIndex*3]=curvature.x;
        model->vertMaxCurvatures[vertexIndex*3+1]=curvature.y;
        model->vertMaxCurvatures[vertexIndex*3+2]=curvature.z;

        model->vertMinCurvatures[vertexIndex*3]=minCurvature.x;
        model->vertMinCurvatures[vertexIndex*3+1]=minCurvature.y;
        model->vertMinCurvatures[vertexIndex*3+2]=minCurvature.z;
    }

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
        vec3<float> eprime = vec3<float>(0,0,0);
        for (int adjVert=0;adjVert<adjacentVertices.size();adjVert++)
        {
            int adjIndex = adjacentVertices[adjVert];
            vec3<float> alpha;
            alpha.x = model->vertMaxCurvatures[adjIndex*3];
            alpha.y = model->vertMaxCurvatures[adjIndex*3+1];
            alpha.z = model->vertMaxCurvatures[adjIndex*3+2];


            eprime.x += model->vertMaxCurvatures[adjIndex*3];
            eprime.y += model->vertMaxCurvatures[adjIndex*3+1];
            eprime.z += model->vertMaxCurvatures[adjIndex*3+2];
        }
        eprime=eprime/((float)adjacentVertices.size());


        eprime.normalize();
        eprime = eprime-(eprime.dot(normal))*normal;
        eprime.normalize();
        model->vertCurvatures[i*3]=eprime.x;
        model->vertCurvatures[i*3+1]=eprime.y;
        model->vertCurvatures[i*3+2]=eprime.z;

    }



    GLfloat* vertices = model->vertices;
    GLMtriangle* triangles = model->triangles;
    for (int i=0;i<model->numtriangles;i++)
    {
        GLMtriangle tri = model->triangles[i];






        vec3<float> v1,v2,v3;
        v1.x+= vertices[triangles[i].vindices[0]*3];
        v1.y+= vertices[triangles[i].vindices[0]*3+1];
        v1.z+= vertices[triangles[i].vindices[0]*3+2];

        v2.x+= vertices[triangles[i].vindices[1]*3];
        v2.y+= vertices[triangles[i].vindices[1]*3+1];
        v2.z+= vertices[triangles[i].vindices[1]*3+2];

        v3.x+=vertices[triangles[i].vindices[2]*3];
        v3.y+= vertices[triangles[i].vindices[2]*3+1];
        v3.z+= vertices[triangles[i].vindices[2]*3+2];

        vec3<float> normal = (v2-v1).cross(v3-v1);
        normal.normalize();

        float cx=0;
        float cy=0;
        float cz=0;
        for (int vi =0;vi<3;vi++)
        {
            int vIndex = tri.vindices[vi];
            cx += model->vertCurvatures[vIndex*3];
            cy += model->vertCurvatures[vIndex*3+1];
            cz += model->vertCurvatures[vIndex*3+2];
        }
        cx=cx/3.0;
        cy=cy/3.0;
        cz=cz/3.0;
        vec3<float> c = vec3<float>(cx,cy,cz);
        c.normalize();
        c = c-(c.dot(normal))*normal;
        c.normalize();
        model->triCurvatures[i*3]=c.x;
        model->triCurvatures[i*3+1]=c.y;
        model->triCurvatures[i*3+2]=c.z;

    }
}


