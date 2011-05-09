#include "MPSandbox.h"
#include "drawengine.h"
#include <QGLShaderProgram>
#include <QGLShader>
#include <QHash>
#include <QQueue>
#include <QPainter>
#include "tamgenerator.h"
#include "Sphere.h"
#include "assert.h"
#include <LappedUtils.h>
#include "meshoperator.h"
using std::cout;
using std::endl;
using std::cin;
using std::vector;
using std::stringstream;
MPSandbox::MPSandbox(DrawEngine* parent,QHash<QString, QGLShaderProgram *>* shad, QHash<QString, GLuint>* tex, QHash<QString, Model>* mod) : Shot(parent,shad,tex,mod)
{

}

MPSandbox::~MPSandbox()
{

}


//In begin, initialize things that could not have been initialized beforehand
//(gl state)
void MPSandbox::begin()
{
    QImage img;
    img.load("/home/mprice/Desktop/Patch/PatchMask.png");
    LappedUtils* lu = new LappedUtils();
    MeshOperator* mo = new MeshOperator();
    GLMmodel* mod = models_->value("cube").model;
    mo->calculateCurvatures(mod);

    polyHull* pHull = lu->getPolyHull(&img,6);

    //DEBUG DRAW OUTPUT
    /*QPainter patr(&img);
    patr.setPen(Qt::black);
    QList<vert2d*>* vLi = pHull->verts;
    for(int i=0;i<vLi->size();i++)
    {
        vert2d* v = vLi->at(i);
        patr.drawLine(v->e1->v1->x,v->e1->v1->y,v->e1->v2->x,v->e1->v2->y);
        patr.drawLine(v->e2->v1->x,v->e2->v1->y,v->e2->v2->x,v->e2->v2->y);
    }
    QList<edge2d*>* eLi = pHull->edges;
    for(int i=0;i<eLi->size();i++)
    {
        edge2d* e = eLi->at(i);
        patr.drawLine(e->v1->x, e->v1->y, e->v2->x, e->v2->y);
    }
    patr.setPen(Qt::blue);
    for(int i=0;i<vLi->size();i++)
    {
        vert2d* v = vLi->at(i);
        assert(v->nedges==2);
        patr.drawPoint(v->x,v->y);
    }


    for(int i=0; i<50; i++)
    {
    float s0,t0,s1,t1;
    s0 = rand()%100/100.0;
    t0 = rand()%100/100.0;
    s1 = rand()%100/100.0;
    t1 = rand()%100/100.0;


    //(0,107.52) ---- (32,108.8)
    //s0 = 0.0;
    //t0 = 1.0- (107.52 / img.height());
    //s1 = 32.0/img.width();
    //t1 = 1.0-(108.8 / img.height());


    cout<<"*******************************************************"<<endl;
    cout<<"("<<s0*img.width() << "," << (1.0-t0)*img.height() <<") ---- (" << s1*img.width() <<","<< (1.0-t1)*img.height() <<")"<<endl;
    if(pHull->isectHullUV(s0,t0,s1,t1))
    {patr.setPen(Qt::red); cout<<"red"<<endl;}
    else
    {patr.setPen(Qt::green); cout<<"green"<<endl;}
    cout<<"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"<<endl;

    patr.drawLine(img.width()*s0, (1.0-t0)*img.height(), img.width()*s1, (1.0-t1)*img.height());
    }*/


    QList<LappedPatch*>* LP = lu->generatePatches(mod,pHull);
    cout<<"got patches.  patches: "<<LP->size()<<endl;
    cout<<"patch 1 tris: "<<LP->at(0)->tris->size()<<endl;
    lu->vizualizePatch(LP->at(0),&img);
    img.save("/home/mprice/Desktop/Patch/Collision.png","PNG");

}



//called every frame before draw.
void MPSandbox::update()
{
    m_framesElapsed++;
}

void MPSandbox::draw()
{
    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
    glEnable(GL_TEXTURE_2D);
}
