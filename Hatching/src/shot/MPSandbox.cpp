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
    lu = new LappedUtils();
    MeshOperator* mo = new MeshOperator();
    mod = models_->value("evenodder").model;
    mo->calculateCurvatures(mod);

    polyHull* pHull = lu->getPolyHull(&img,10);

    //DEBUG DRAW OUTPUT
    /*
    QPainter patr(&img);
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


    //cout<<"*******************************************************"<<endl;
   // cout<<"("<<s0*img.width() << "," << (1.0-t0)*img.height() <<") ---- (" << s1*img.width() <<","<< (1.0-t1)*img.height() <<")"<<endl;
    if(pHull->isectHullUV(s0,t0,s1,t1))
    {patr.setPen(Qt::red); cout<<"red"<<endl;}
    else
    {patr.setPen(Qt::green); cout<<"green"<<endl;}
    //cout<<"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"<<endl;

    patr.drawLine(img.width()*s0, (1.0-t0)*img.height(), img.width()*s1, (1.0-t1)*img.height());
    }
img.save("/home/mprice/Desktop/Patch/Collision.png","PNG");*/

    LP = lu->generatePatches(mod,pHull);
    //cout<<"got patches.  patches: "<<LP->size()<<endl;
    //cout<<"patch 1 tris: "<<LP->at(0)->tris->size()<<endl;
    //lu->printPatchTri2d(LP->at(0)->tris->at(0));

    /*
    lu->vizualizePatch(LP->at(0),&img);
    img.save("realPatch.png","PNG");
    img.load("PatchMask.png");
    lu->vizualizePatch(LP->at(1),&img);
    img.save("realPatch2.png","PNG");
    img.load("PatchMask.png");
    lu->vizualizePatch(LP->at(2),&img);
    img.save("realPatch3.png","PNG");
    img.load("PatchMask.png");
    lu->vizualizePatch(LP->at(3),&img);
    img.save("realPatch4.png","PNG");
    img.load("PatchMask.png");
    lu->vizualizePatch(LP->at(4),&img);
    img.save("realPatch5.png","PNG");
    img.load("PatchMask.png");
*/










    TAMGenerator tgen;
    int TAMtones = 6;
    int TAMsizes = 4;
    int TAMmaxw = 256;
    QImage** imgTAM = tgen.ldImgTAM("../Hatching/src/TAM/",TAMtones,TAMsizes);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHT0);
    float lightpos[4];
    lightpos[0]=0.f;
    lightpos[1]=-2.f;
    lightpos[2]=-1.f;
    lightpos[3]=0.f;
    glLightfv(GL_LIGHT0,GL_POSITION,lightpos);

    glMatrixMode(GL_TEXTURE);
    //glScalef(4,4,4);
    glTranslatef(0.5, 0.5, 0);
    glRotatef(90.0, 0, 0, 1);
    glTranslatef(-0.5, -0.5, 0);
    glMatrixMode(GL_MODELVIEW);

    glEnable(GL_TEXTURE_2D);

    GLuint tid[7];
    glGenTextures(7,tid);

    glClearColor(1.0,1.0,1.0,1.0);

    for(int i=0;i<6;i++)
    {
    glActiveTexture(GL_TEXTURE0+i);
    glBindTexture(GL_TEXTURE_2D,tid[i]);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 3);
    glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);


    QImage tex = QGLWidget::convertToGLFormat(*imgTAM[3*TAMtones+i]);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB, TAMmaxw, TAMmaxw, 0, GL_BGRA, GL_UNSIGNED_BYTE, tex.bits());
    QImage tex2 = QGLWidget::convertToGLFormat(*imgTAM[2*TAMtones+i]);
    glTexImage2D(GL_TEXTURE_2D, 1, GL_RGBA32F_ARB, tex2.width(), tex2.width(), 0, GL_BGRA, GL_UNSIGNED_BYTE, tex2.bits());
    QImage tex3 = QGLWidget::convertToGLFormat(*imgTAM[1*TAMtones+i]);
    glTexImage2D(GL_TEXTURE_2D, 2, GL_RGBA32F_ARB, tex3.width(), tex3.width(), 0, GL_BGRA, GL_UNSIGNED_BYTE, tex3.bits());
    QImage tex4 = QGLWidget::convertToGLFormat(*imgTAM[0*TAMtones+i]);
    glTexImage2D(GL_TEXTURE_2D, 3, GL_RGBA32F_ARB, tex4.width(), tex4.width(), 0, GL_BGRA, GL_UNSIGNED_BYTE, tex4.bits());

    }

    //our alpha blob!
    glActiveTexture(GL_TEXTURE0+6);
    glBindTexture(GL_TEXTURE_2D,tid[6]);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 3);
    glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
    QImage atex = QGLWidget::convertToGLFormat(img);
    glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA32F_ARB, img.width(), img.height(), 0, GL_BGRA, GL_UNSIGNED_BYTE, atex.bits());


    //SHADER
    QGLShaderProgram* shd = new QGLShaderProgram();
    shd->addShaderFromSourceFile(QGLShader::Vertex,"../Hatching/src/shaders/hatch.vert");
    shd->addShaderFromSourceFile(QGLShader::Fragment,"../Hatching/src/shaders/hatch.frag");
    shd->link();
    shd->bind();
    shd->setUniformValue("tone0",0);
    shd->setUniformValue("tone1",1);
    shd->setUniformValue("tone2",2);
    shd->setUniformValue("tone3",3);
    shd->setUniformValue("tone4",4);
    shd->setUniformValue("tone5",5);
    shd->setUniformValue("alph",6);
    //shd->release();
}



//called every frame before draw.
void MPSandbox::update()
{
    m_framesElapsed++;
}

void MPSandbox::draw()
{

    float lightpos[4];
    lightpos[0]=1.2f;
    lightpos[1]=-2.f * sin(m_framesElapsed/20.0);
    lightpos[2]=-1.f;
    lightpos[3]=0.f;
    glLightfv(GL_LIGHT0,GL_POSITION,lightpos);
/*
    if(m_framesElapsed%100==0)
    {glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
        if(m_framesElapsed%200==0)
        {
            glEnable(GL_DEPTH_TEST);
        }
        else
        {
            glDisable(GL_DEPTH_TEST);
        }
    }
*/
    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHTING);
    //glEnable(GL_BLEND);
    //glDisable(GL_DEPTH_TEST);
    //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_TEXTURE_2D);
    lu->drawFromPatches(LP, mod);
    //lu->DrawSinglePatch(LP,mod,(m_framesElapsed/20)%LP->size());
    //glTranslatef(0.5,0,0);
    //glCallList(models_->value("teapot").idx);
}
