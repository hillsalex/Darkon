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
using std::cout;
using std::endl;
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
    /*
    polyHull* pHull = lu->getPolyHull(&img,0);

    //DEBUG DRAW OUTPUT
    QPainter patr(&img);
    patr.setPen(Qt::black);
    QList<vert2d*>* vLi = pHull->verts;
    for(int i=0;i<vLi->size();i++)
    {
        vert2d* v = vLi->at(i);
        patr.drawLine(v->e1->v1->x,v->e1->v1->y,v->e1->v2->x,v->e1->v2->y);
        patr.drawLine(v->e2->v1->x,v->e2->v1->y,v->e2->v2->x,v->e2->v2->y);
    }
    patr.setPen(Qt::red);
    for(int i=0;i<vLi->size();i++)
    {
        vert2d* v = vLi->at(i);
        assert(v->nedges==2);
        patr.drawPoint(v->x,v->y);
    }
        patr.end();
    img.save("/home/mprice/Desktop/Patch/LappedUtiloutput.png","PNG");*/


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
