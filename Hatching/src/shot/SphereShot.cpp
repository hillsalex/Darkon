

#include "SphereShot.h"
#include "drawengine.h"
#include <QGLShaderProgram>
#include <QGLShader>
#include <QHash>
#include "tamgenerator.h"
#include "Sphere.h"
#include "assert.h"
using std::cout;
using std::endl;
using std::vector;
using std::stringstream;
SphereShot::SphereShot(DrawEngine* parent,QHash<QString, QGLShaderProgram *>* shad, QHash<QString, GLuint>* tex, QHash<QString, Model>* mod) : Shot(parent,shad,tex,mod)
{
     sph = new Sphere(25,25);
}



//In begin, initialize things that could not have been initialized beforehand
//(gl state)
void SphereShot::begin()
{
    TAMGenerator tgen;
    int TAMtones = 6;
    int TAMsizes = 4;
    int TAMmaxw = 256;
    //uchar** TAM = tgen.genTAM(TAMsizes,TAMtones,TAMmaxw);
    //tgen.saveTAM("/home/mprice/Desktop/TAM/",TAM,TAMtones,TAMsizes,TAMmaxw);
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
   glScalef(4,4,4);
   glTranslatef(0.5, 0.5, 0);
   glRotatef(90.0, 0, 0, 1);
   glTranslatef(-0.5, -0.5, 0);
   glMatrixMode(GL_MODELVIEW);

   glEnable(GL_TEXTURE_2D);

   GLuint tid[6];
   glGenTextures(6,tid);

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
}

//called every frame before draw.
void SphereShot::update()
{
    m_framesElapsed++;
}

void SphereShot::draw()
{
    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
    glEnable(GL_TEXTURE_2D);

    float lightpos[4];
    lightpos[0]=1.2f;
    lightpos[1]=-2.f * sin(m_framesElapsed/10.0);
    lightpos[2]=-1.f;
    lightpos[3]=0.f;
    glLightfv(GL_LIGHT0,GL_POSITION,lightpos);

   sph->draw();

   glMatrixMode(GL_MODELVIEW);
   glTranslatef(2,0,0);
   glRotatef(45,0,0,1);

   sph->draw();

}
