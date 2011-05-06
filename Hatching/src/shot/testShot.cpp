#include "testShot.h"
#include "drawengine.h"
#include "meshoperator.h"
using std::cout;
using std::endl;

testShot::testShot(DrawEngine* parent,QHash<QString, QGLShaderProgram *>* shad, QHash<QString, GLuint>* tex, QHash<QString, Model>* mod) : Shot(parent,shad,tex,mod)
{
    //lasts 150 frames
    m_lifespan = 150;
    m_operator =  new MeshOperator();
}

testShot::~testShot()
{

}

//In begin, initialize things that could not have been initialized beforehand
//(gl state)
void testShot::begin()
{

    QList<QString> keys = models_->keys();
    foreach (QString k,keys)
    {
        GLMmodel* model = models_->value(k).model;
        m_operator->calculateCurvatures(model);
    }
}

//called every frame before draw.
void testShot::update()
{
    m_framesElapsed++;
    if(m_framesElapsed>m_lifespan)
    {
        //end shot and call engine
    }
}

//draw!
extern "C"{
    extern void APIENTRY glActiveTexture(GLenum);
}
void testShot::draw()
{
    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glMatrixMode(GL_MODELVIEW);
 glActiveTexture(GL_TEXTURE0);
 QList<QString> keys = models_->keys();
 glBegin(GL_LINES);

 GLMmodel* model = models_->value("nail").model;
 GLfloat* curvatures = model->curvatures;
 GLfloat* vertices = model->vertices;
 for (int i=0;i<model->numvertices;i++)
 {
     glVertex3f(vertices[i*3], vertices[i*3+1],vertices[i*3+2]); // origin of the line
     glVertex3f(vertices[i*3]+curvatures[i*3], vertices[i*3+1]+curvatures[i*3+1],vertices[i*3+2]+curvatures[i*3+2]); // ending point of the line
 }

 glEnd( );
 glBindTexture(GL_TEXTURE_CUBE_MAP, textures_->value("cube_map_1"));
 shader_programs_->value(NAIL_SHADER)->bind();
 shader_programs_->value(NAIL_SHADER)->setUniformValue("CubeMap",0);
 shader_programs_->value(NAIL_SHADER)->setUniformValue("eyept",m_engine->camera_.eye.x, m_engine->camera_.eye.y, m_engine->camera_.eye.z);
glCallList(models_->value("nail").idx);
shader_programs_->value("reflect")->release();
}


