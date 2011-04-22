#include "testShot.h"
#include "drawengine.h"
using std::cout;
using std::endl;

testShot::testShot(DrawEngine* parent,QHash<QString, QGLShaderProgram *>* shad, QHash<QString, GLuint>* tex, QHash<QString, Model>* mod) : Shot(parent,shad,tex,mod)
{
    //lasts 150 frames
    m_lifespan = 150;
}

testShot::~testShot()
{

}

//In begin, initialize things that could not have been initialized beforehand
//(gl state)
void testShot::begin()
{
glShadeModel(GL_SMOOTH);
glClearColor(0.8f,0.8f,0.8f,0.0f);
glEnable(GL_LIGHTING);
//Make some lights
glEnable(GL_LIGHT0);
float lightpos[4];
lightpos[0]=0.f;
lightpos[1]=2.f;
lightpos[2]=1.f;
lightpos[3]=0.f;
glLightfv(GL_LIGHT0,GL_POSITION,lightpos);
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
 glBindTexture(GL_TEXTURE_CUBE_MAP, textures_->value("cube_map_1"));
 shader_programs_->value(NAIL_SHADER)->bind();
 shader_programs_->value(NAIL_SHADER)->setUniformValue("CubeMap",0);
 shader_programs_->value(NAIL_SHADER)->setUniformValue("eyept",m_engine->camera_.eye.x, m_engine->camera_.eye.y, m_engine->camera_.eye.z);
glCallList(models_->value("nail").idx);
shader_programs_->value("reflect")->release();
}


