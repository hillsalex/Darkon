#include "MPSandbox.h"
#include "drawengine.h"
#include "tamgenerator.h"
#include <sstream>
#include "assert.h"
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
    TAMGenerator tgen;
    int TAMtones = 6;
    int TAMsizes = 4;
    int TAMmaxw = 256;
    uchar** TAM = tgen.genTAM(TAMsizes,TAMtones,TAMmaxw);
    tgen.saveTAM("/home/mprice/Desktop/TAM/",TAM,TAMtones,TAMsizes,TAMmaxw);

}

//called every frame before draw.
void MPSandbox::update()
{
    m_framesElapsed++;
}

void MPSandbox::draw()
{
    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
}
