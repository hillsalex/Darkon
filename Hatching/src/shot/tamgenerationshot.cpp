#include "tamgenerationshot.h"


TamGenerationShot::TamGenerationShot(DrawEngine* parent,QHash<QString, QGLShaderProgram *>* shad, QHash<QString, GLuint>* tex, QHash<QString, Model>* mod) : Shot(parent,shad,tex,mod)
{

}



void TamGenerationShot::update()
{

}


void TamGenerationShot::draw()
{

}

void TamGenerationShot::begin()
{
    int sizes = 4;
    int tones = 6;
    int maxwidth = 256;
TAMGenerator tgen;
QImage** TAM = tgen.genImgTAM(sizes, tones, maxwidth);
tgen.saveTAM("/home/mprice/Desktop/newTAM/", TAM, 6, 4, 256);
}
