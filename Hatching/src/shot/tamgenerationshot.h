#ifndef TAMGENERATIONSHOT_H
#define TAMGENERATIONSHOT_H
#include "tamgenerator.h"
#include "shot.h"
class TamGenerationShot : public Shot
{
public:
    TamGenerationShot(DrawEngine* parent,QHash<QString, QGLShaderProgram *>* shad, QHash<QString, GLuint>* tex, QHash<QString, Model>* mod);

    //In begin, initialize things that could not have been initialized beforehand
    //(gl state)
    void begin();

    //called every frame before draw.
    void update();

    //draw!
    void draw();
};

#endif // TAMGENERATIONSHOT_H
