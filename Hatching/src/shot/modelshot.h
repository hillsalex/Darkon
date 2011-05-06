#ifndef MODELSHOT_H
#define MODELSHOT_H
#include <shot.h>
#include "meshoperator.h"

class modelShot : public Shot
{
public:
    modelShot(DrawEngine* parent,QHash<QString, QGLShaderProgram *>* shad, QHash<QString, GLuint>* tex, QHash<QString, Model>* mod);

    ~modelShot();

    //In begin, initialize things that could not have been initialized beforehand
    //(gl state)
    void begin();

    //called every frame before draw.
    void update();

    MeshOperator* m_operator;
    //draw!
    void draw();
};

#endif // MODELSHOT_H
