#ifndef TESTSHOT_H
#define TESTSHOT_H
#include <shot.h>
#include "meshoperator.h"

class testShot : public Shot
{
public:
    testShot(DrawEngine* parent,QHash<QString, QGLShaderProgram *>* shad, QHash<QString, GLuint>* tex, QHash<QString, Model>* mod);

    ~testShot();

    //In begin, initialize things that could not have been initialized beforehand
    //(gl state)
    void begin();

    //called every frame before draw.
    void update();

    MeshOperator* m_operator;
    //draw!
    void draw();
};

#endif // TESTSHOT_H
