#ifndef MODELSHOT_H
#define MODELSHOT_H
#include <shot.h>
#include "meshoperator.h"
#include "math/CS123Algebra.h"

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
    //draw!
    void draw();
protected:

    void renderNormal(const Vector4 &vertex, const Vector4 &direction);
            MeshOperator* m_operator;

};

#endif // MODELSHOT_H
