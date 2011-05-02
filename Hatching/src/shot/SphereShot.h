#ifndef SPHERESHOT_H
#define SPHERESHOT_H
#include <Sphere.h>
#include <shot.h>
#include <math/CS123Algebra.h>

class SphereShot : public Shot
{
public:
    SphereShot(DrawEngine* parent,QHash<QString, QGLShaderProgram *>* shad, QHash<QString, GLuint>* tex, QHash<QString, Model>* mod);
    //In begin, initialize things that could not have been initialized beforehand
    //(gl state)
    void begin();

    //called every frame before draw.
    void update();

    //draw!
    void draw();

    Sphere* sph;
};

#endif // SPHERESHOT_H
