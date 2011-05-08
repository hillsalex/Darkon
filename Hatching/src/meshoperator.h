#ifndef MESHOPERATOR_H
#define MESHOPERATOR_H

#include "glm/glm.h"
#include "CS123Common.h"
#include <qgl.h>

class MeshOperator
{
public:
    MeshOperator();
    void calculateCurvatures(GLMmodel* model);
};

#endif // MESHOPERATOR_H
