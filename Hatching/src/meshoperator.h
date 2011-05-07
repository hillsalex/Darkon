#ifndef MESHOPERATOR_H
#define MESHOPERATOR_H

#include "glm/glm.h"
#include "CS123Common.h"
#include <qgl.h>

class MeshOperator
{
public:
    MeshOperator();
    static void calculateCurvatures(GLMmodel* model);
};

#endif // MESHOPERATOR_H
