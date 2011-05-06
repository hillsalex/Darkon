#include <glm/glm.h>

#include "LappedUtils.h" //note to self: watch out for circular reference
#include <QSet>


#ifndef SURFACEPATCH_H
#define SURFACEPATCH_H

class SurfacePatch
{
public:
    SurfacePatch();
    QSet<vert2d*> verts;

};

#endif // SURFACEPATCH_H
