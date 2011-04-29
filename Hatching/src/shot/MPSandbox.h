#ifndef MPSANDBOX_H
#define MPSANDBOX_H

#include <shot.h>
#include <math/CS123Algebra.h>
class TAMstroke;


struct Stroke
{
    float x0,y0;
    float length;
    float x1,y1;
    float angle; //angle CCW from positive x in degrees
    Stroke(){}
};

class MPSandbox : public Shot
{
public:
    MPSandbox(DrawEngine* parent,QHash<QString, QGLShaderProgram *>* shad, QHash<QString, GLuint>* tex, QHash<QString, Model>* mod);

    ~MPSandbox();

    //In begin, initialize things that could not have been initialized beforehand
    //(gl state)
    void begin();

    //called every frame before draw.
    void update();

    //draw!
    void draw();

    float* brushKernel(int radius);
    void drawStroke(uchar* data, int width, float* kern, int strokerad, Stroke* stroke);
    float avgVal(uchar* data, int width);
    uchar** genTAM(int sizes,int tones,int maxwidth);
};



#endif // MPSANDBOX_H
