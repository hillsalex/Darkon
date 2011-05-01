#ifndef TAMGENERATOR_H
#define TAMGENERATOR_H
#include <QString>
struct Stroke
{
    float x0,y0;
    float length;
    float x1,y1;
    float angle; //angle CCW from positive x in degrees
    Stroke(){}
};

class TAMGenerator
{
public:
    TAMGenerator();

    float* brushKernel(int radius);
    void drawStroke(uchar* data, int width, float* kern, int strokerad, Stroke* stroke);
    float avgVal(uchar* data, int width);
    uchar** genTAM(int sizes,int tones,int maxwidth);
    void saveTAM(QString path, uchar** TAM, int tones, int sizes, int maxw);
};

#endif // TAMGENERATOR_H
