#ifndef TAMGENERATOR_H
#define TAMGENERATOR_H
#include <QString>
#include <QImage>
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
    void drawqStroke(QImage* img,Stroke* stroke);
    float avgVal(uchar* data, int width);
    float qavgVal(QImage* img);
    uchar** genTAM(int sizes,int tones,int maxwidth);
    QImage** genImgTAM(int sizes, int tones, int maxwidth);
    QImage** ldImgTAM(QString path, int tones, int sizes);

    void saveTAM(QString path, uchar** TAM, int tones, int sizes, int maxw);
    void saveTAM(QString path, QImage** TAM, int tones, int sizes, int maxw);
};

#endif // TAMGENERATOR_H
