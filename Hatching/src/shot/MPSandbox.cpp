#include "MPSandbox.h"
#include "drawengine.h"
using std::cout;
using std::endl;

MPSandbox::MPSandbox(DrawEngine* parent,QHash<QString, QGLShaderProgram *>* shad, QHash<QString, GLuint>* tex, QHash<QString, Model>* mod) : Shot(parent,shad,tex,mod)
{

}

MPSandbox::~MPSandbox()
{

}

struct ARGB
{
    uchar a,r,g,b;
    ARGB(uchar _a, uchar _r, uchar _g, uchar _b) : a(_a),r(_r),g(_g),b(_b){}
    ARGB(uchar _r, uchar _g, uchar _b):a(255),r(_r),g(_g),b(_b){}
    ARGB():a(255),r(0),g(0),b(0){}
    ARGB(float _r, float _g, float _b):a(255),r(_r*255+0.5),g(_g*255+0.5),b(_b*255+0.5){}
};




//In begin, initialize things that could not have been initialized beforehand
//(gl state)
void MPSandbox::begin()
{
    int stroker = 1;
    float* kern = this->brushKernel(stroker);

    int imgw = 32;//always square
    uchar* imgdata = new uchar[4 * imgw * imgw];
    for(int i=0; i<4*imgw*imgw;i++)
        if(i%4==3)//BGRA
            imgdata[i]=255;
        else
            imgdata[i]=255;

    Stroke* s = new Stroke();
    s->angle=0;
    s->length=0.5;
    s->x0 = 0.5;
    s->y0 = 0.5;
    drawStroke(imgdata, imgw, kern, stroker,s);


    int TAMtones = 6;
    int TAMsizes = 4;
    uchar** TAM = genTAM(4,6,256);


    QImage* img = new QImage(imgdata,imgw,imgw,QImage::Format_ARGB32);//ACTUALLY BGRA?
    img->save("/home/mprice/Desktop/TAMTest","PNG");
    delete img;
}

//stroke: x0,y0,length,

void MPSandbox::drawStroke(uchar* data, int width, float* kern, int rad, Stroke* stroke)
{
    double ang = stroke->angle;
    double x0 = stroke->x0;
    double y0 = stroke->y0;
    double slen = stroke->length;
    double vX1 = slen*cos(ang * M_PI/180.0 );
    double vY1 = slen*sin(ang * M_PI/180.0);

//number of increments
    int ninc = 2*slen*width/rad;

    for(int i=0; i<=ninc; i++)
    {
        double t = slen*i/ninc;
        int x = (int)((x0 + t*cos(ang * M_PI/180.0))*width) % width;
        int y = (int)((y0 + t*sin(ang * M_PI/180.0))*width) % width;

        for(int r=0;r<rad*2+1;r++)
        {
            for(int c=0;c<rad*2+1;c++)
            {
                float alph = kern[r*(2*rad+1)+c];

                int _x = x-rad+c % width;
                int _y = y-rad+r % width;
                {
                    data[4 * (_y * width + _x) + 0] = alph*0 + (1.0-alph)*data[4 * (_y * width + _x) + 0];
                    data[4 * (_y * width + _x) + 1] = alph*0 + (1.0-alph)*data[4 * (_y * width + _x) + 1];
                    data[4 * (_y * width + _x) + 2] = alph*0 + (1.0-alph)*data[4 * (_y * width + _x) + 2];
                }
            }
        }
    }
}

//assumes ARGB greyscale input so uses one color channel
float MPSandbox::avgVal(uchar *data, int width)
{
    float val=0;
    for(int r=0;r<width;r++)
    {
        for(int c=0;c<width;c++)
        {
            val+=data[4*(r*width+c)+2];
        }
    }
    return val/(width*width);
}

uchar** MPSandbox::genTAM(int sizes,int tones,int maxwidth)
{
//Allocate TAM as 2d array of pointers to images(uchar arrays)
//organized as appears in paper
uchar** TAM = new uchar*[tones * sizes];
for(int r=0;r<tones;r++)
    {
    int w=maxwidth;
    for(int c=sizes-1;c>=0;c--)
        {
            TAM[r * tones + c] = new uchar[4*w];
            w/=2;
        }
    }

//for tone level
    //for size
        //while average tone<desired tone
            //1000 times:
                //generate candidate stroke
                //for each size determine goodness of stroke, sum
                //normalize by length
            //determine and add best stroke, delete rest
            //update average tone
}


float* MPSandbox::brushKernel(int radius)
{
    int w = radius*2+1;
    int h = w;
    float* kern = new float[w*h];
    for(int r=0;r<w;r++)
    {
        for(int c=0; c<h;c++)
        {
            int dx=c-radius;
            int dy=r-radius;
            float dist = sqrt((float)(dx*dx+dy*dy));
            if(dist<radius*1.0+0.0000001)
                //linear falloff w/fudging
                kern[r*(radius*2+1)+c]=1.0-.99*dist/radius;
            else
                kern[r*(radius*2+1)+c]=0.0;
            //cout<<kern[r*(radius*2+1)+c]<<"\t";
        }
        //cout<<endl;
    }
    return kern;
}


//called every frame before draw.
void MPSandbox::update()
{
    m_framesElapsed++;
}

void MPSandbox::draw()
{
    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

}
