#include "tamgenerator.h"
#include <QImage>
#include <sstream>
#include <iostream>
#include <QPainter>
#include "math.h"
#include "assert.h"
using std::cout;
using std::endl;
using std::vector;
using std::stringstream;
TAMGenerator::TAMGenerator()
{
}



void TAMGenerator::saveTAM(QString path, uchar** TAM, int tones, int sizes, int maxw)
{cout<<"Saving TAM in "<<path.toStdString()<<"...";
    for(int t=0; t<tones; t++)
    {
        for(int sz=0;sz<sizes;sz++)
        {
            stringstream ss;
            ss<<path.toStdString()<<"TAM_"<<t<<"_"<<sz;
            int cursize = maxw >> (sizes-sz-1);
            QImage* curimg = new QImage(TAM[sz*tones+t],cursize,cursize,QImage::Format_ARGB32);
            curimg->save(QString::fromStdString(ss.str()),"PNG");
            delete curimg;
        }
    }cout<<"done."<<endl;
}
void TAMGenerator::saveTAM(QString path, QImage** TAM, int tones, int sizes, int maxw)
{cout<<"Saving TAM in "<<path.toStdString()<<"...";
    for(int t=0; t<tones; t++)
    {
        for(int sz=0;sz<sizes;sz++)
        {
            stringstream ss;
            ss<<path.toStdString()<<"TAM_"<<t<<"_"<<sz;
            int cursize = maxw >> (sizes-sz-1);
            QImage* curimg = TAM[tones*sz+t];
            curimg->save(QString::fromStdString(ss.str()),"PNG");
            delete curimg;
        }
    }cout<<"done."<<endl;
}

QImage** TAMGenerator::ldImgTAM(QString path, int tones, int sizes)
{
    QImage** toreturn = new QImage*[tones*sizes];
    for(int t=0;t<tones;t++)
    {
        for(int s=0;s<sizes;s++)
        {
            stringstream ss;
            ss<<path.toStdString()<<"TAM_"<<t<<"_"<<s;
            QImage* curtam = new QImage(QString::fromStdString(ss.str()));
            toreturn[s*tones+t] = curtam;
        }
    }
    return toreturn;
}

void TAMGenerator::drawqStroke(QImage* img,Stroke* stroke)
{
    double ang = stroke->angle;
    int x0 = stroke->x0 * img->width();
    int y0 = stroke->y0 * img->height();

    double slen = stroke->length;
    int dX1 = slen*cos(ang * M_PI/180.0) * img->width();
    int dY1 = slen*sin(ang * M_PI/180.0) * img->height();
    int x2 = x0+dX1;
    int y2 = y0+dY1;
    //first just draw line x0y0 --- (x2,y2) (it'll cut off fine)
    //then draw line (x2%w, y2%h) ---- (x2%w - dX1 , y2%h - dY1)
    QPainter* ptr = new QPainter(img);



    ptr->setPen(QPen(
            QBrush(QColor(100,100,100)), 3.0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin
            ));

    if(x0==x2 && y0 == y2)
    {
        ptr->drawPoint(x0,y0);
        //ptr->drawEllipse(x0,y0,8,8);
    }
    else
    {
        ptr->drawLine(x0,y0,x2,y2);
        ptr->drawLine(x2%img->width(), y2%img->height(), x2%img->width() - dX1, y2%img->height() - dY1);
        //ptr->end();
    }

    ptr->setPen(QPen(
            QBrush(QColor(50,50,50)), 2.0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin
            ));

    if(x0==x2 && y0 == y2)
    {
        ptr->drawPoint(x0,y0);
        //ptr->drawEllipse(x0,y0,8,8);
    }
    else
    {
        ptr->drawLine(x0,y0,x2,y2);
        ptr->drawLine(x2%img->width(), y2%img->height(), x2%img->width() - dX1, y2%img->height() - dY1);
        //ptr->end();
    }

    //PASS 2
    ptr->setPen(QPen(
            QBrush(Qt::black), 1.0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin
            ));

    if(x0==x2 && y0 == y2)
    {
        ptr->drawPoint(x0,y0);
        //ptr->drawEllipse(x0,y0,8,8);
    }
    else
    {
        ptr->drawLine(x0,y0,x2,y2);
        ptr->drawLine(x2%img->width(), y2%img->height(), x2%img->width() - dX1, y2%img->height() - dY1);
        //ptr->end();
    }
    ptr->end();


    delete ptr;

}

//stroke: x0,y0,length,
void TAMGenerator::drawStroke(uchar* data, int width, float* kern, int rad, Stroke* stroke)
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
//returns float 0 - 255???
float TAMGenerator::avgVal(uchar *data, int width)
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

float TAMGenerator::qavgVal(QImage* img)
{
    return avgVal(img->bits(), img->width());
}

QImage** TAMGenerator::genImgTAM(int sizes, int tones, int maxwidth)
{
    QImage** TAM = new QImage*[tones*sizes];
    QList<Stroke*>** strokeLists = new QList<Stroke*>*[tones*sizes];
    //Allocate TAM
    for(int r=0; r<tones;r++)
    {
        int w = maxwidth;
        for(int c=sizes-1;c>=0;c--)
        {
            TAM[c*tones+r] = new QImage(w,w,QImage::Format_ARGB32);
            strokeLists[c*tones+r] = new QList<Stroke*>();
            w/=2;
        }
    }

    for(int tlev=0; tlev<tones; tlev++)
    {
        float curtone = (tlev+1)*(1.0/(tones+1));
        for(int slev=0; slev<sizes; slev++)
        {
            int cursize = maxwidth >> (sizes-slev-1);
            QImage* curtam = TAM[slev*tones+tlev];

            if(tlev==0)//if first column, set to white initially
                memset(curtam->bits(),255,sizeof(uchar)*cursize*cursize*4);
            else//otherwise copy memory from guy on the left
                memcpy(curtam->bits(), TAM[slev*tones+tlev-1]->bits(),sizeof(uchar)*cursize*cursize*4);

            //draw unique strokes from all guys on top
            for(int i=0; i<slev; i++)
            {
                QList<Stroke*>* strokestoadd = strokeLists[i*tones+tlev];
                for(int sid=0; sid<strokestoadd->size();sid++)
                {
                    drawqStroke(curtam,strokestoadd->at(sid));
                }
            }

            //then generate and draw strokes necessary to add tone
            float av = qavgVal(curtam);
            QImage** guesses = new QImage*[sizes];
            int _w = cursize;
            for(int _s=slev; _s>=0; _s--)
            {
                guesses[_s] = new QImage(_w,_w,QImage::Format_ARGB32);
                _w/=2;
            }

           // QImage* guess = new QImage(cursize,cursize, QImage::Format_ARGB32);
            int kkk=0;
            while(av > 255.0 - 255*curtone*2)
            {
                float maxgoodness = 0;
                Stroke* bestStroke;
                for(int i=0; i<1; i++)//# trials
                {kkk++; //cout<<kkk<<endl;
                    int saz = cursize;
                    float goodness =0;

                    Stroke* s = new Stroke();
                    if(tlev<2)
                        s->angle = rand()%7 -3;
                    else
                        s->angle = rand()%7 - 3 + 90;
                    s->length = (rand()%8+3)/10.0;
                    s->x0 = (rand()%1000)/1000.0;
                    s->y0 = (rand()%1000)/1000.0;

                    for(int _sz=slev; _sz>=0; _sz--)
                    {
                        memcpy(guesses[_sz]->bits(), TAM[_sz*tones+tlev]->bits(),sizeof(uchar)*saz*saz*4);
                        saz/=2;
                        drawqStroke(guesses[_sz],s);
                        float gav = qavgVal(guesses[_sz]);
                        goodness += qavgVal(TAM[_sz*tones+tlev]) - gav;
                    }
                        goodness /= s->length;
                        if(goodness>=maxgoodness)
                        {
                            bestStroke = s;
                            maxgoodness = goodness;
                        }
                        else
                        {
                            delete s;
                        }
                }
                    QList<Stroke*>* curStrLst = strokeLists[slev*tones+tlev];
                    curStrLst->append(bestStroke);
                    drawqStroke(curtam,bestStroke);
                    av = qavgVal(curtam);
            }
            cout<<"yay"<<endl;
            for(int _s=slev; _s>=0; _s--)
            {
                delete guesses[_s];
            }
            delete[] guesses;
        }
    }
    return TAM;
}

uchar** TAMGenerator::genTAM(int sizes,int tones,int maxwidth)
{
//Allocate TAM as 2d array of pointers to images(uchar arrays)

/*    light small   #   #   #   #   dark small
 *         #        #   #   #   #   #
 *         #        #   #   #   #   #
 *      light big   #   #   #   #   dark big
 */
uchar** TAM = new uchar*[tones * sizes];
vector<Stroke*>** strokeLists = new vector<Stroke*>*[tones*sizes];

int stroker = 1;
float* kern = this->brushKernel(stroker);

//ALLOCATE TAM
for(int r=0;r<tones;r++)
    {
    int w=maxwidth;
    for(int c=sizes-1;c>=0;c--)
        {
            TAM[c * tones +r] = new uchar[4*w*w];
            strokeLists[c*tones+r] = new vector<Stroke*>();
            w/=2;
        }
    }

    for(int tlev=0; tlev<tones; tlev++)
        {
        float curtone = (tlev+1)*(1.0/(tones+1));
        for(int slev=0;slev<sizes;slev++)
            {
                int cursize = maxwidth >> (sizes-slev-1);
                uchar* curtam = TAM[slev*tones+tlev];

                if(tlev==0)//if first column, set to white initially
                    memset(curtam,255,sizeof(uchar)*cursize*cursize*4);
                else//otherwise copy memory from guy on the left
                    memcpy(curtam,TAM[slev*tones+tlev-1],sizeof(uchar)*cursize*cursize*4);

                cout<<"about to draw strokes on top"<<endl;
                //draw unique strokes from all guys on top
                for(int i=0;i<slev;i++)
                {
                    vector<Stroke*>* strokestoadd = strokeLists[i*tones+tlev];
                    cout<<"size is: "<<strokestoadd->size()<<endl;
                    for(int sid=0; sid<strokestoadd->size(); sid++)
                    {   cout<<"drawing stroke "<<sid<<endl;
                        drawStroke(curtam,cursize,kern,stroker,strokestoadd->at(sid));
                    }
                }
                cout<<"done drawing strokes on top"<<endl;

                //then generate and draw strokes necessary to add tone
                float av = avgVal(curtam,cursize);
                while(av>255.0 - 255*curtone)
                {
                    cout<<"making stroke"<<endl;
                    Stroke* s = new Stroke();
                    cout<<"assigning stroke properties"<<endl;
                    s->angle = rand()%7 -3;
                    s->length = (rand()%8+3)/10.0;
                    s->x0 = (rand()%1000)/1000.0;
                    s->y0 = (rand()%1000)/1000.0;

                    vector<Stroke*>* curStrLst = strokeLists[slev*tones+tlev];
                    cout<<"done"<<endl;
                    cout<<"stroke list at (" << tlev <<","<< slev <<") "<< curStrLst <<" size: "<<curStrLst->size()<<endl;
                    cout<<"["; for(int bl=0;bl<curStrLst->size();bl++){cout<<(curStrLst->at(bl))->x0<<",";} cout<<"]"<<endl;
                    curStrLst->push_back(s);
                    cout<<"inserted stroke of x0 "<<s->x0<<endl;
                    drawStroke(curtam,cursize,kern,stroker,s);
                    cout<<"..."<<endl;
                    av = avgVal(curtam,cursize);
                }
            }
        }
    return TAM;
}


float* TAMGenerator::brushKernel(int radius)
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
        }
    }
    return kern;
}
