#ifndef LAPPEDUTILS_H
#define LAPPEDUTILS_H
#include <QList>
#include <QImage>
#include <iostream>
#include <assert.h>
#include <CS123Common.h>
using std::cout;
using std::endl;

struct vert2d;
struct edge2d
{
    vert2d* v1;
    vert2d* v2;
    edge2d(vert2d* _v1, vert2d* _v2):v1(_v1),v2(_v2){}
};
struct vert2d
{
    int x,y;
    edge2d* e1;
    edge2d* e2;
    int nedges;
    vert2d():nedges(0){}
    vert2d(int _x,int _y):x(_x),y(_y),nedges(0){}
    edge2d* otherEdge(edge2d* in){if(in==e1){return e2;}else{return e1;}}
    void addEdge(edge2d* e){if(nedges==0){nedges++; e1 = e;}else if(nedges==1 && e1!=e){nedges++; e2=e;}else if(e1!=e && e2!=e){cout<<"UHOHMORETHANTWOEDGES"<<nedges<<endl;}}
    vert2d* otherAdjVert(vert2d* vin)
    {
        //e1 is old edge
        if(e1->v1 == vin || e1->v2 == vin)
        {
            if(e2->v1==this)
                return e2->v2;
            return e2->v1;
        }//e2 is old edge
        else if(e2->v1 == vin || e2->v2 == vin)
        {
            if(e1->v1==this)
                return e1->v2;
            return e1->v1;
        }
        else
        {
            cout<<"otherAdjVert nonAdjProblem: ("<<x<<","<<y<<") received ("<<vin->x<<","<<vin->y<<")"<<endl;
            print();
            vin->print();
        }
    }
    void replaceEdge(edge2d* _old, edge2d* _new)
    {
        if(e1==_old)
            e1=_new;
        else if(e2==_old)
            e2=_new;
        else
        {cout<<"TRIED TO REPLACE NONEXISTENT EDGE"<<endl;print();}
    }
    void print()
    {
        cout<<"("<<x<<","<<y<<") e1:{("<<e1->v1->x<<","<<e1->v1->y<<") ("<<e1->v2->x<<","<<e1->v2->y<<")} e2:{("<<e2->v1->x<<","<<e2->v1->y<<") ("<<e2->v2->x<<","<<e2->v2->y<<")}"<<endl;
    }
};

struct polyHull
{
    QList<vert2d*>* verts;
    QList<edge2d*>* edges;
    polyHull(QList<vert2d*>* vs, QList<edge2d*>* es):verts(vs),edges(es){}
};


class LappedUtils
{
public:
    LappedUtils();
    int ccw(vert2d* p1,vert2d* p2,vert2d* p3);
    void printEdge(edge2d* e){cout<<"{("<<e->v1->x<<","<<e->v1->y<<") ("<<e->v2->x<<","<<e->v2->y<<")}"<<endl;}
    polyHull* getPolyHull(QImage* blob,int iterations);

};

#endif // LAPPEDUTILS_H
