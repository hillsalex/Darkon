#ifndef LAPPEDUTILS_H
#define LAPPEDUTILS_H
#include <QList>
#include <QImage>
#include <iostream>
#include <assert.h>
#include <math/CS123Algebra.h>
#include <glm.h>
#include <QHash>
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
    int imgw;
    int imgh;
    polyHull(QList<vert2d*>* vs, QList<edge2d*>* es):verts(vs),edges(es){}
    //returns true if the given segment intersects the given edge
    bool isectEdge(edge2d* e, float x0, float y0, float x1, float y1)
    {  // cout<<"isectEdge "<<x0<<","<<y0<<" "<<x1<<","<<y1<<" edge2d: "<<e->v1->x<<","<<e->v1->y <<") (" << e->v2->x << "," << e->v2->y << ")" << endl;
        if(x0==x1 && e->v1->x == e->v2->x)return false;//both vertical
        if(y0==y1 && e->v1->y == e->v2->y)return false;//both horizontal
        int isecx,isecy;
        if(x0==x1)
        {//given segment is vertical
            int Me = (e->v2->y - e->v1->y)/(e->v2->x - e->v1->x);
            int Be = e->v1->y - Me * e->v1->x;
            isecx = x0;
            isecy = Me * isecx + Be;
        }
        else if(e->v1->x == e->v2->x)
        {//edge is vertical
            int Ms = (y1-y0)/(x1-x0);
            int Bs = y0 - Ms * x0;
            isecx = e->v1->x;
            isecy = Ms * isecx + Bs;
        }
        else
        {//otherwise no special case
            int Ms = (y1-y0)/(x1-x0);
            int Me = (e->v2->y - e->v1->y)/(e->v2->x - e->v1->x);
            if(Ms==Me)return false;
            int Be = e->v1->y - Me * e->v1->x;
            int Bs = y0 - Ms * x0;
            isecx = (Bs - Be)/(Me - Ms);
            isecy = Me * isecx + Be;
        }
        int xmin,xmax,ymin,ymax;
        //first test edge bounds
        if(e->v1->x<=e->v2->x){xmin=e->v1->x; xmax=e->v2->x;}else{xmin=e->v2->x; xmax=e->v1->x;}
        if(e->v1->y<=e->v2->y){ymin=e->v1->y; ymax=e->v2->y;}else{ymin=e->v2->y; ymax=e->v1->y;}
        if(isecx<xmin || isecx>xmax || isecy<ymin || isecy>ymax)return false;
        //then segment bounds
        if(x0<=x1){xmin=x0; xmax=x1;}else{xmin=x1; xmax=x0;}
        if(y0<=y1){ymin=y0; ymax=y1;}else{ymin=y1; ymax=y0;}
        if(isecx<xmin || isecx>xmax || isecy<ymin || isecy>ymax)return false;
        return true;
    }
    //returns true if the given segment intersects any edge in the hull
    bool isectAnyEdge(float x0, float y0, float x1, float y1)
    {
        for(int i=0;i<edges->size();i++)
        {//just check intersection for every edge
            if(isectEdge(edges->at(i),x0,y0,x1,y1))return true;
        }return false;
    }
    //returns true if the point is in the interior of the hull
    bool isInteriorPt(float x, float y)
    {
        //cout<<"chkinterior: "<<x<<","<<y<<endl;
        //use segment (0,0) --> (x,y)
        int numIntersections=0;
        for(int i=0;i<edges->size();i++)
        {//just check intersection for every edge
            if(isectEdge(edges->at(i),0,0,x,y))
            {//cout<<"isect ("<< edges->at(i)->v1->x <<","<<edges->at(i)->v1->y <<") (" << edges->at(i)->v2->x << "," << edges->at(i)->v2->y << ")" << endl;
                numIntersections++;
            }
        }
        //cout<<"numIntersections: "<<numIntersections<<endl;
        return numIntersections%2==1;
    }
    //returns true if the segment intersects the hull
    bool isectHull (float x0,float y0,float x1, float y1)
    {//A segment intersects the hull if A: it intersects an edge of the hull or B: one point is inside the hull
        return(isectAnyEdge(x0,y0,x1,y1) || isInteriorPt(x0,y0));
    }
    //basically the same as above, but first converts the UV coords into texture space coords (ints)
    bool isectHullUV(float s0, float t0, float s1, float t1)
    {   bool a,b,c;
        a = isectAnyEdge(s0*imgw,(1.0-t0)*imgh, s1*imgw, (1.0-t1)*imgh);
        b = isInteriorPt(s0*imgw,(1.0-t0)*imgh);
        c = isInteriorPt(s1*imgw,(1.0-t1)*imgh);
        if(a)cout<<"edge intersected"<<endl;
        if(b)cout<<"v0 inside"<<endl;
        if(c)cout<<"v1 inside"<<endl;
        return a||b||c;
    }

    void print()
    {   cout<<"PRINTING HULL"<<endl;
        cout<<"verts: ("<<verts->size()<<")"<<endl;
        for(int i=0; i<verts->size(); i++)
        {
            cout << "(" << verts->at(i)->x << "," << verts->at(i)->y <<")"<<endl;
            verts->at(i)->print();
        }
        cout<<"edges: ("<<edges->size()<<")"<<endl;
        for(int i=0; i<edges->size(); i++)
            cout<<"("<< edges->at(i)->v1->x <<","<< edges->at(i)->v1->y <<") --- (" <<  edges->at(i)->v2->x <<","<<  edges->at(i)->v2->y << ")"<<endl;
    }
};

struct PatchTri;
struct PatchVert
{
    //the texture coordinates (we assign these)
    float s,t;//really should be outside since vary perpatch


    //the index of the x coordinate in the GLM vertex array
    int GLMidx;
    //the coordinate values
    //float x,y,z;
    Vector4 pos;
    //the triangles that contain this vertex
    QList<PatchTri*>* tris;
    PatchVert(){}
};

struct PatchEdge;
struct PatchTri
{
    //the patchverts we made
    PatchVert* v0;
    PatchVert* v1;
    PatchVert* v2;
    //edges
    PatchEdge* e01;
    PatchEdge* e12;
    PatchEdge* e20;

    //the GLM triangle (maybe some redundant information but w/e)
    _GLMtriangle* GLMtri;
    Vector4 tangent;

    PatchVert* otherVert(PatchVert* _A, PatchVert* _B)
    {
        if(v0!= _A && v0!= _B)
            return v0;
        if(v1!= _A && v1!= _B)
            return v1;
        if(v2!= _A && v2!= _B)
            return v2;
        return NULL;
    }
};

struct PatchEdge
{
    PatchVert* v0;
    PatchVert* v1;
    int ntris;
    PatchTri* t1;
    PatchTri* t2;
    PatchEdge(PatchVert* _v0, PatchVert* _v1):v0(_v0),v1(_v1){ntris=0;}
    void addTri(PatchTri* t){if(ntris==0){ntris++;t1=t;}else if(ntris==1){ntris++;t2=t;}else{cout<<"TRIED TO ADD MORE THAN TWO TRIANGLES TO AN EDGE"<<endl;}}
};

struct LappedPatch
{
    PatchTri* seed;
    QList<PatchTri*>* tris;
    QHash<PatchVert*, vec2<float> >* uvs;
};

class LappedUtils
{
public:
    LappedUtils();
    int ccw(vert2d* p1,vert2d* p2,vert2d* p3);
    void printEdge(edge2d* e){cout<<"{("<<e->v1->x<<","<<e->v1->y<<") ("<<e->v2->x<<","<<e->v2->y<<")}"<<endl;}
    polyHull* getPolyHull(QImage* blob,int iterations);
    QList<LappedPatch*>* generatePatches(GLMmodel* model, polyHull* polyhull);
    vec2<float> estimateUV(PatchVert* A, PatchVert* B, PatchVert* C, vec2<float> Ast, vec2<float> Bst);
    void assignSeedUV(PatchTri* seed, vec2<float> &v0st, vec2<float> &v1st, vec2<float> &v2st);
    void vizualizePatch(LappedPatch* patch, QImage* img);
    void drawEdgeFromUV(QImage* img, QPainter* patr, vec2<float> v0, vec2<float>v1);
};

#endif // LAPPEDUTILS_H
