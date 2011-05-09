#include "LappedUtils.h"
#include <QQueue>
#include <QSet>
#include "assert.h"

LappedUtils::LappedUtils()
{
}

//from the wikipedia article for graham scan, I was lazy
/* Three points are a counter-clockwise turn if ccw > 0, clockwise if
 ccw < 0, and collinear if ccw = 0 because ccw is a determinant that
 gives the signed area of the triangle formed by p1, p2 and p3.*/
int LappedUtils::ccw(vert2d* p1,vert2d* p2,vert2d* p3)
{
    float f =(p2->x - p1->x)*(p3->y - p1->y) - (p2->y - p1->y)*(p3->x - p1->x);
    if(EQ(f,0.0))
        return 0;
    else if(f>0)
        return 1;
    else
        return -1;
}

//Assuming A and B have correct UV coords and world space coords, return an estimate of C's UV coords
vec2<float> LappedUtils::estimateUV(PatchVert* A, PatchVert* B, PatchVert* C, vec2<float> Ast, vec2<float> Bst)
{


    Vector4 AC(C->pos.x - A->pos.x, C->pos.y - A->pos.y, C->pos.z - A->pos.z,0);
    Vector4 AB(B->pos.x - A->pos.x, B->pos.y - A->pos.y, B->pos.z - A->pos.z,0);

    Vector4 AB_r = AB * getRotMat(Vector4(A->pos.x,A->pos.y,A->pos.z,0),AC.cross(AB),M_PI/2.0);
    float x = AC.dot(AB) / AB.getMagnitude2();
    float y = AC.dot(AB_r) / AB.getMagnitude2();

    //these are UVS
    vec2<float> ABp,AB_rp;
    ABp.x = Bst.x - Ast.x;
    ABp.y = Bst.y - Ast.y;
    AB_rp.x = ABp.y;
    AB_rp.y = -ABp.x;



    vec2<float> Cp; //Ap,Bp
   // Ap.x = Ast.x; Ap.y = Ast.y; Bp.x = B->s; Bp.y = B->t;
    Cp = Ast + x*ABp + y*AB_rp;

    //cout<<"AB: "<<AB<<" AC: "<<AC<<" AB_r: "<<AB_r<<" ABp: "<<ABp<<" AB_rp: "<<AB_rp<<" x: "<<x<<" y: "<<y<<" Ap:"<<Ap<<" Bp: "<<Bp<<" Cp: "<<Cp<<endl;
    //cout<<"AC dot AB: "<<AC.dot(AB)<<" AB dot AC: "<<AB.dot(AC)<<" mag AB: "<<AB.getMagnitude()<<endl;
    return Cp;
}


void printVector4(Vector4* v)
{

    cout << '[' << v->x << " , " << v->y << " , " << v->z << ',' << v->w << ']' << endl;
}


void LappedUtils::assignSeedUV(PatchTri* seed, vec2<float> &v0st, vec2<float> &v1st, vec2<float> &v2st)
{
    float scale = 0.5;
    Vector4 A,B,C;
    A = seed->v0->pos;
    B = seed->v1->pos;
    C = seed->v2->pos;
    A.w = 1;
    B.w = 1;
    C.w = 1;

    Vector4 ctr = (A+B+C)/3.0;
    ctr.w=0;
    Vector4 Ap, Bp, Cp, Tp;
    Matrix4x4 transMat = getTransMat(-ctr);
    Ap = transMat*A;
    Bp = transMat*B;
    Cp = transMat*C;

    Ap.w=0;
    Bp.w=0;
    Cp.w=0;

    Vector4 norm = ((Bp-Ap).cross(Cp-Ap)).getNormalized();
    printVector4(norm);
    double angle = acos(norm.dot(Vector4(0,0,1.0,0)));
    cout << "angle between normal and plane z=1: " << angle << endl;
    Matrix4x4 rotMat = getRotMat(Vector4(0,0,0,0),norm.cross(Vector4(0,0,1.0,0)).getNormalized(),angle);
    Ap.w = 1;
    Bp.w = 1;
    Cp.w = 1;

    Ap = rotMat*Ap;
    Bp = rotMat*Bp;
    Cp = rotMat*Cp;
    Tp = rotMat*seed->tangent;
    Ap.w = 0;
    Bp.w = 0;
    Cp.w = 0;


    double avgLen = max(max((Ap-Cp).getMagnitude(),(Ap-Bp).getMagnitude()),(Cp-Bp).getMagnitude());
    Ap = Ap*scale/avgLen;
    Bp = Bp*scale/avgLen;
    Cp = Cp*scale/avgLen;
    Ap.w = 1;
    Bp.w = 1;
    Cp.w = 1;
    Tp.w = 0;

    seed->tangent = Tp;
    Tp.normalize();
    double tanAngle = acos(Vector4(0,1,0,0).dot(Tp));
    Matrix4x4 tanMat = getRotMat(Vector4(0,0,0,0),Vector4(0,1,0,0).cross(Tp),tanAngle);
    Ap = tanMat*Ap;
    Bp = tanMat*Bp;
    Cp = tanMat*Cp;
    Ap.w = 1;
    Bp.w = 1;
    Cp.w = 1;

    transMat = getTransMat(Vector4(.5,.5,0,0));
    Ap = transMat*Ap;
    Bp = transMat*Bp;
    Cp = transMat*Cp;


    v0st.x = Ap.x;
    v0st.y = Ap.y;
    v1st.x = Bp.x;
    v1st.y = Bp.y;
    v2st.x = Cp.x;
    v2st.y = Cp.y;
}

polyHull* LappedUtils::getPolyHull(QImage* blob,int iterations)
{
    QImage img = *blob;
    int w = img.width();
    int h = img.height();

    uchar* imgdata = img.bits();

    QQueue<vert2d*>* vQu = new QQueue<vert2d*>();
    QList<vert2d*>* vLi = new QList<vert2d*>();//just for debugging/drawing
    QList<edge2d*>* eLi = new QList<edge2d*>();
    QHash<QPair<vert2d*,vert2d*>,edge2d*>* edgesMade = new QHash<QPair<vert2d*,vert2d*>,edge2d*>();
    QHash<QPair<int,int>,vert2d*>* vertsMade = new QHash<QPair<int,int>,vert2d*>();

    vert2d* vseed;
    for(int x=1;x<w-1;x++)
    {bool br=false;
        for(int y=1;y<h-1;y++)
        {
            if(imgdata[4*(y*w+x)]<255 && (imgdata[4*((y+1)*w+x)]==255 || imgdata[4*((y-1)*w+x)]==255 || imgdata[4*(y*w+(x+1))]==255 || imgdata[4*(y*w+(x-1))]==255))
            {
                vseed = new vert2d(x,y);
                vQu->append(vseed);
                vLi->append(vseed);
                vertsMade->insert(qMakePair(x,y),vseed);
                br=true;
                break;
            }
        }if(br)break;
    }

    while(!vQu->empty())
    {
        vert2d* v = vQu->dequeue();
        //cout<<"dequeued ("<<v->x<<","<<v->y<<")"<<endl;
        for(int x2=v->x-1;x2<=v->x+1;x2++)
        {
            for(int y2=v->y-1;y2<=v->y+1;y2++)
            {
                if(x2==v->x && y2==v->y){continue;}
                if(imgdata[4*(y2*w+x2)]<255 && (imgdata[4*((y2+1)*w+x2)]==255 || imgdata[4*((y2-1)*w+x2)]==255 || imgdata[4*(y2*w+(x2+1))]==255 || imgdata[4*(y2*w+(x2-1))]==255))
                {
                    vert2d* v2;
                    if(!vertsMade->contains(qMakePair(x2,y2)))
                    {
                        v2 = new vert2d(x2,y2);
                        vQu->append(v2);
                        vLi->append(v2);
                        vertsMade->insert(qMakePair(x2,y2),v2);
                    }
                    else
                    {
                        v2 = vertsMade->value(qMakePair(x2,y2));
                    }
                    if(!edgesMade->contains(qMakePair(v,v2)))
                    {
                        edge2d* e = new edge2d(v,v2);
                        edgesMade->insert(qMakePair(v,v2),e);
                        edgesMade->insert(qMakePair(v2,v),e);
                        eLi->append(e); //cout<<"157- "; printEdge(e);
                        v->addEdge(e);
                        v2->addEdge(e);
                    }
                }
            }
        }
    }


    //first, remove all collinear points
    vert2d* vs;
    vert2d* ve;
    vert2d* vnew;

    vs=vseed;
    //cout<<"vs: ";vs->print();

    if(vs->e1->v1==vs)
        ve = (vs->e1)->v2;
    else
        ve = (vs->e1)->v1;
    //cout<<"ve: ";ve->print();

    vnew=ve->otherAdjVert(vs);
    //cout<<"new: ";vnew->print();

    do
    {
        while(ccw(vs,ve,vnew)==0)
        {
            edge2d* e = new edge2d(vs,vnew);
            edgesMade->insert(qMakePair(vs,vnew),e);
            edgesMade->insert(qMakePair(vnew,vs),e);
            eLi->append(e); //cout<<"192- "; printEdge(e);

            edge2d* e1 = edgesMade->value(qMakePair(vs,ve));
            //cout<<"e1: ";printEdge(e1);
            edge2d* e2 = edgesMade->value(qMakePair(ve,vnew));
            //cout<<"e2: ";printEdge(e2);

            vert2d* temp = vnew->otherAdjVert(ve);

            vs->replaceEdge(e1,e);
            vnew->replaceEdge(e2,e);

            //remove e1,e2 from edgesMade and delete
            //cout<<"removing e1:"<<endl;printEdge(e1);
            edgesMade->remove(qMakePair(e1->v1,e1->v2));
            edgesMade->remove(qMakePair(e1->v2,e1->v1));
            edgesMade->remove(qMakePair(e2->v1,e2->v2));
            edgesMade->remove(qMakePair(e2->v2,e2->v1));
            eLi->removeAll(e1);
            eLi->removeAll(e2);
            delete e1;
            delete e2;

            //remove ve from vLi, vertsMade and delete
            vLi->removeAll(ve);
            vertsMade->remove(qMakePair(ve->x,ve->y));
            delete ve;


            ve=vnew;
            vnew = temp;
        }
        vs = ve;
        ve = vnew;
        vnew = ve->otherAdjVert(vs);
    }while(vs!=vseed);

    //get the leftmost vert2d
    assert(vertsMade->contains(qMakePair(vseed->x,vseed->y)));
    vert2d* pve1;
    vert2d* pve2;

    //at this point vs--->ve should be in the CLOCKWISE direction
    //one simplification pass:

        //get four points to work on, and next one to start from
            //if one of the points other than the first was vseed, set done flag
            //arrange pointers, etc for our four points

   for(int i=0;i<iterations;i++)
   {
       cout<<"beginning. edgecount: "<<eLi->size()<<endl;

       if(vseed->e1->v1 == vseed)
           pve1=vseed->e1->v2;
       else
           pve1=vseed->e1->v1;
       pve2 = vseed->otherAdjVert(pve1);

       vs=vseed;
       if(pve1->y<pve2->y)
           ve=pve1;
       else if(pve1->y>pve2->y)
           ve=pve2;
       else
           cout<<"that thing that you said would never happen happened"<<endl;

    vert2d* A = vseed;
    vert2d* B = ve;
    vert2d* C = B->otherAdjVert(A);
    vert2d* D = C->otherAdjVert(B);
    vert2d* nextA = D->otherAdjVert(C);


    bool passDoneFlag=false;
    while(!passDoneFlag)
    {
    //we have our four points A-->B-->C-->D
        if(ccw(A,B,C)==1)
        {//angle1 convex
            if(ccw(B,C,D)==1)
            {//both angles convex /'''\
            //ultimately remove B and C, first project AB and CD to intersection
                int isecx, isecy;
                if(B->x==A->x)
                {
                 int mCD = (D->y - C->y)/(D->x - C->x);
                 int bCD = C->y - mCD*C->x;
                 isecx = B->x;
                 isecy = mCD*B->x+bCD;
                }
                else if(D->x==C->x)
                {
                 int mAB = (B->y - A->y)/(B->x - A->x);
                 int bAB = A->y - mAB*A->x;
                 isecx = D->x;
                 isecy = mAB*D->x+bAB;
                }
                else
                {
                int mAB = (B->y - A->y)/(B->x - A->x);
                int mCD = (D->y - C->y)/(D->x - C->x);
                int bAB = A->y - mAB*A->x;
                int bCD = C->y - mCD*C->x;
                if(mAB-mCD==0)
                {//not sure when this stupid case happens
                    isecx =-100;
                    isecy =-100;
                }
                else
                {
                isecx = (bCD - bAB)/(mAB - mCD);
                isecy = mAB * isecx + bAB;
                }
                }

                if(isecx>1 && isecx<img.width()-1 && isecy>1 && isecy<img.width()-1)
                {
                vert2d* VN = new vert2d(isecx,isecy);
                vertsMade->insert(qMakePair(isecx,isecy),VN);
                vLi->append(VN);

                edge2d* eAVN = new edge2d(A,VN);
                edgesMade->insert(qMakePair(A,VN),eAVN);
                edgesMade->insert(qMakePair(VN,A),eAVN);
                eLi->append(eAVN);// cout<<"317- "; printEdge(eAVN);
                edge2d* eVND = new edge2d(VN,D);
                edgesMade->insert(qMakePair(D,VN),eVND);
                edgesMade->insert(qMakePair(VN,D),eVND);
                eLi->append(eVND); //cout<<"321- "; printEdge(eVND);

                VN->addEdge(eAVN);
                VN->addEdge(eVND);

                edge2d* AB = edgesMade->value(qMakePair(A,B));
                edge2d* CD = edgesMade->value(qMakePair(C,D));
                A->replaceEdge(AB,eAVN);
                D->replaceEdge(CD,eVND);

                edgesMade->remove(qMakePair(A,B));
                edgesMade->remove(qMakePair(B,A));
                edgesMade->remove(qMakePair(C,D));
                edgesMade->remove(qMakePair(D,C));
                eLi->removeAll(AB);
                eLi->removeAll(CD);
                vLi->removeAll(B);
                vLi->removeAll(C);
                vertsMade->remove(qMakePair(B->x,B->y));
                vertsMade->remove(qMakePair(C->x,C->y));
                delete AB;
                delete CD;
                delete B;
                delete C;
                }



            }
            else
            {// convex-concave  /\/
             //remove C
             edge2d* e = new edge2d(B,D);
             edgesMade->insert(qMakePair(B,D),e);
             edgesMade->insert(qMakePair(D,B),e);
             eLi->append(e); //cout<<"356- "; printEdge(e);
             edge2d* oldBe = edgesMade->value(qMakePair(B,C));
             edge2d* oldDe = edgesMade->value(qMakePair(C,D));
             B->replaceEdge(oldBe,e);
             D->replaceEdge(oldDe,e);
             edgesMade->remove(qMakePair(B,C));
             edgesMade->remove(qMakePair(C,B));
             edgesMade->remove(qMakePair(C,D));
             edgesMade->remove(qMakePair(D,C));
             eLi->removeAll(oldBe);
             eLi->removeAll(oldDe);
             delete oldBe;
             delete oldDe;
             vLi->removeAll(C);
             vertsMade->remove(qMakePair(C->x,C->y));
             delete C;
            }

        }
        else if(ccw(A,B,C)==-1)
        {
            if(ccw(B,C,D)==0)
            {
                edge2d* e = new edge2d(B,D);
                edgesMade->insert(qMakePair(B,D),e);
                edgesMade->insert(qMakePair(D,B),e);
                eLi->append(e); //cout<<"382- "; printEdge(e);
                edge2d* oldBe = edgesMade->value(qMakePair(B,C));
                edge2d* oldDe = edgesMade->value(qMakePair(C,D));
                B->replaceEdge(oldBe,e);
                D->replaceEdge(oldDe,e);
                edgesMade->remove(qMakePair(B,C));
                edgesMade->remove(qMakePair(C,B));
                edgesMade->remove(qMakePair(C,D));
                edgesMade->remove(qMakePair(D,C));
                eLi->removeAll(oldBe);
                eLi->removeAll(oldDe);
                delete oldBe;
                delete oldDe;
                vLi->removeAll(C);
                vertsMade->remove(qMakePair(C->x,C->y));
                delete C;
            }
            else
            {
             //concave-convex \/\  OR
             //both angles concave \__/
             //either way: remove B
                edge2d* e = new edge2d(A,C);
                edgesMade->insert(qMakePair(A,C),e);
                edgesMade->insert(qMakePair(C,A),e);
                eLi->append(e);// cout<<"407- "; printEdge(e);
                edge2d* oldAe = edgesMade->value(qMakePair(A,B));
                edge2d* oldCe = edgesMade->value(qMakePair(B,C));
                A->replaceEdge(oldAe,e);
                C->replaceEdge(oldCe,e);
                edgesMade->remove(qMakePair(A,B));
                edgesMade->remove(qMakePair(B,A));
                edgesMade->remove(qMakePair(B,C));
                edgesMade->remove(qMakePair(C,B));
                eLi->removeAll(oldAe);
                eLi->removeAll(oldCe);
                delete oldAe;
                delete oldCe;
                vLi->removeAll(B);
                vertsMade->remove(qMakePair(B->x,B->y));
                delete B;
            }
        }
        else
        {//a,b,c straight.  remove b
            edge2d* e = new edge2d(A,C);
            edgesMade->insert(qMakePair(A,C),e);
            edgesMade->insert(qMakePair(C,A),e);
            eLi->append(e); //cout<<"430- "; printEdge(e);
            edge2d* oldAe = edgesMade->value(qMakePair(A,B));
            edge2d* oldCe = edgesMade->value(qMakePair(B,C));
            A->replaceEdge(oldAe,e);
            C->replaceEdge(oldCe,e);
            edgesMade->remove(qMakePair(A,B));
            edgesMade->remove(qMakePair(B,A));
            edgesMade->remove(qMakePair(B,C));
            edgesMade->remove(qMakePair(C,B));
            eLi->removeAll(oldAe);
            eLi->removeAll(oldCe);
            delete oldAe;
            delete oldCe;
            vLi->removeAll(B);
            vertsMade->remove(qMakePair(B->x,B->y));
            delete B;
        }
    //advance!
        //cout<<"advancing"<<endl;
        A = nextA;
        B = nextA->otherAdjVert(D);
        C = B->otherAdjVert(A);
        D = C->otherAdjVert(B);
        nextA = D->otherAdjVert(C);
        //cout<<"done advancing"<<endl;
        if(A==vseed || B==vseed || C==vseed || D==vseed)
        {
            passDoneFlag=true;
        }
    }
}


   delete eLi;
   eLi = new QList<edge2d*>();

   QSet<edge2d*>* visitedE = new QSet<edge2d*>();
   for(int i=0; i<vLi->size(); i++)
   {
       vert2d* vv = vLi->at(i);
       if(!visitedE->contains(vv->e1))
       {
           eLi->append(vv->e1);
           visitedE->insert(vv->e1);
       }
       if(!visitedE->contains(vv->e2))
       {
           eLi->append(vv->e2);
           visitedE->insert(vv->e2);
       }
   }




   polyHull* pHull = new polyHull(vLi,eLi);

   //pHull->print();

   pHull->imgh = h;
   pHull->imgw = w;


   return pHull;
}

//OKAY LETS GET REAL HERE
//for now lets just try to create a single patch
QList<LappedPatch*>* LappedUtils::generatePatches(GLMmodel* model, polyHull* polyhull)
{

//*************REFERENCE PSEUDOCODE*************************************
//PREPROCESSING:
    //Convert entire model into our format.  I know this is goofy but its only done once and makes it easier.

    //For each triangle:
        //make a PatchTri
            //populate vertices - make sure not a repeat by hashing
                //add this triangle to the vert ;)
            //populate edges - make sure not a repeat by hashing
                //add this triangle to the edge ;)


//choose seed triangle, get center point
    //for now, tangent equals vo--->v1  (should be curvature!)
//assign UVS to PatchTri s.t. centered at .5,.5 and tangent is aligned with texture //HOW???? NEED BITANGENT OR NO?  HOW TO DETERMINE SCALE?
    //enqueue edges of triangle

//make our patchlist to return!

//WHILE MESH IS NOT COVERED:

//make visitedSets of edges, verts, tris
//mark everything from seed as visited

//while Q not empty
    //dequeue edge e
        //find other triangle in edge
        //if tri is not in patch
            //if e isects hull  //see struct polyHull::isectUV
                //if homeomorphic to a disc (new vert not in patch OR only 1 edge not in patch)
                  //if new vert already in patch
                      //just add triangle to patch!
                          //enqueue new edges of triangle
                          //and mark them as visited, as well as new tri itself
                  //else need to add vert:
                      //estimate parametrization of new vert like so:
                          //for each triangle containing vert and an already-mapped edge
                              //stick a similar triangle onto that edge, and see where the corresponding location of vert is
                          //take average of those guys as UV for new vert
                      //then add triangle like above: enq new edges, mark any new edges, vert, tri as visited
   //patch is done, delete visitedSets, continue to next patch

    //mesh is covered (or we stopped making patches for whatever reason

    //return patch list!

    //*****************OKAY LETS ACTUALLY IMPLEMENT IT******************************
    QList<LappedPatch*>* PatchList = new QList<LappedPatch*>();
    //PREPROCESSING
    GLMtriangle* glmtris = model->triangles;
    GLfloat* glmverts = model->vertices;
    QHash<GLuint,PatchVert*>* vertsMade = new QHash<GLuint,PatchVert*>();
    QHash<QPair<PatchVert*,PatchVert*>,PatchEdge*>* edgesMade = new QHash<QPair<PatchVert*,PatchVert*>,PatchEdge*>();
    PatchTri** trisMade = new PatchTri*[model->numtriangles];

    for(int i=0;i<model->numtriangles;i++)
    {
        PatchTri* pt = new PatchTri();
        trisMade[i] = pt;
        pt->GLMtri = &glmtris[i];

        //make/find PatchVerts for new triangle
        for(int vi=0;vi<3;vi++)
        {
            PatchVert* v;
            if(vertsMade->contains(pt->GLMtri->vindices[vi]))
            {
                 v = vertsMade->value(pt->GLMtri->vindices[vi]);
                v->tris->append(pt);
            }
            else
            {
                v = new PatchVert();
                v->GLMidx = pt->GLMtri->vindices[vi];
                v->pos.x = glmverts[v->GLMidx];
                v->pos.y = glmverts[v->GLMidx+1];
                v->pos.z = glmverts[v->GLMidx+2];
                v->tris->append(pt);
            }
            switch(vi)
            {
            case 0:
                pt->v0=v;
                break;
            case 1:
                pt->v1=v;
                break;
            case  2:
                pt->v2=v;
                break;
            }
        }

        QPair<PatchVert*,PatchVert*> e01, e12, e20;
        e01 = qMakePair(pt->v0,pt->v1);
        e12 = qMakePair(pt->v1,pt->v2);
        e20 = qMakePair(pt->v2,pt->v0);
        //should triangles have pointers to their edges?
        if(edgesMade->contains(e01))
        {
            PatchEdge* e = edgesMade->value(e01);
            e->addTri(pt);
            pt->e01 = e;
        }
        else
        {
            PatchEdge* e = new PatchEdge(pt->v0, pt->v1);
            edgesMade->insert(e01,e);
            edgesMade->insert(qMakePair(pt->v1,pt->v0),e);
            e->addTri(pt);
            pt->e01=e;
        }
        if(edgesMade->contains(e12))
        {
            PatchEdge* e = edgesMade->value(e12);
            e->addTri(pt);
            pt->e12=e;
        }
        else
        {
            PatchEdge* e = new PatchEdge(pt->v1, pt->v2);
            edgesMade->insert(e12,e);
            edgesMade->insert(qMakePair(pt->v2,pt->v1),e);
            e->addTri(pt);
            pt->e12=e;
        }
        if(edgesMade->contains(e20))
        {
            PatchEdge* e = edgesMade->value(e20);
            e->addTri(pt);
            pt->e20=e;
        }
        else
        {
            PatchEdge* e = new PatchEdge(pt->v2, pt->v0);
            edgesMade->insert(e20,e);
            edgesMade->insert(qMakePair(pt->v0,pt->v2),e);
            e->addTri(pt);
            pt->e20=e;
        }
    }//OK TRIANGLES MADE


    //WHILE MESH IS NOT COVERED******
    for(int wtf=0; wtf<1; wtf++)//***
    {                           //***
    //WHILE MESH IS NOT COVERED******


    //choose seed triangle, get center point
        //for now, tangent equals vo--->v1  (should be curvature!)
    //assign UVS to PatchTri s.t. centered at .5,.5 and tangent is aligned with texture //HOW???? NEED BITANGENT OR NO?  HOW TO DETERMINE SCALE?
        //enqueue edges of triangle

   //DONT DELETE THESE THINGS!
    //hash table of uvs for this specific patch!
    QHash<PatchVert*,vec2<float> >* UVs = new QHash<PatchVert*, vec2<float> >();
    //list of Tris for this specific patch!
    QList<PatchTri*> PTris = *(new QList<PatchTri*>());
    //seed for this specific patch!
    PatchTri* seed = trisMade[rand()%model->numtriangles];

    vec2<float> seedUV0, seedUV1, seedUV2;
    assignSeedUV(seed, seedUV0, seedUV1, seedUV2);
    UVs->insert(seed->v0,seedUV0);
    UVs->insert(seed->v1,seedUV1);
    UVs->insert(seed->v2,seedUV2);

    QQueue<PatchEdge*> edgeQ = *(new QQueue<PatchEdge*>());
    edgeQ.enqueue(seed->e01);
    edgeQ.enqueue(seed->e12);
    edgeQ.enqueue(seed->e20);

    //make visitedSets of edges, verts, tris
    //mark everything from seed as visited
    QSet<PatchVert*> vertsInPatch = *(new QSet<PatchVert*>());
    QSet<PatchEdge*> edgesInPatch = *(new QSet<PatchEdge*>());
    QSet<PatchTri*> trisInPatch = *(new QSet<PatchTri*>());
    vertsInPatch.insert(seed->v0);
    vertsInPatch.insert(seed->v1);
    vertsInPatch.insert(seed->v2);
    edgesInPatch.insert(seed->e01);
    edgesInPatch.insert(seed->e12);
    edgesInPatch.insert(seed->e20);
    trisInPatch.insert(seed);
    PTris.append(seed);

    //while Q not empty
        //dequeue edge e
            //find other triangle in edge
            //if tri is not in patch
                //if e isects hull  //see struct polyHull::isectUV
                    //if homeomorphic to a disc (new vert not in patch OR only 1 edge not in patch)
                      //if new vert already in patch
                          //just add triangle to patch!
                              //enqueue new edges of triangle
                              //and mark them as visited, as well as new tri itself
                      //else need to add vert:
                          //estimate parametrization of new vert like so:
                              //for each triangle containing vert and an already-mapped edge
                                  //stick a similar triangle onto that edge, and see where the corresponding location of vert is
                              //take average of those guys as UV for new vert
                          //then add triangle like above: enq new edges, mark any new edges, vert, tri as visited
       //patch is done, delete visitedSets, continue to next patch

    while(!edgeQ.isEmpty())
    {
        PatchEdge* e = edgeQ.dequeue();
            PatchTri* otherTri;
            if(trisInPatch.contains(e->t1))
            {
                if(e->ntris==2 && !trisInPatch.contains(e->t2))
                {
                    otherTri = e->t2;
                }
                else
                {
                 cout<<"both tris already in patch or edge only has one tri"<<endl;
                }
            }
            else if(e->ntris>=1)
            {
               otherTri = e->t1;
            }
            else
            {
                cout<<"edge somehow has no triangles"<<endl;
                continue;
            }

            //if e isects hull
            vec2<float> uv1 = UVs->value(e->v0);
            vec2<float> uv2 = UVs->value(e->v1);
            if(polyhull->isectHullUV(uv1.x,uv1.y,uv2.x,uv2.y))
            {
                PatchEdge* newEdge;//only relevant if theres only one
                PatchVert* newvert = otherTri->otherVert(e->v0, e->v1);
                //if homeomorphic to a disk
                bool newvertinpatch = vertsInPatch.contains(newvert);
                int numNewEdgesInPatch=0;
                if(edgesInPatch.contains(otherTri->e01)){
                    numNewEdgesInPatch++;
                    newEdge=otherTri->e01;}
                if(edgesInPatch.contains(otherTri->e12)){
                    numNewEdgesInPatch++;
                    newEdge=otherTri->e12;}
                if(edgesInPatch.contains(otherTri->e20)){
                    numNewEdgesInPatch++;
                    newEdge=otherTri->e20;}
                if( !newvertinpatch || numNewEdgesInPatch==1)
                {
                    if(newvertinpatch)
                    {//just add the triangle/edges!
                        edgesInPatch.insert(newEdge);
                        trisInPatch.insert(otherTri);
                        PTris.append(otherTri);
                        edgeQ.enqueue(newEdge);
                    }
                    else
                    {//gotta add the new vertex
                        //estimate parametrization of new vert like so:
                            //for each triangle containing vert and an already-mapped edge
                                //stick a similar triangle onto that edge, and see where the corresponding location of vert is
                            //take average of those guys as UV for new vert
                        //then add triangle like above: enq new edges, mark any new edges, vert, tri as visited

                        vec2<float> sumUVs;
                        int numUVs=0;

                        for(int ti=0; ti<newvert->tris->size(); ti++)
                        {
                            PatchTri* curtri = newvert->tris->at(ti);
                            PatchEdge* mappedE;
                            PatchVert *_A,*_B;//args to estimateUV
                            if(edgesInPatch.contains(curtri->e01))
                            {
                                assert(newvert == curtri->v2);
                                mappedE = curtri->e01;
                                _A = curtri->v0;
                                _B = curtri->v1;
                            }
                            else if(edgesInPatch.contains(curtri->e12))
                            {
                                assert(newvert == curtri->v0);
                                mappedE = curtri->e12;
                                _A = curtri->v1;
                                _B = curtri->v2;
                            }
                            else if(edgesInPatch.contains(curtri->e20))
                            {
                                assert(newvert == curtri->v1);
                                mappedE = curtri->e20;
                                _A = curtri->v2;
                                _B = curtri->v0;
                            }
                            else
                                continue;
                            numUVs++;
                            sumUVs = sumUVs + estimateUV(_A,_B,newvert, UVs->value(_A), UVs->value(_B));
                        }
                        //sumUVs is actually now the UVs we want
                        sumUVs = sumUVs / (1.0+numUVs);
                        //add all the crap
                        UVs->insert(newvert,sumUVs);
                        vertsInPatch.insert(newvert);
                        trisInPatch.insert(otherTri);
                        if(!edgesInPatch.contains(otherTri->e01))
                            {edgesInPatch.insert(otherTri->e01);
                             edgeQ.enqueue(otherTri->e01);}
                        if(!edgesInPatch.contains(otherTri->e12))
                            {edgesInPatch.insert(otherTri->e12);
                             edgeQ.enqueue(otherTri->e12);}
                        if(!edgesInPatch.contains(otherTri->e20))
                        {edgesInPatch.insert(otherTri->e20);
                         edgeQ.enqueue(otherTri->e20);}
                        PTris.append(otherTri);
                    }
                }
            }//endif isectHull
    }//endwhile Q!Empty
    //PTris, UVs, seed should be sufficient to define patch!
    LappedPatch* newPatch = new LappedPatch();
    newPatch->tris = &PTris;
    newPatch->seed = seed;
    newPatch->uvs = UVs;
    PatchList->append(newPatch);


    //delete temporary stupid stuff
    delete &edgeQ;
    delete &edgesInPatch;
    delete &vertsInPatch;
    delete &trisInPatch;


    //ENDWHILE MESH NOT COVERED*
    }                       //**
    //ENDWHILE MESH NOT COvERED*

    return PatchList;
}

void LappedUtils::vizualizePatch(LappedPatch* patch, QImage* img)
{
    QPainter patr(img);
    patr.setPen(Qt::black);
    QHash<PatchVert*, vec2<float> >* UVs = patch->uvs;

    for(int i=0;i<patch->tris->size();i++)
    {
        PatchTri* pt = patch->tris->at(i);
        drawEdgeFromUV(img, &patr, UVs->value(pt->v0), UVs->value(pt->v1));
        drawEdgeFromUV(img, &patr, UVs->value(pt->v1), UVs->value(pt->v2));
        drawEdgeFromUV(img, &patr, UVs->value(pt->v2), UVs->value(pt->v0));
    }
    //draw seed in red because why not
    patr.setPen(Qt::red);
    PatchTri* pt = patch->seed;
    drawEdgeFromUV(img, &patr, UVs->value(pt->v0), UVs->value(pt->v1));
    drawEdgeFromUV(img, &patr, UVs->value(pt->v1), UVs->value(pt->v2));
    drawEdgeFromUV(img, &patr, UVs->value(pt->v2), UVs->value(pt->v0));
    patr.end();
}
void LappedUtils::drawEdgeFromUV(QImage* img, QPainter* patr, vec2<float> v0, vec2<float>v1)
{
    int x0,y0,x1,y1;
    x0 = v0.x * img->width();
    y0 = (1.0-v0.y) * img->height();
    x1 = v1.x * img->width();
    y1 = (1.0-v1.y) * img->height();
    patr->drawLine(x0,y0,x1,y1);
}
