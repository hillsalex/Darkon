#include "LappedUtils.h"
#include <QQueue>
#include <QHash>

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


polyHull* LappedUtils::getPolyHull(QImage* blob,int iterations)
{
    QImage img = *blob;
    //QImage img;
    //img.load("/home/mprice/Desktop/Patch/PatchMask.png");
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
                        eLi->append(e);
                        //cout<<"adding edge ("<<v->x<<","<<v->y<<")->("<<v2->x<<","<<v2->y<<")"<<endl;
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
            eLi->append(e);
            //cout<<"e: ";printEdge(e);
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

   for(int i=0;i<6;i++)
   {
       cout<<"beginning. vseed: "<<vseed<<endl;

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
    cout<<"here?"<<endl;
    vert2d* C = B->otherAdjVert(A);
    cout<<"yup"<<endl;
    vert2d* D = C->otherAdjVert(B);
    vert2d* nextA = D->otherAdjVert(C);
    cout<<"done beginning"<<endl;


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
                eLi->append(eAVN);
                edge2d* eVND = new edge2d(VN,D);
                edgesMade->insert(qMakePair(D,VN),eVND);
                edgesMade->insert(qMakePair(VN,D),eVND);
                eLi->append(eVND);

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
             eLi->append(e);
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
                eLi->append(e);
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
                eLi->append(e);
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
            eLi->append(e);
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
        {cout<<"setting that flag"<<endl;
            passDoneFlag=true;
        }
    }
}

   //DEBUG RENDER OUTPUT
    /*
    QPainter patr(&img);
    patr.setPen(Qt::black);
    for(int i=0;i<vLi->size();i++)
    {
        vert2d* v = vLi->at(i);
        patr.drawLine(v->e1->v1->x,v->e1->v1->y,v->e1->v2->x,v->e1->v2->y);
        patr.drawLine(v->e2->v1->x,v->e2->v1->y,v->e2->v2->x,v->e2->v2->y);
    }
    patr.end();
    for(int i=0;i<vLi->size();i++)
    {
        vert2d* v = vLi->at(i);
        int _idx = v->y * w + v->x;
        assert(v->nedges==2);
        imgdata[4*_idx]=0;
        imgdata[4*_idx+1]=0;
        imgdata[4*_idx+2]=255;
    }
    img.save("/home/mprice/Desktop/Patch/output.png","PNG");
    */

   polyHull* pHull = new polyHull(vLi,eLi);
   return pHull;
}
