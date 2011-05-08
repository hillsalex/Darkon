#include "modelshot.h"
#include "drawengine.h"
#include "meshoperator.h"
#include <time.h>
using std::cout;
using std::endl;

QString modelstring = "mug";


modelShot::modelShot(DrawEngine* parent,QHash<QString, QGLShaderProgram *>* shad, QHash<QString, GLuint>* tex, QHash<QString, Model>* mod) : Shot(parent,shad,tex,mod)
{
    //lasts 150 frames
    m_lifespan = 150;
    m_operator =  new MeshOperator();
}

modelShot::~modelShot()
{

}

void modelShot::renderNormal(const Vector4 &vertex, const Vector4 &direction)
{
    Vector4 normalDirection = direction.getNormalized();

    // Draw a normal with a fixed length of 0.15
    glBegin(GL_LINES);
    glVertex3dv(vertex.data);
    glVertex3dv((vertex + normalDirection * 0.05).data);
    glEnd();

    // End the normal with an axis-aligned billboarded triangle (billboarding means always rotating
    // to face the camera, and axis-aligned means it can only rotate around the axis of the normal)
    Vector4 eye;
    eye.x=m_engine->camera_.eye.x;
    eye.y=m_engine->camera_.eye.y;
    eye.z=m_engine->camera_.eye.z;
    Vector4 triangleVector = direction.cross(eye - vertex);
    if (triangleVector.getMagnitude2() > 1.0e-6f)
    {
        triangleVector = triangleVector.getNormalized() * 0.01;
        glBegin(GL_TRIANGLES);
        glVertex3dv((vertex + normalDirection * 0.1).data);
        glVertex3dv((vertex + normalDirection * 0.05 - triangleVector).data);
        glVertex3dv((vertex + normalDirection * 0.05 + triangleVector).data);
        glEnd();
    }
}




//In begin, initialize things that could not have been initialized beforehand
//(gl state)
void modelShot::begin()
{
    GLMmodel* model = models_->value(modelstring).model;
    m_operator->calculateCurvatures(model);
/*

    QList<QString> keys = models_->keys();
    foreach (QString k,keys)
    {
        GLMmodel* model = models_->value(k).model;
        m_operator->calculateCurvatures(model);
    }*/
}

//called every frame before draw.
void modelShot::update()
{
    m_framesElapsed++;
    if(m_framesElapsed>m_lifespan)
    {
        //end shot and call engine
    }
}

//draw!
extern "C"{
    extern void APIENTRY glActiveTexture(GLenum);
}
void modelShot::draw()
{
    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glMatrixMode(GL_MODELVIEW);
    glActiveTexture(GL_TEXTURE0);
    QList<QString> keys = models_->keys();
    bool drawVertCurvatures = false;
    bool drawVertMinCurvatures = false;
    bool drawFaceCurvatures = true;
    float scaleFactor = 1/20.0;




    GLMmodel* model = models_->value(modelstring).model;
    GLfloat* curvatures = model->vertMaxCurvatures;
    GLfloat* mincurvatures = model->vertMinCurvatures;
    GLfloat* vertices = model->vertices;
    GLMtriangle* triangles = model->triangles;
    GLfloat* triCurvatures = model->triCurvatures;
    if (drawVertCurvatures)
    {
        glBegin(GL_LINES);
        for (int i=0;i<model->numvertices;i++)
        {
            Vector4 pos;
            pos.x = vertices[i*3];
            pos.y = vertices[i*3+1];
            pos.z = vertices[i*3+2];
            Vector4 dir;
            dir.x = curvatures[i*3];
            dir.y = curvatures[i*3+1];
            dir.z = curvatures[i*3+2];
            glColor3f(1.0,0,0);
            renderNormal(pos,dir);
            glColor3f(1,1,1);
            //glVertex3f(vertices[i*3], vertices[i*3+1],vertices[i*3+2]); // origin of the line
            //glVertex3f(vertices[i*3]+curvatures[i*3]*scaleFactor, vertices[i*3+1]+curvatures[i*3+1]*scaleFactor,vertices[i*3+2]+curvatures[i*3+2]*scaleFactor); // ending point of the line
        }
        glEnd( );
    }
    if (drawVertMinCurvatures)
    {
        glBegin(GL_LINES);
        for (int i=0;i<model->numvertices;i++)
        {
            Vector4 pos;
            pos.x = vertices[i*3];
            pos.y = vertices[i*3+1];
            pos.z = vertices[i*3+2];
            Vector4 dir;
            dir.x = mincurvatures[i*3];
            dir.y = mincurvatures[i*3+1];
            dir.z = mincurvatures[i*3+2];
            glColor3f(0,0,1.0);
            renderNormal(pos,dir);
            glColor3f(1,1,1);
            //glVertex3f(vertices[i*3], vertices[i*3+1],vertices[i*3+2]); // origin of the line
            //glVertex3f(vertices[i*3]+curvatures[i*3]*scaleFactor, vertices[i*3+1]+curvatures[i*3+1]*scaleFactor,vertices[i*3+2]+curvatures[i*3+2]*scaleFactor); // ending point of the line
        }
        glEnd( );
    }

    if (drawFaceCurvatures)
    {
        glBegin(GL_LINES);
        for (int i=0;i<model->numtriangles;i++)
        {

                float x=0;
                float y=0;
                float z=0;
                x+= vertices[triangles[i].vindices[0]*3];
                y+= vertices[triangles[i].vindices[0]*3+1];
                z+= vertices[triangles[i].vindices[0]*3+2];

                x+= vertices[triangles[i].vindices[1]*3];
                y+= vertices[triangles[i].vindices[1]*3+1];
                z+= vertices[triangles[i].vindices[1]*3+2];

                x+=vertices[triangles[i].vindices[2]*3];
                y+= vertices[triangles[i].vindices[2]*3+1];
                z+= vertices[triangles[i].vindices[2]*3+2];

                //cout << "Xverts: " << vertices[triangles[i].vindices[0]*3] <<
                //        ',' << vertices[triangles[i].vindices[1]*3] << ','
                //        << vertices[triangles[i].vindices[2]*3] << endl;

                x=x/3.0;
                y=y/3.0;
                z=z/3.0;
                //cout << "x: " << x << endl;

                Vector4 pos;
                pos.x = x;
                pos.y = y;
                pos.z = z;
                Vector4 dir;
                dir.x = triCurvatures[i*3];
                dir.y = triCurvatures[i*3+1];
                dir.z = triCurvatures[i*3+2];
                glColor3f(0,1,0);
                renderNormal(pos,dir);
                glColor3f(1,1,1);

                //glVertex3f(x, y,z); // origin of the line
                //glVertex3f(x+triCurvatures[i*3]*scaleFactor, y+triCurvatures[i*3+1]*scaleFactor,z+triCurvatures[i*3+2]*scaleFactor); // origin of the line
                //glVertex3f(x+curve[0]/20.0, y+curve[1]/20.0,z+curve[2]/20.0); // ending point of the line
                //cout << curve[0] << ',' << curve[1] << ',' << curve[2] << endl;

        }
        glEnd( );
    }


    glBindTexture(GL_TEXTURE_CUBE_MAP, textures_->value("cube_map_1"));
    shader_programs_->value(NAIL_SHADER)->bind();
    shader_programs_->value(NAIL_SHADER)->setUniformValue("CubeMap",0);
    shader_programs_->value(NAIL_SHADER)->setUniformValue("eyept",m_engine->camera_.eye.x, m_engine->camera_.eye.y, m_engine->camera_.eye.z);
    glCallList(models_->value(modelstring).idx);
    shader_programs_->value("reflect")->release();
}


