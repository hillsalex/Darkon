/**
  A simple OpenGL drawing engine.

  @author psastras
**/

#include "drawengine.h"

#ifdef WIN32
#define GL_GLEXT_PROTOTYPES
#define GL_GLEXT_LEGACY
#endif

#include <glm/glm.h>
#include <qgl.h>
#include <QKeyEvent>
#include <QGLContext>
#include <QHash>

#ifdef WIN32
#include "glext.h"
#endif

#include <QGLShaderProgram>
#include <QQuaternion>
#include <QVector3D>
#include <QString>
#include <GL/glu.h>
#include <iostream>
#include <QFile>
#include <QGLFramebufferObject>
#include <GL/glext.h>
#include <QFileDialog>

//***********************************
//INCLUDE SHOTS HERE AS YOU MAKE THEM
#include <testShot.h>
#include <MPSandbox.h>
#include <SphereShot.h>
//***********************************

using std::cout;
using std::endl;

extern "C"{
    extern void APIENTRY glActiveTexture (GLenum);
}

/**
  @paragraph DrawEngine ctor.  Expects a Valid OpenGL context and the viewport's current
  width and height.  Initializes the draw engine.  Loads models,textures,shaders,
  and allocates framebuffers.  Also sets up OpenGL to begin drawing.

  @param context The current OpenGL context this drawing engine is associated with.
  Probably should be the context from the QGLWidget.

  @param w The viewport width used to allocate void drawUnitRope();the correct framebuffer size.
  @param h The viewport heigh used to alloacte the correct framebuffer size.

**/
DrawEngine::DrawEngine(const QGLContext *context,int w,int h, GLWidget* widget) : context_(context) {
    //initialize ogl settings
    glEnable(GL_TEXTURE_2D);
    glFrontFace(GL_CCW);
    glDisable(GL_DITHER);
    glDisable(GL_LIGHTING);
    glShadeModel(GL_FLAT);
    glClearColor(0.0f,0.0f,0.0f,0.0f);

    //init member variables
    previous_time_ = 0.0f;
    camera_.center.x = 0.f,camera_.center.y = 0.f,camera_.center.z = 0.f;
    camera_.eye.x = 0.f,camera_.eye.y = 0.0f,camera_.eye.z = 2.f;
    camera_.up.x = 0.f,camera_.up.y = 1.f,camera_.up.z = 0.f;
    camera_.near = 0.1f,camera_.far = 100.f;
    camera_.fovy = 60.f;


    load_shaders();
    load_models();
    load_textures();
    create_curvatures();
    create_fbos(w,h);


    /****************************************/
    //  SHOT ARRAY
    m_shots = new QList<Shot*>();
    m_curShot =NULL;

    frameNumber = 0;
    m_widget = widget;

    //m_shots->append(  new testShot(this, &shader_programs_, &textures_, &models_));
    m_shots->append(  new SphereShot(this, &shader_programs_, &textures_, &models_));
    //m_shots->append(  new MPSandbox(this, &shader_programs_, &textures_, &models_));

    m_shots->at(m_curShot)->begin();
    /****************************************/
}

/**
  @paragraph Dtor
**/
DrawEngine::~DrawEngine() {
    foreach(QGLShaderProgram *sp,shader_programs_)
        delete sp;
    foreach(QGLFramebufferObject *fbo,framebuffer_objects_)
        delete fbo;
    foreach(GLuint id,textures_)
        ((QGLContext *)(context_))->deleteTexture(id);
    foreach(Model m,models_)
        glmDelete(m.model);

    for(int i=0;i<m_shots->size();i++)
    {
        delete m_shots->at(i);
    }
    delete m_shots;

}

/**
  @paragraph Loads models used by the program.  Called by the ctor once upon
  initialization.
**/
void DrawEngine::load_models() {

    models_["nail"].model = glmReadOBJ(  "../Hatching/src/models/nail.obj"  );
    glmUnitize(models_["nail"].model);
    models_["nail"].idx = glmList(models_["nail"].model,GLM_SMOOTH);


}
/**
  @paragraph Loads shaders used by the program.  Caleed by the ctor once upon
  initialization.
**/


void DrawEngine::load_shaders() {
    shader_programs_["reflect"] = new QGLShaderProgram(context_);
    shader_programs_["reflect"]->addShaderFromSourceFile(QGLShader::Vertex,"../Hatching/src/shaders/reflect.vert");
    shader_programs_["reflect"]->addShaderFromSourceFile(QGLShader::Fragment,"../Hatching/src/shaders/reflect.frag");
   shader_programs_["reflect"]->link();
}
/**
  @paragraph Loads textures used by the program.  Called by the ctor once upon
  initialization.
**/
void DrawEngine::load_textures() {
   // cout << "Loading textures..." << endl;
    /** CUBE MAP EXAMPLE **/
    QList<QFile *> fileList;
    fileList.append(new QFile("../Hatching/src/textures/astra/posx.jpg"));
    fileList.append(new QFile("../Hatching/src/textures/astra/negx.jpg"));
    fileList.append(new QFile("../Hatching/src/textures/astra/posy.jpg"));
    fileList.append(new QFile("../Hatching/src/textures/astra/negy.jpg"));
    fileList.append(new QFile("../Hatching/src/textures/astra/posz.jpg"));
    fileList.append(new QFile("../Hatching/src/textures/astra/negz.jpg"));
    textures_["cube_map_1"] = load_cube_map(fileList);
    /** NORMAL TEXTURE EXAMPLE **/
    //textures_["rope_occlusion"] = load_texture(":/textures/ropeOcc.png");
    //*********************************************************************


}


void DrawEngine::create_curvatures()
{
    foreach (Model m, models_)
    {
        GLMmodel* model = m.model;
        QHash<GLint,QSet<GLint>* > adjacencyMatrix;
        QHash<GLint, GLint> vertexNormalComps;
        for(int i=0;i<model->numtriangles;i++) //Create adjacency matrix and vertex-normal mapping

        {
            GLMtriangle tri = model->triangles[i];
            GLint v0 = tri.vindices[0];
            GLint v1 = tri.vindices[1];
            GLint v2 = tri.vindices[2];

            if (adjacencyMatrix.contains(v0))
            {
                adjacencyMatrix.value(v0)->insert(v1);
                adjacencyMatrix.value(v0)->insert(v2);
            }
            else
            {
                QSet<GLint>* set = new QSet<GLint>;
                set->insert(v1);
                set->insert(v2);
                adjacencyMatrix.insert(v0,set);
            }

            if (adjacencyMatrix.contains(v1))
            {
                adjacencyMatrix.value(v1)->insert(v0);
                adjacencyMatrix.value(v1)->insert(v2);
            }
            else
            {
                QSet<GLint>* set = new QSet<GLint>;
                set->insert(v0);
                set->insert(v2);
                adjacencyMatrix.insert(v1,set);
            }

            if (adjacencyMatrix.contains(v2))
            {
                adjacencyMatrix.value(v2)->insert(v1);
                adjacencyMatrix.value(v2)->insert(v0);
            }
            else
            {
                QSet<GLint>* set = new QSet<GLint>;
                set->insert(v1);
                set->insert(v0);
                adjacencyMatrix.insert(v2,set);
            }
            if (!vertexNormalComps.contains(v0)) vertexNormalComps.insert(v0,tri.nindices[0]);
            if (!vertexNormalComps.contains(v1)) vertexNormalComps.insert(v1,tri.nindices[1]);
            if (!vertexNormalComps.contains(v2)) vertexNormalComps.insert(v2,tri.nindices[2]);
        }
        QList<GLint> keys = adjacencyMatrix.keys();
        for (int i=0;i<keys.size();i++)
        {
            GLint vertexIndex = keys[i]; //Index of current vertex
            GLint normalIndex = vertexNormalComps.value(vertexIndex); //Normal index of current vertex
            vec3<float> normal; //Normal of current index
            vec3<float> vert; //Position of current vertex
            normal.x = model->normals[normalIndex*3];
            normal.y = model->normals[normalIndex*3+1];
            normal.z = model->normals[normalIndex*3+2];
            vert.x = model->vertices[vertexIndex*3];
            vert.y = model->vertices[vertexIndex*3+1];
            vert.z = model->vertices[vertexIndex*3+2];

            QList<GLint> adjacentVertices = adjacencyMatrix.value(vertexIndex)->values(); //Adjacent vertices
            GLint v0index = adjacentVertices[0]; //first adjacent vertex, for calculating basis
            vec3<float> v0;
            //v0.x = model->vertices[v0*3];
            //v0.y = model->vertices[v0*3+1];
            //v0.z = model->vertices[v0*3+2];
            //Translate v0 to vert as origin, project to uv plane
            vec3<float> u = v0-vert;
            u.normalize();
            u = u - normal * u.dot(normal);
            //v = normal.cross(u);


        }
    }
}
/**
  @paragraph Creates the intial framebuffers for drawing.  Called by the ctor once
  upon initialization.

  @todo Finish filling this in.

  @param w:    the viewport width
  @param h:    the viewport height
**/
void DrawEngine::create_fbos(int w,int h) {
    //Allocate the main framebuffer object for rendering the scene to
    //This needs a depth attachment.
   /* framebuffer_objects_["fbo_0"] = new QGLFramebufferObject(w,h,QGLFramebufferObject::Depth,
                                                             GL_TEXTURE_2D,GL_RGB16F_ARB);*/
}
/**
  @paragraph Reallocates all the framebuffers.  Called when the viewport is
  resized.

  @param w:    the viewport width
  @param h:    the viewport height
**/
void DrawEngine::realloc_framebuffers(int w,int h) {
    foreach(QGLFramebufferObject *fbo,framebuffer_objects_)  {
        const QString &key = framebuffer_objects_.key(fbo);
        QGLFramebufferObjectFormat format = fbo->format();
        delete fbo;
        framebuffer_objects_[key] = new QGLFramebufferObject(w,h,format);
    }
}



/**
  @paragraph Should render one frame at the given elapsed time in the program.
  Assumes that the GL context is valid when this method is called.

  @todo Finish filling this in

  @param time: the current program time in milliseconds
  @param w:    the viewport width
  @param h:    the viewport height

**/
void DrawEngine::draw_frame(float time,int w,int h) {
    fps_ = 1000.f / (time - previous_time_),previous_time_ = time;
    m_w = w;
    m_h = h;
    perspective_camera(w,h);    

    m_shots->at(m_curShot)->update();
    m_shots->at(m_curShot)->draw();
}

void DrawEngine::orthogonal_camera()
{
    int w = m_w;
    int h = m_h;
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0,static_cast<float>(w),static_cast<float>(h),0.f,-1.f,1.f);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}




void DrawEngine::endShot()
{
    /**
      RESET THE FUCKING CAMERA
    **/
    camera_.center.x = 0.f,camera_.center.y = 0.f,camera_.center.z = 0.f;
    camera_.eye.x = 0.f,camera_.eye.y = 0.0f,camera_.eye.z = 2.f;
    camera_.up.x = 0.f,camera_.up.y = 1.f,camera_.up.z = 0.f;
    camera_.near = 0.1f,camera_.far = 100.f;
    camera_.fovy = 60.f;

    if(m_curShot>=m_shots->size()-1)
    {
    //some sort of end action...
    }
    else
    {
        m_curShot++;
        m_shots->at(m_curShot)->begin();
    }
}

/*
 Will be useful so keep around.
 */
/**
  @paragraph Draws a textured quad. The texture most be bound and unbound
  before and after calling this method - this method assumes that the texture
  has been bound before hand.

  @param w: the width of the quad to draw
  @param h: the height of the quad to draw
  @param flip: flip the texture vertically

**/
void DrawEngine::textured_quad(int w,int h,bool flip) {
    glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_CLAMP_TO_EDGE);
    glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_CLAMP_TO_EDGE);
    glBegin(GL_QUADS);
    glTexCoord2f(0.0f,flip ? 1.0f : 0.0f);
    glVertex2f(0.0f,0.0f);
    glTexCoord2f(1.0f,flip ? 1.0f : 0.0f);
    glVertex2f(w,0.0f);
    glTexCoord2f(1.0f,flip ? 0.0f : 1.0f);
    glVertex2f(w,h);
    glTexCoord2f(0.0f,flip ? 0.0f : 1.0f);
    glVertex2f(0.0f,h);
    glEnd();
}

/**
  @paragraph Called to switch to the perspective OpenGL camera.
  Used to render the scene regularly with the current camera parameters.

  @param w: the viewport width
  @param h: the viewport height

**/
void DrawEngine::perspective_camera(int w,int h) {
    float ratio = w / static_cast<float>(h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(camera_.fovy,ratio,camera_.near,camera_.far);
    gluLookAt(camera_.eye.x,camera_.eye.y,camera_.eye.z,
              camera_.center.x,camera_.center.y,camera_.center.z,
              camera_.up.x,camera_.up.y,camera_.up.z);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void DrawEngine::perspective_camera() {
    int w= m_w;
    int h = m_h;
    float ratio = w / static_cast<float>(h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(camera_.fovy,ratio,camera_.near,camera_.far);
    gluLookAt(camera_.eye.x,camera_.eye.y,camera_.eye.z,
              camera_.center.x,camera_.center.y,camera_.center.z,
              camera_.up.x,camera_.up.y,camera_.up.z);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

/**
  @paragraph Called to switch to an orthogonal OpenGL camera.
  Useful for rending a textured quad across the whole screen.

  @param w: the viewport width
  @param h: the viewport height

**/
void DrawEngine::orthogonal_camera(int w,int h) {
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0,static_cast<float>(w),static_cast<float>(h),0.f,-1.f,1.f);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

/**
  @paragraph Called when the viewport has been resized. Needs to
  resize the camera perspective and reallocate the framebuffer
  sizes.

  @param w: the viewport width
  @param h: the viewport height

**/
void DrawEngine::resize_frame(int w,int h) {
    glViewport(0,0,w,h);
    m_w =w;
    m_h = h;
    realloc_framebuffers(w,h);      
}

/**
  @paragraph Called by GLWidget when the mouse is dragged.  Rotates the camera
  based on mouse movement.

  @param p0: the old mouse position
  @param p1: the new mouse position
**/
void DrawEngine::mouse_drag_event(float2 p0,float2 p1) {
    int dx = p1.x - p0.x,dy = p1.y - p0.y;
    QQuaternion qq = QQuaternion::fromAxisAndAngle(0, 1, 0, -dx / 5.0);
    QVector3D qv3 = qq.rotatedVector(QVector3D(camera_.eye.x, camera_.eye.y,
                                               camera_.eye.z));
    qq = QQuaternion::fromAxisAndAngle(qq.rotatedVector(QVector3D(1, 0, 0)), dy / 5.0);
    qv3 = qq.rotatedVector(qv3);
    camera_.eye.x = qv3.x(), camera_.eye.y = qv3.y(), camera_.eye.z = qv3.z();
}

/**
  @paragraph Called by GLWidget when the mouse wheel is turned. Zooms the camera in
  and out.

  @param dx: The delta value of the mouse wheel movement.
**/
void DrawEngine::mouse_wheel_event(int dx) {
    if((camera_.center - camera_.eye).getMagnitude() > .5 || dx < 0)
        camera_.eye += (camera_.center - camera_.eye).getNormalized() * dx * .005;
}

/**
  @paragraph Loads the cube map into video memory.

  @param files: a list of files containing the cube map images (should be length
  six) in order.
  @return The assigned OpenGL id to the cube map.
**/
GLuint DrawEngine::load_cube_map(QList<QFile *> files) {
    GLuint id;
    glGenTextures(1,&id);
    glBindTexture(GL_TEXTURE_CUBE_MAP,id);
    for(unsigned i = 0; i < 6; ++i) {
        QImage image,texture;
        image.load(files[i]->fileName());
        image = image.mirrored(false,true);
        texture = QGLWidget::convertToGLFormat(image);
        texture = texture.scaledToWidth(1024,Qt::SmoothTransformation);
        glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i,3,3,texture.width(),texture.height(),0,GL_RGBA,GL_UNSIGNED_BYTE,texture.bits());
        gluBuild2DMipmaps(GL_TEXTURE_CUBE_MAP_POSITIVE_X +i, 3, texture.width(), texture.height(), GL_RGBA, GL_UNSIGNED_BYTE, texture.bits());
    }
    glTexParameteri(GL_TEXTURE_CUBE_MAP,GL_TEXTURE_MIN_FILTER,GL_NEAREST_MIPMAP_NEAREST);
    glTexParameteri(GL_TEXTURE_CUBE_MAP,GL_TEXTURE_MAG_FILTER,GL_NEAREST_MIPMAP_NEAREST);
    glBindTexture(GL_TEXTURE_CUBE_MAP,0);
    return id;
}


/**
  Loads a texture into video memory.

  Param: a Qfile, texture
**/

GLuint DrawEngine::load_texture(QString name)
{

 QImage image,texture;
 GLuint id;
 image.load(name);
 image = image.mirrored(false,true);
 //image = (*img).mirrored(false,true);
 texture = QGLWidget::convertToGLFormat(image);


 glGenTextures(1,&id);

 glBindTexture(GL_TEXTURE_2D,id);

 gluBuild2DMipmaps(GL_TEXTURE_2D, 3, texture.width(), texture.height(), GL_RGBA, GL_UNSIGNED_BYTE, texture.bits());
 glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
 glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
 glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
 glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

 return id;
}


/**
  @paragraph Called when a key has been pressed in the GLWidget.

  @param event: The key press event associated with the current key press.
  **/
void DrawEngine::key_press_event(QKeyEvent *event) {
    switch(event->key()) {

    }
}

