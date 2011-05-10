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
#include <QVector>

//including solver things.
#include "src/common.h"
#include "src/SparseMatrix.h"
#include "src/LinearSolver.h"

//***********************************
//INCLUDE SHOTS HERE AS YOU MAKE THEM
#include <testShot.h>
#include <MPSandbox.h>
#include <SphereShot.h>
#include <modelshot.h>
#include <BYSandbox.h>
//***********************************
#define T(x) (model->triangles[(x)])

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
    //create_curvatures();
    create_fbos(w,h);


    /****************************************/
    //  SHOT ARRAY
    m_shots = new QList<Shot*>();
    m_curShot =NULL;

    frameNumber = 0;
    m_widget = widget;

    //m_shots->append(  new testShot(this, &shader_programs_, &textures_, &models_));
    //m_shots->append(  new modelShot(this, &shader_programs_, &textures_, &models_));
    //m_shots->append(  new SphereShot(this, &shader_programs_, &textures_, &models_));
    m_shots->append(  new MPSandbox(this, &shader_programs_, &textures_, &models_));

    m_shots->at(m_curShot)->begin();
    //BYSandbox *blah = new BYSandbox(this);
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

    models_["buddha"].model = glmReadOBJ(  "../Hatching/src/models/buddha.obj"  );
    glmUnitize(models_["buddha"].model);
    models_["buddha"].idx = glmList(models_["buddha"].model,GLM_SMOOTH);

    models_["teapot"].model = glmReadOBJ(  "../Hatching/src/models/teapot.obj"  );
    glmUnitize(models_["teapot"].model);
    models_["teapot"].idx = glmList(models_["teapot"].model,GLM_SMOOTH);

    models_["ring"].model = glmReadOBJ(  "../Hatching/src/models/ring.obj"  );
    glmUnitize(models_["ring"].model);
    models_["ring"].idx = glmList(models_["ring"].model,GLM_SMOOTH);

    models_["cylinder"].model = glmReadOBJ(  "../Hatching/src/models/Cylinder.obj"  );
    glmUnitize(models_["cylinder"].model);
    models_["cylinder"].idx = glmList(models_["cylinder"].model,GLM_SMOOTH);

    models_["sphere"].model = glmReadOBJ(  "../Hatching/src/models/Sphere.obj"  );
    glmUnitize(models_["sphere"].model);
    models_["sphere"].idx = glmList(models_["sphere"].model,GLM_SMOOTH);

    models_["torus"].model = glmReadOBJ(  "../Hatching/src/models/Torus.obj"  );
    glmUnitize(models_["torus"].model);
    models_["torus"].idx = glmList(models_["torus"].model,GLM_SMOOTH);

    models_["spherehd"].model = glmReadOBJ(  "../Hatching/src/models/Sphere_highpoly.obj"  );
    glmUnitize(models_["spherehd"].model);
    models_["spherehd"].idx = glmList(models_["spherehd"].model,GLM_SMOOTH);

    models_["conehd"].model = glmReadOBJ(  "../Hatching/src/models/Cone2_highpoly.obj"  );
    glmUnitize(models_["conehd"].model);
    models_["conehd"].idx = glmList(models_["conehd"].model,GLM_SMOOTH);

    models_["mug"].model = glmReadOBJ(  "../Hatching/src/models/Mug.obj"  );
    glmUnitize(models_["mug"].model);
    models_["mug"].idx = glmList(models_["mug"].model,GLM_SMOOTH);

    models_["torushd"].model = glmReadOBJ(  "../Hatching/src/models/Torus_highpoly.obj"  );
    glmUnitize(models_["torushd"].model);
    models_["torushd"].idx = glmList(models_["torushd"].model,GLM_SMOOTH);

    models_["cube"].model = glmReadOBJ(  "../Hatching/src/models/cube.obj"  );
    glmUnitize(models_["cube"].model);
    models_["cube"].idx = glmList(models_["cube"].model,GLM_SMOOTH);

    models_["bird"].model = glmReadOBJ(  "../Hatching/src/models/Durbird.obj"  );
    glmUnitize(models_["bird"].model);
    glmFacetNormals(models_["bird"].model);
    glmVertexNormals(models_["bird"].model,90);
    models_["bird"].idx = glmList(models_["bird"].model,GLM_SMOOTH);

    //models_["hourglass"].model = glmReadOBJ(  "../Hatching/src/models/hourglass.obj"  );
    //glmUnitize(models_["hourglass"].model);
    //models_["hourglass"].idx = glmList(models_["hourglass"].model,GLM_SMOOTH);

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

void DrawEngine::orientTexture(GLMmodel* model, QVector<int> patch) {
    //debugging: see if you can get this to be just the number of vertices in the vector -> hashmap, maybe.
    SparseMatrix H00 = SparseMatrix(model->numvertices*2, model->numvertices*2);

    int idx_ax, idx_ay, idx_bx, idx_by, idx_cx, idx_cy;
    double alpha, beta;

    //calculating f0.
    double f0[model->numvertices*2];

    //entering in values for H00 for all vertices.
    for (int tri = 0; tri < patch.size(); tri++) {
        GLMtriangle ctri = T(patch.at(tri));

        //idx locations in H
        idx_ax = ctri.vindices[0]*2;
        idx_ay = ctri.vindices[0]*2+1;
        idx_bx = ctri.vindices[1]*2;
        idx_by = ctri.vindices[1]*2+1;
        idx_cx = ctri.vindices[2]*2;
        idx_cy = ctri.vindices[2]*2+1;

        //CONVERT METHOD TO 3D to get alpha and beta. oops. o.O
        //maybe assume coplanar, and just do distance from barycenter?
        alpha = 1.f;
        beta  = 1.f;

        //ax ax = alpha^2
        H00.addValue(idx_ax, idx_ax, H00.getValue(idx_ax, idx_ax) + alpha*alpha);
        //ax bx = alpha*beta
        H00.addValue(idx_ax, idx_bx, H00.getValue(idx_ax, idx_bx) + alpha*beta);
        //ax cx = - alpha^2 - alpha*beta
        H00.addValue(idx_ax, idx_cx, H00.getValue(idx_ax, idx_cx) - alpha*alpha - alpha*beta);
        //ay ay = alpha^2
        H00.addValue(idx_ay, idx_ay, H00.getValue(idx_ay, idx_ay) + alpha*alpha);
        //ay by = alpha*beta
        H00.addValue(idx_ay, idx_by, H00.getValue(idx_ay, idx_by) + alpha*beta);
        //ay cy = - alpha^2 - alpha*beta
        H00.addValue(idx_ay, idx_cy, H00.getValue(idx_ay, idx_cy) - alpha*alpha - alpha*beta);
        //bx ax = alpha*beta
        H00.addValue(idx_bx, idx_ax, H00.getValue(idx_bx, idx_ax) + alpha*beta);
        //bx bx = beta^2
        H00.addValue(idx_bx, idx_bx, H00.getValue(idx_bx, idx_bx) + alpha*beta);
        //bx cx = -alpha*beta - beta^2
        H00.addValue(idx_bx, idx_cx, H00.getValue(idx_bx, idx_cx) - beta*beta - alpha*beta);
        //by ay = alpha*beta
        H00.addValue(idx_by, idx_ay, H00.getValue(idx_by, idx_ay) + alpha*beta);
        //by by = beta^2
        H00.addValue(idx_by, idx_by, H00.getValue(idx_by, idx_by) + alpha*beta);
        //by cy = -alpha*beta - beta^2
        H00.addValue(idx_by, idx_cy, H00.getValue(idx_by, idx_cy) - beta*beta - alpha*beta);
        //cx ax = - alpha^2 - alpha*beta
        H00.addValue(idx_cx, idx_ax, H00.getValue(idx_cx, idx_ax) - alpha*alpha - alpha*beta);
        //cx bx = -alpha*beta - beta^2
        H00.addValue(idx_cx, idx_bx, H00.getValue(idx_cx, idx_bx) - beta*beta - alpha*beta);
        //cx cx = (alpha+beta)^2
        H00.addValue(idx_cx, idx_cx, H00.getValue(idx_cx, idx_cx) + (alpha+beta)*(alpha+beta));
        //cy ay = - alpha^2 - alpha*beta
        H00.addValue(idx_cy, idx_ay, H00.getValue(idx_cy, idx_ay) - alpha*alpha - alpha*beta);
        //cy by = -alpha*beta - beta^2
        H00.addValue(idx_cy, idx_by, H00.getValue(idx_cy, idx_by) - beta*beta - alpha*beta);
        //cy cy = (alpha+beta)^2
        H00.addValue(idx_cy, idx_cy, H00.getValue(idx_cy, idx_cy) + (alpha+beta)*(alpha+beta));

        //filling in f0.
        f0[idx_ax] += -2*alpha; //-2alpha*s_x
        f0[idx_ay] += -2*alpha; //-2alpha*t_y
        f0[idx_bx] += -2*beta;  //-2beta*s_x
        f0[idx_by] += -2*beta;  //-2beta*t_y
        f0[idx_cx] += -2*alpha-2*beta; //-2*alpha*s_x-2*beta*s_x
        f0[idx_cy] += -2*alpha-2*beta; //-2*alpha*t_y-2*beta*t_y
    }

    //making H01, H10
    SparseMatrix H10 = SparseMatrix(2,model->numvertices*2);
    SparseMatrix H01 = SparseMatrix(model->numvertices*2,2);

    GLMtriangle seed = T(patch.at(0));

    //adding seed barycenter = 1/3(phi(A) + phi(B) + phi(C))
    H10.addValue(0, seed.vindices[0]*2,   0.166);
    H10.addValue(1, seed.vindices[0]*2+1, 0.166);
    H10.addValue(0, seed.vindices[1]*2,   0.166);
    H10.addValue(1, seed.vindices[1]*2+1, 0.166);
    H10.addValue(0, seed.vindices[2]*2,   0.166);
    H10.addValue(1, seed.vindices[2]*2+1, 0.166);

    H01.addValue(seed.vindices[0]*2,   0, 0.166);
    H01.addValue(seed.vindices[0]*2+1, 1, 0.166);
    H01.addValue(seed.vindices[1]*2,   0, 0.166);
    H01.addValue(seed.vindices[1]*2+1, 1, 0.166);
    H01.addValue(seed.vindices[2]*2,   0, 0.166);
    H01.addValue(seed.vindices[2]*2+1, 1, 0.166);

    //Hprime.
    SparseMatrix *Hprime = new SparseMatrix(H00 + H00.getTranspose());
    LinearSolver *m_orientSolver = new LinearSolver(Hprime);

    //D
    SparseMatrix *D = new SparseMatrix(H01 + H10.getTranspose());

    //calculating D*q
    double DQf[model->numvertices*2];
    double q[2];
    D->multiply(DQf, q, model->numvertices*2, 1);

    //calculate -D*q-f0
    for (int row = 0; row < D->getM(); row++) {
        DQf[row] = -DQf[row] - f0[row];
    }

    double unconstrained[model->numvertices*2];
    m_orientSolver->solve(DQf, unconstrained);

}

double2 relativeCoord(double2 A, double2 B, double2 C) {
    double2 eCA = A-C;
    double2 eBA = A-B;
    double2 eBAs = double2(-eBA.y, eBA.x);
    double d = eBA.x*eBA.x+eBA.y*eBA.y;
    double x = (eCA.x*eBA.x+eCA.y*eBA.y) / d;
    double y = (eCA.x*eBAs.x+eCA.y*eBAs.y) / d;

    return double2(x,y);
}

/* baryCoord: Barycentric Coordinates
 * To express phi(T) as a linear combination of points phi(A), phi(B), phi(C),
 * phi(T) = alpha * phi(A) + beta * phi(B) + gamma * phi(C)
 * alpha + beta + gamma = 0;
 * assuming T is coplanar with A, B, C
 */
double2 baryCoord(double2 A, double2 B, double2 C, double2 T) {
    double det   = 1/((A.x-C.x)*(B.y-C.y)-(A.y-C.y)*(B.x-C.y));
    double alpha = ((B.y-C.y)*T.x + (-B.x+C.x)*T.y)*det;
    double beta  = ((-A.y+C.y)*T.x + (A.x-C.x)*T.y)*det;
    return double2(alpha, beta);
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
        camera_.eye += (camera_.center - camera_.eye).getNormalized() * dx * .00005;
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

