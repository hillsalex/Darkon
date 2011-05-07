#-------------------------------------------------
#
# Project created by QtCreator 2011-04-21T21:12:39
#
#-------------------------------------------------

QT       += core gui opengl

TARGET = Hatching
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp\
        src/LinearSolver.cpp\
src/SparseMatrix.cpp\
src/SparseMatrix.inl\
        src/drawengine.cpp\
        src/glwidget.cpp\
        src/shot/shot.cpp\
        src/shot/testShot.cpp\
        src/glm/glm.cpp\
        src/glm/targa.cpp\
        src/math/CS123Matrix.cpp\
        src/math/CS123Matrix.inl\
        src/math/CS123Vector.inl\
    src/shot/MPSandbox.cpp \
    src/shot/tamgenerator.cpp \
    Sphere.cpp \
    src/shot/SphereShot.cpp \
    src/meshoperator.cpp \
    LappedUtils.cpp \
    SurfacePatch.cpp \
    src/alglib/solvers.cpp \
    src/alglib/linalg.cpp \
    src/alglib/ap.cpp \
    src/alglib/alglibmisc.cpp \
    src/alglib/alglibinternal.cpp




HEADERS  += mainwindow.h\
        src/drawengine.h\
        src/LinearSolver.h\
src/SparseMatrix.h\
src/matrix.h\
        src/glwidget.h\
        src/shot/shot.h\
        src/shot/testShot.h\
        src/glm/glm.h\
        src/glm/targa.h\
        src/math/testMatrix.h \
        src/math/CS123Vector.h \
        src/math/CS123Matrix.h \
        src/math/CS123Algebra.h \
        src/lib/CS123Common.h\
        glext.h \
    src/shot/MPSandbox.h \
    src/shot/tamgenerator.h \
    Sphere.h \
    src/shot/SphereShot.h \
    src/meshoperator.h \
    LappedUtils.h \
    SurfacePatch.h \
    src/alglib/solvers.h \
    src/alglib/linalg.h \
    src/alglib/ap.h \
    src/alglib/alglibmisc.h \
    src/alglib/alglibinternal.h \
    src/alglib/stdafx.h


FORMS += mainwindow.ui
INCLUDEPATH += src src/shot src/lib src/glm \
            /course/cs224/lib/umfpack/include
DEPENDPATH += src src/shot src/lib src/glm
LIBS += -L/course/cs032/lib/OpenCV/release/lib \
        -lcv -lcxcore -lcvaux -lhighgui \
        -L/course/cs224/lib/umfpack \
        -lumfpack -lamd -lblas -lcerbla \
