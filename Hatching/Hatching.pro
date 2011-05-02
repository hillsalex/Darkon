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
    src/shot/SphereShot.cpp




HEADERS  += mainwindow.h\
        src/drawengine.h\
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
    src/shot/SphereShot.h

FORMS += mainwindow.ui
INCLUDEPATH += src src/shot src/lib src/glm
DEPENDPATH += src src/shot src/lib src/glm
