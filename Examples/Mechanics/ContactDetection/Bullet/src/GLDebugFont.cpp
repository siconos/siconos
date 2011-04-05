/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2009 Erwin Coumans  http://bulletphysics.org

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it freely,
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

#include "GLDebugFont.h"


#ifdef _WIN32//for glut.h
#include <windows.h>
#endif

//think different
#if defined(__APPLE__) && !defined (VMDMESA)
#include <TargetConditionals.h>
#if (defined (TARGET_OS_IPHONE) && TARGET_OS_IPHONE) || (defined (TARGET_IPHONE_SIMULATOR) && TARGET_IPHONE_SIMULATOR)
#import <OpenGLES/ES1/gl.h>
#define glOrtho glOrthof
#else
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#endif
#else



#ifdef _WINDOWS
#include <windows.h>
#include <GL/gl.h>
#include <GL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glut.h>
#endif
#endif

#include <stdio.h>
#include <string.h> //for memset

extern unsigned char sFontData[];
static bool sTexturesInitialized = false;

static GLuint sTexture = -1;
static int sScreenWidth = -1;
static int sScreenHeight = -1;


void GLDebugResetFont(int screenWidth, int screenHeight)
{

  if ((sScreenWidth == screenWidth) && (sScreenHeight == screenHeight))
    return;

  sScreenWidth = screenWidth;
  sScreenHeight = screenHeight;

  if (!sTexturesInitialized)
  {
    sTexturesInitialized = true;
    glGenTextures(1, &sTexture);
    glBindTexture(GL_TEXTURE_2D, sTexture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, 3, 256 , 256 , 0, GL_RGB, GL_UNSIGNED_BYTE, &sFontData[0]);
  }

  printf("generating font at resolution %d,%d\n", screenWidth, screenHeight);

}

#define USE_ARRAYS 1

void  GLDebugDrawStringInternal(int x, int y, const char* string, const btVector3& rgb)
{
  GLDebugDrawStringInternal(x, y, string, rgb, true, 10);
}

void  GLDebugDrawStringInternal(int x, int y, const char* string, const btVector3& rgb, bool enableBlend, int spacing)
{

  if (!sTexturesInitialized)
  {
    GLDebugResetFont(sScreenWidth, sScreenHeight);
  }
  if (strlen(string))
  {

    glColor4f(rgb.getX(), rgb.getY(), rgb.getZ(), 1.f);
    float cx;
    float cy;

    glMatrixMode(GL_TEXTURE);
    glLoadIdentity();

    glDisable(GL_TEXTURE_GEN_S);
    glDisable(GL_TEXTURE_GEN_T);
    glDisable(GL_TEXTURE_GEN_R);

    glEnable(GL_TEXTURE_2D);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);
    glDepthFunc(GL_LEQUAL);

    if (enableBlend)
    {
      glEnable(GL_BLEND);
    }
    else
    {
      glDisable(GL_BLEND);
    }
    glEnable(GL_DEPTH_TEST);
    glBindTexture(GL_TEXTURE_2D, sTexture);
    glDisable(GL_DEPTH_TEST);
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();

    glOrtho(0, sScreenWidth, 0, sScreenHeight, -1, 1);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glTranslatef(btScalar(x), btScalar(sScreenHeight - y), btScalar(0));

#if USE_ARRAYS

    glDisableClientState(GL_COLOR_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_TEXTURE_COORD_ARRAY);
#endif

    GLfloat verts[] =
    {
      0.0f, 1.0f, 0.0f,
      -1.0f, -1.0f, 0.0f,
      1.0f, -1.0f, 0.0f,
      0.f, 0.f, 0.f
    };

    GLfloat uv_texcoords[] =
    {
      0, 0,
      0, 0,
      0, 0,
      0, 0
    };
    verts[0] = 0;
    verts[1] = 0;
    verts[2] = 0;
    verts[3] = 16 - 1;
    verts[4] = 0;
    verts[5] = 0;
    verts[6] = 16 - 1;
    verts[7] = 16 - 1;
    verts[8] = 0;
    verts[9] = 0;
    verts[10] = 16 - 1;
    verts[11] = 0;

    for (int i = 0; i < int (strlen(string)); i++)
    {
      char ch = string[i] - 32;
      if (ch >= 0)
      {
        cx = float(ch % 16) * btScalar(1. / 16.f);
        cy = float(ch / 16) * btScalar(1. / 16.f);

        uv_texcoords[0] = cx;
        uv_texcoords[1] = btScalar(1 - cy - 1. / 16.f);
        uv_texcoords[2] = btScalar(cx + 1. / 16.f);
        uv_texcoords[3] = btScalar(1 - cy - 1. / 16.f);
        uv_texcoords[4] = btScalar(cx + 1. / 16.f);
        uv_texcoords[5] = btScalar(1 - cy);
        uv_texcoords[6] = cx;
        uv_texcoords[7] = btScalar(1 - cy);
#if USE_ARRAYS
        glTexCoordPointer(2, GL_FLOAT, 0, uv_texcoords);
        glVertexPointer(3, GL_FLOAT, 0, verts);
        glDrawArrays(GL_TRIANGLE_FAN, 0, 4);
#else
        glBegin(GL_QUADS);
        glTexCoord2f(cx, 1 - cy - 1. / 16.f);

        glVertex2i(0, 0);
        glTexCoord2f(cx + 1. / 16.f, 1 - cy - 1. / 16.f);

        glVertex2i(16 - 1, 0);
        glTexCoord2f(cx + 1. / 16.f, 1 - cy);

        glVertex2i(16 - 1, 16 - 1);
        glTexCoord2f(cx, 1 - cy);

        glVertex2i(0, 16 - 1);
        glEnd();
#endif

        glTranslatef(spacing, 0, 0);
      }
    }

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
#if 1
    glEnable(GL_DEPTH_TEST);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);
    glDepthFunc(GL_LEQUAL);
    glDisable(GL_BLEND);
    glDisable(GL_TEXTURE_2D);

    glMatrixMode(GL_TEXTURE);
    glLoadIdentity();
    glScalef(btScalar(0.025), btScalar(0.025), btScalar(0.025));
#endif
    glMatrixMode(GL_MODELVIEW);
#if USE_ARRAYS
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_TEXTURE_COORD_ARRAY);
#endif
    //glDisable(GL_TEXTURE_2D);
  }
}

void  GLDebugDrawString(int x, int y, const char* string)
{

  btVector3 rgb(1, 1, 1);
  GLDebugDrawStringInternal(x, y, string, rgb);
}


unsigned char sFontData[] =
{
  // emty see bullet src if needed
};
