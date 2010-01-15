/* Siconos-sample version 3.1.0, Copyright INRIA 2005-2009.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 *
 */

#include "SpheresViewer.hpp"

using namespace qglviewer;

void SpheresViewer::init()
{
  // viewer and scene state
  restoreStateFromFile();
  setSceneRadius(100.);
  showEntireScene();
  setGridIsDrawn();

  // color
  setBackgroundColor(QColor(200, 200, 200));
  setForegroundColor(QColor(0, 0, 0));

  // camera
  camera()->setPosition(Vec(1., -1.0, 0.5));
  //camera()->lookAt(sceneCenter());
  //camera()->setType(Camera::ORTHOGRAPHIC);
  camera()->showEntireScene();
  setManipulatedFrame(camera()->frame());

  // help screen
  help();

  // siconos setup
  Siconos_.reset(new Spheres());
  Siconos_->init();

  int i;
  DSIterator itDS;

  DynamicalSystemsGraph::VIterator dsi, dsend;
  boost::tie(dsi, dsend) = GETALLDS(Siconos_)->vertices();

  for (i = 0; dsi != dsend; ++i, ++dsi)
  {

    boost::shared_ptr<SetDrawing> setdrawing(new SetDrawing(*this, i));
    GETALLDS(Siconos_)->bundle(*dsi)->accept(setdrawing);
  }

  bodydraw_.reset(new BodyDraw(*this));

  sbodydraw_.reset(new SelectedBodyDraw(*this));

  lastSelected_ = -1;

  myMouseBehavior_ = false;

  startAnimation();
}


void SpheresViewer::draw()
{

  int i, mrow;

  char qs[6];

  DSIterator itDS;
  SP::DynamicalSystemsSet involvedDS;
  SP::Interaction interaction;
  SP::Relation relation;


  for (i = 0; i < GETNDS(Siconos_); i++)
  {
    if (drawings_[i]->selected())
    {
      drawings_[i]->getDS()->accept(sbodydraw_);
    }
    else
    {
      drawings_[i]->getDS()->accept(bodydraw_);
    }
  }

  glColor3f(.45, .45, .45);
  glLineWidth(1.);
  drawGrid(100, 200);
  setGridIsDrawn();

  glColor3f(.1, .1, .3);
  drawVec(-100, 0, 100, 0);
  drawVec(0, -100, 0, 100);

  glColor3f(0, 0, 1);
  glLineWidth(4.);

  for (unsigned int i = 0 ; i < Siconos_->getPlansPtr()->size(0) ; ++i)
  {
    double A = (*Siconos_->getPlansPtr())(i, 0);
    double B = (*Siconos_->getPlansPtr())(i, 1);
    double C = (*Siconos_->getPlansPtr())(i, 2);
    double D = (*Siconos_->getPlansPtr())(i, 3);

    double N2 = A * A + B * B + C * C;

    glPushMatrix();
    glLineWidth(1.);
    glTranslatef(-D * A / N2, -D * B / N2, -D * C / N2);
    glMultMatrixd(Quaternion(Vec(0., 0., 1.), Vec(A, B, C)).matrix());
    drawGrid(100, 200);
    glPopMatrix();

  }

  glColor3f(.1, .1, .1);
  glLineWidth(4.);
  QGLViewer::drawArrow(qglviewer::Vec(0, 0, .1), qglviewer::Vec(1, 0, .1), .01, 3);
  QGLViewer::drawArrow(qglviewer::Vec(0, 0, .1), qglviewer::Vec(0, 1, .1), .01, 3);

  glLineWidth(1.);
  for (i = -100; i <= 100; i += 5)
  {
    sprintf(qs, "%d", i);
    //    print((float)i,-.8,qs,small_text);
    //print(-.8,(float)i,qs,small_text);
    drawVec((float)i, -.2, (float)i, .2);
    drawVec(-.2, (float)i, .2, (float)i);
  }
  for (i = -100; i <= 100; i++)
  {
    drawVec((float)i, -.1, (float)i, .1);
    drawVec(-.1, (float)i, .1, (float)i);
  }
}


void SpheresViewer::animate()
{

  // THREAD or PROCESS ?
  //clock_gettime(CLOCK_THREAD_CPUTIME_ID,&ts1);

  Siconos_->compute();
  Siconos_->compute();
  //saveSnapshot();

  //clock_gettime(CLOCK_THREAD_CPUTIME_ID,&ts2);

  //timeSiconos = ts2.tv_nsec - ts1.tv_nsec;
  //std::cout << timeSiconos << std::endl;
}


void SpheresViewer::mousePressEvent(QMouseEvent* e)
{
  if ((e->button() == Qt::LeftButton) && (e->modifiers() == Qt::NoButton))
    myMouseBehavior_ = true;
  else
    QGLViewer::mousePressEvent(e);
}

void SpheresViewer::mouseMoveEvent(QMouseEvent *e)
{

  float m_coor[3];
  float c_coor[3] = { 0. , 0., 0.};
  float coor[3] = { 0. , 0., 0.};

  if (selectedName() >= 0)
  {

    SP::Drawing selected_drawing = drawings_[selectedName()];
    SP::LagrangianDS DS = selected_drawing->getDS();
    double mass = DS->mass()->getValue(0, 0);

    QPoint pixel;
    qglviewer::Vec orig;
    qglviewer::Vec dir;

    lastSelected_ = selectedName();
    m_coor[0] = (float) e->x();
    m_coor[1] = (float) e->y();
    m_coor[2] = 0;

    pixel.setX(e->x());
    pixel.setY(camera()->screenHeight() - e->y());

    //     camera()->convertClickToLine(pixel,orig,dir);

    camera()->getUnprojectedCoordinatesOf(m_coor, m_coor, NULL);

    coor[0] = m_coor[0];
    coor[1] = m_coor[1];
    coor[2] = m_coor[2];

    // toward camera
    DS->fExt()->setValue(0, (coor[0] - DS->q()->getValue(0))*mass * 10);
    DS->fExt()->setValue(1, (coor[1] - DS->q()->getValue(1))*mass * 10);
    DS->fExt()->setValue(2, (coor[2] - DS->q()->getValue(2))*mass * 10);
  }

  QGLViewer::mouseMoveEvent(e);
}

void SpheresViewer::mouseReleaseEvent(QMouseEvent* e)
{
  if (myMouseBehavior_)
    myMouseBehavior_ = false;

  if (lastSelected_ >= 0)
  {
    drawings_[lastSelected_]->restoreFExt();
    drawings_[lastSelected_]->nextSelection();
    lastSelected_ = -1;
    setSelectedName(-1);
    setManipulatedFrame(camera()->frame());
  };

  QGLViewer::mouseReleaseEvent(e);
}


QString SpheresViewer::helpString() const
{
  QString text("<h2>S i m p l e V i e w e r</h2>");
  text += "Use the mouse to move the camera around the object. ";
  text += "You can respectively revolve around, zoom and translate with the three mouse buttons. ";
  text += "Left and middle buttons pressed together rotate around the camera view direction axis<br><br>";
  text += "Pressing <b>Alt</b> and one of the function keys (<b>F1</b>..<b>F12</b>) defines a camera keyFrame. ";
  text += "Simply press the function key again to restore it. Several keyFrames define a ";
  text += "camera path. Paths are saved when you quit the application and restored at next start.<br><br>";
  text += "Press <b>F</b> to display the frame rate, <b>A</b> for the world axis, ";
  text += "<b>Alt+Return</b> for full screen mode and <b>Control+S</b> to save a snapshot. ";
  text += "See the <b>Keyboard</b> tab in this window for a complete shortcut list.<br><br>";
  text += "Double clicks automates single click actions: A left button double click aligns the closer axis with the camera (if close enough). ";
  text += "A middle button double click fits the zoom of the camera and the right button re-centers the scene.<br><br>";
  text += "A left button double click while holding right button pressed defines the camera <i>Revolve Around Point</i>. ";
  text += "See the <b>Mouse</b> tab and the documentation web pages for details.<br><br>";
  text += "Press <b>Escape</b> to exit the viewer.";
  return text;
}

