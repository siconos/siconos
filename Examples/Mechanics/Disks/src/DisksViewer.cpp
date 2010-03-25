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

#include "DisksViewer.hpp"

using namespace qglviewer;

void DisksViewer::init()
{
  // viewer and scene state
  restoreStateFromFile();
  setSceneRadius(100.);
  showEntireScene();

  // color
  setBackgroundColor(QColor(200, 200, 200));
  setForegroundColor(QColor(0, 0, 0));

  // 2D setup
  glDisable(GL_LIGHTING);
  camera()->setPosition(Vec(0.0, 0.0, 1.0));
  camera()->lookAt(sceneCenter());
  camera()->setType(Camera::ORTHOGRAPHIC);
  camera()->showEntireScene();
  constraint_.reset(new WorldConstraint());
  constraint_->setRotationConstraintType(AxisPlaneConstraint::FORBIDDEN);
  camera()->frame()->setConstraint(&*constraint_);

  // help screen
  help();

  // siconos setup
  Siconos_.reset(new Disks());

  Siconos_->init();

  int i;
  DSIterator itDS;

  DynamicalSystemsGraph::VIterator dsi, dsend;
  boost::tie(dsi, dsend) = GETALLDS(Siconos_)->vertices();
  for (i = 0; dsi != dsend; ++i, ++dsi)
  {

    boost::shared_ptr<SetDrawing> setdrawing(new SetDrawing(*this, i));
    GETALLDS(Siconos_)->bundle(*dsi)->acceptSP(setdrawing);
  }


  bodydraw_.reset(new BodyDraw(*this));

  sbodydraw_.reset(new SelectedBodyDraw(*this));

  lastSelected_ = -1;

  initUCircle();

  myMouseBehavior_ = false;

  startAnimation();
}




void DisksViewer::draw()
{

  int i, mrow;
  clock_t lcptime;

  char qs[6];

  float lbd, lbdmax, w, x1, y1, x2, y2, d, r;


  DSIterator itDS;
  SP::DynamicalSystemsSet involvedDS;
  SP::UnitaryRelationsGraph I1 = Siconos_->model()->simulation()->indexSet(1);
  SP::Interaction interaction;
  SP::Relation relation;


  // calibration
  UnitaryRelationsGraph::VIterator ui, uiend;
  for (boost::tie(ui, uiend) = I1->vertices(); ui != uiend; ++ui)
  {
    lbdmax = fmax(I1->bundle(*ui)->interaction()->lambdaOld(1)->getValue(0), lbdmax);
  }

  for (boost::tie(ui, uiend) = I1->vertices(); ui != uiend; ++ui)
  {
    interaction = I1->bundle(*ui)->interaction();
    relation = interaction->relation();

    lbd = interaction->lambdaOld(1)->getValue(0);

    // screen width of interaction
    w = lbd / (2 * fmax(lbdmax, 1.)) + .03;

    involvedDS = interaction->dynamicalSystems();

    // disk/disk
    itDS = involvedDS->begin();

    SP::DynamicalSystem d1 = *itDS;
    SP::DynamicalSystem d2;
    if (involvedDS->size() == 2)
      d2 = *++itDS;
    else
      d2 = d1;

    boost::shared_ptr<LambdaFirst> lambda1(new LambdaFirst());

    boost::shared_ptr<LambdaSecond> lambda2(new LambdaSecond(*this, w));

    relation->acceptSP(lambda1);
    d1->acceptSP(lambda1);
    lambda1->acceptSP(lambda2);
    d2->acceptSP(lambda2);


    /*else
      }
      if (involvedDS->size() == 1) {

         lbd=interaction->getLambdaOldPtr(1)->getValue(0);

         itDS=involvedDS->begin();

         x1 = AS(LagrangianDS )(*itDS)->q()->getValue(0);
         y1 = AS(LagrangianDS )(*itDS)->q()->getValue(1);
         r =  AS(Disk )(*itDS)->getRadius();

         if ( relation->getSubType() == RELATION::ScleronomousR ) {

           // disk/plan
           x2 = AS(LagrangianScleronomousR)(relation)->getJacHPtr(0)->getValue(0,0);
           y2 = AS(LagrangianScleronomousR)(relation)->getJacHPtr(0)->getValue(0,1);
         } else {

           // unknown relation, do nothing
           x2 = x1;
           y2 = y1;
         }

         glColor3f(0.,0.,0.);
         d = sqrt(x2*x2+y2*y2);
         if (AS(DiskQGL)(*itDS)->drawing()->shape() == SHAPE::DISK)
         drawRec(x1,y1,x1-r*x2/d,y1-r*y2/d,w);
         else if ((AS(DiskQGL)(*itDS)->drawing()->shape() == SHAPE::CIRCLE))
           drawRec(x1-r*x2/d,y1-r*y2/d,x1-r*x2/d,y1-r*y2/d,w);
           }*/
  }

  for (i = 0; i < GETNDS(Siconos_); i++)
  {
    if (drawings_[i]->selected())
    {
      drawings_[i]->getDS()->acceptSP(sbodydraw_);
    }
    else
    {
      drawings_[i]->getDS()->acceptSP(bodydraw_);
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

  if (Siconos_->plans())
  {
    for (unsigned int i = 0 ; i < Siconos_->plans()->size(0) ; ++i)
    {
      double A = (*Siconos_->plans())(i, 0);
      double B = (*Siconos_->plans())(i, 1);
      double C = (*Siconos_->plans())(i, 2);
      double xc = (*Siconos_->plans())(i, 3);
      double yc = (*Siconos_->plans())(i, 4);
      double w = fmin(1e10, (*Siconos_->plans())(i, 5));
      double H = hypot(A, B);

      if (w == 0) w = 1e10;

      //      assert ( fabs(A*xc + B*yc + C) <= std::numeric_limits<double>::epsilon() );

      drawVec(xc, yc, xc - 0.5 * w * B / H, yc + 0.5 * w * A / H);
      drawVec(xc, yc, xc + 0.5 * w * B / H, yc - 0.5 * w * A / H);
    }
  }


  if (Siconos_->movingPlans())
  {
    double time = Siconos_->model()->currentTime();
    for (unsigned int i = 0 ; i < Siconos_->movingPlans()->size1() ; ++i)
    {
      double A = (*Siconos_->movingPlans())(i, 0)(time);
      double B = (*Siconos_->movingPlans())(i, 1)(time);
      double C = (*Siconos_->movingPlans())(i, 2)(time);
      double w = 1e10;
      double H = hypot(A, B);
      double xc, yc;

      if (fabs(C) > std::numeric_limits<double>::epsilon())
      {
        if (A == 0)
          // By+C=0
        {
          yc = -C / B;
        }
        else if (B == 0)
          // Ax+C=0
        {
          xc = -C / A;
        }
        else
          // Ax+By+C=0
        {
          if (xc != 0)
            yc = - (A * xc + C) / B;
          else
            xc = - (B * yc + C) / A;
        }
      }

      drawVec(xc, yc, xc - 0.5 * w * B / H, yc + 0.5 * w * A / H);
      drawVec(xc, yc, xc + 0.5 * w * B / H, yc - 0.5 * w * A / H);
    }
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


void DisksViewer::animate()
{

  struct timespec ts1, ts2;

  Siconos_->compute();
  Siconos_->compute();
  //saveSnapshot();

}


void DisksViewer::mousePressEvent(QMouseEvent* e)
{
  if ((e->button() == Qt::LeftButton) && (e->modifiers() == Qt::NoButton))
    myMouseBehavior_ = true;
  else
    QGLViewer::mousePressEvent(e);
}

void DisksViewer::mouseMoveEvent(QMouseEvent *e)
{

  float coor[3];

  if (selectedName() >= 0)
  {
    lastSelected_ = selectedName();
    coor[0] = (float) e->x();
    coor[1] = (float) e->y();
    coor[2] = 0;

    camera()->getUnprojectedCoordinatesOf(coor, coor, NULL);

    drawings_[selectedName()]->getDS()->fExt()
    ->setValue(0, ((float)coor[0] - drawings_[selectedName()]->getDS()->q()->getValue(0))*drawings_[selectedName()]->getDS()->mass()->getValue(0, 0) * 30);
    drawings_[selectedName()]->getDS()->fExt()
    ->setValue(1, ((float)coor[1] - drawings_[selectedName()]->getDS()->q()->getValue(1))*drawings_[selectedName()]->getDS()->mass()->getValue(0, 0) * 30);
  }

  QGLViewer::mouseMoveEvent(e);
}

void DisksViewer::mouseReleaseEvent(QMouseEvent* e)
{
  if (myMouseBehavior_)
    myMouseBehavior_ = false;

  if (lastSelected_ >= 0)
  {
    drawings_[lastSelected_]->restoreFExt();
    drawings_[lastSelected_]->nextSelection();
  };

  QGLViewer::mouseReleaseEvent(e);
}


QString DisksViewer::helpString() const
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


