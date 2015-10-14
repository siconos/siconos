// Copyright (C) INRIA 1999-2005
// 
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 2 as published
// by the Free Software Foundation.
// 
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
// Public License for more details.
// 
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
//%
// @file LagrangianDynamics/Complete/Visu.scilab
// @author Pierre-Brice Wieber
//
// Affiliation(s): INRIA, team BIPOP
//
// Email(s): Pierre-Brice.Wieber@inria.fr
// 
// @brief Draw the visualisation 3D
// 

// Description:
// 
// Modifications:
// $Log$
// Revision 1.6  2005/10/28 08:14:14  billet
// Moving from old to new graphic style and removing of windows scilab visualization bug
//
// Revision 1.5  2005/05/04 13:13:15  billet
// Correction and changes for doxygen documentation
//
// Revision 1.4  2005/05/03 11:45:17  billet
// Comments for doxygen documentation
//
// Revision 1.3  2005/05/02 12:03:13  wieber
// Added the possibility to print a Message above the figure and prepared the new figure_style mode, still not active because buggy in Scilab 3.0.
//
// Revision 1.2  2005/04/07 15:37:29  wieber
// Drawing contact forces on different solids through the use of the variable ContactSolids.
//
// Revision 1.1.1.1  2005/02/08 13:08:22  rpissard
// version start HuMAnS
//
// 


function [] = VisuPlay(flag, gwin)
///////////////////////////////////////////
// Draws a 3D stick figure of the system following the trajectory Q1. 
// If LAMBDA is available, also draws the contact forces with a 
// ratio of 1000N for 1m.
///////////////////////////////////////////
  
  global QVISU;
  global LAMBDAVISU;
  
  ///////////////////////////////////////////
  // Compute the bounding box.
  ///////////////////////////////////////////
  mmin = [%inf, %inf, %inf];
  mmax = -mmin;
  for k = 1:size(QVISU, 2),
    tag = Tags(QVISU(:, k));

    mmin = min([mmin; tag], 'r');
    mmax = max([mmax; tag], 'r');
  end;
  mmin(2) = 0;
  delta_iso = max(mmax-mmin);
  mmin = (mmin+mmax)/2-0.55*delta_iso;
  mmin(2) = 0;
  mmax = mmin+1.1*delta_iso;

  f = scf(gwin);	
  f.figure_name = "3D Stick Figure (HuMAnS)";
  ha = f.children;
  ha.axes_visible = "on";
  ha.box = "on";
  ha.view = "3d";
  f.pixmap = "on";
  ha.data_bounds = [mmin(3), mmin(1), mmin(2); mmax(3), mmax(1), mmax(2)];	
  
  if flag~=0  then
    delete(ha.children);
  else
    f.figure_size=[600,600];
    ha.rotation_angles = [80, 30];
  end;

  for k = 1:size(QVISU, 2),
    tag = Tags(QVISU(:, k));
    if k>1 then
      delete(ha.children);
    end;
    
    ///////////////////////////////////////////
    // Draw the center of mass of the system,
    ///////////////////////////////////////////
    
    param3d1([tag($,3)-0.02;tag($,3)+0.02;tag($,3)+0.02;tag($,3)-0.02;tag($,3)-0.02],...
	     [tag($,1)-0.02;tag($,1)+0.02;tag($,1)-0.02;tag($,1)+0.02;tag($,1)-0.02],...
	     list([tag($,2);tag($,2);tag($,2);tag($,2);tag($,2)],21), ...
	     ha.rotation_angles(2), ha.rotation_angles(1),"Z@X@Y", [1,4]);
    
    if exists('LAMBDAVISU', 'all')&(size(LAMBDAVISU, 2)==size(QVISU, 2)) then
      contact = Contact(QVISU(:, k));
      contact = [contact(1:$/3), contact($/3+1:2*$/3), contact(2*$/3+1:$)];
      lambda = [LAMBDAVISU(1:$/3, k), LAMBDAVISU($/3+1:2*$/3, k), LAMBDAVISU(2*$/3+1:$, k)];
      
      /////////////////////////////////////////////
      // draw the sum of the contact forces acting on each contacting solid
      // at the corresponding center of pressure,
      /////////////////////////////////////////////
      for i  = 1:size(ContactSolids, 1),
	ContactPoints = ContactSolids(i, 1):ContactSolids(i, 2);
	Force = lambda(ContactPoints, :);
	if sum(Force(:, 2))>sqrt(%eps) then
	  COP = sum(diag(Force(:, 2))*contact(ContactPoints, :), 'r')/sum(Force(:, 2));
	  Force = COP+sum(Force/1000, 'r');
	  
	  param3d1([COP(3)-0.02;COP(3)+0.02;COP(3)+0.02;COP(3)-0.02;COP(3)-0.02;COP(3);Force(3)],...
		   [COP(1)-0.02;COP(1)+0.02;COP(1)-0.02;COP(1)+0.02;COP(1)-0.02;COP(1);Force(1)],...
		   list([COP(2);COP(2);COP(2);COP(2);COP(2);COP(2);Force(2)],24),...
		   ha.rotation_angles(2), ha.rotation_angles(1), "Z@X@Y", [0,0]);
	end;
      end;
      
      ///////////////////////////////////////////
      // draw the sum of all the contact forces acting on the system at the
      // global center of pressure on the ground (aka the zero moment point),
      ///////////////////////////////////////////
      if sum(lambda(:, 2))>sqrt(%eps) then
	contact = contact(lambda(:, 2)>sqrt(%eps), :);
	lambda = lambda(lambda(:, 2)>sqrt(%eps), :);
	contact = contact-diag(contact(:, 2)./lambda(:, 2))*lambda;
	COP = sum(diag(lambda(:, 2))*contact, 'r')/sum(lambda(:, 2));
	lambda = COP+sum(lambda/1000, 'r');
	
	param3d1([COP(3)-0.02;COP(3)+0.02;COP(3)+0.02;COP(3)-0.02;COP(3)-0.02;COP(3);lambda(3)],...
		 [COP(1)-0.02;COP(1)+0.02;COP(1)-0.02;COP(1)+0.02;COP(1)-0.02;COP(1);lambda(1)],...
		 list([COP(2);COP(2);COP(2);COP(2);COP(2);COP(2);lambda(2)],11),...
		 ha.rotation_angles(2), ha.rotation_angles(1), "Z@X@Y", [0,0]);
      end;
    end;
    

    ///////////////////////////////////////////
    // and finally draw the 3D stick figure of the system.
    ///////////////////////////////////////////
    
    param3d1([tag(StickFigure(:,1),3)'; tag(StickFigure(:,2),3)'], ...
	     [tag(StickFigure(:,1),1)'; tag(StickFigure(:,2),1)'], ...
	     [tag(StickFigure(:,1),2)'; tag(StickFigure(:,2),2)'], ...
	     ha.rotation_angles(2), ha.rotation_angles(1), "Z@X@Y", [0,0]);

    if exists('Message') then
      xtitle(Message);
    end;
    show_pixmap();

  end;

endfunction


//%
// Draw the visualisation 3D of a position or a trajectory with contact
// forces if they exist.
//  @param Q (float matrix, size = NDOF x NSAMPLES) trajectory or position
//  @param LAMBDA (float matrix) contact forces
//  @param Message (string, optionnal) Message to print above the figure
//  @param DuringSimulationFlag (boolean, optionnal) Flag to know if we draw 
//  the figure during the simulation or not. 
//
function [] = Visu(Q, LAMBDA, Message, DuringSimulationFlag)

  global QVISU;
  global LAMBDAVISU;

  QVISU = Q;
  if exists('LAMBDA', 'local')|(~exists('LAMBDA', 'local')&exists('LAMBDA', 'all')) then
    LAMBDAVISU = LAMBDA;
  end;

  if ~exists('DuringSimulationFlag', 'local') then
    DuringSimulationFlag = %F;
  end;
  if DuringSimulationFlag then
    scf(2005);
    clf();
  end

  VisuPlay(0, 2005);

  if (~DuringSimulationFlag) then 
    addmenu(2005, "Replay", list(2, 'VisuPlay'));
  end;

endfunction
