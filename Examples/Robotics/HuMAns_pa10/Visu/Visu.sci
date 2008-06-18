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
// $Log: Visu.sci,v $
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
// Revision 3.0.0.1  2005/02/08 13:08:22  rpissard
// version start HuMAnS
//
// 

//%
// Draw the visualisation 3D of a position or a trajectory with contact forces if they exist.
//  @param Q (float matrix, size = NDOF x NSAMPLES) trajectory or position
//  @param LAMBDA (float matrix) contact forces
//  @param POV (float vector, size = 2) graphics parameter
//  @param Message (string, optionnal) Message to print above the figure
//
function [] = Visu(Q,Q1, LAMBDA, POV, Message)
///////////////////////////////////////////
// Draws a 3D stick figure of the system following the trajectory Q. If
// LAMBDA is available, also draws the contact forces with a ratio of
// 1000N for 1m.
///////////////////////////////////////////

  if ~exists('POV', 'all') then
    POV = [13, 77];
  end;

  ///////////////////////////////////////////
  // Compute the bounding box.
  ///////////////////////////////////////////
  mmin = [%inf, %inf, %inf];
  mmax = -mmin;
  for k = 1:size(Q, 2),
    tag = Tags(Q(:, k));
    mmin = min([mmin; tag], 'r');
    mmax = max([mmax; tag], 'r');
  end;
  mmin(2) = 0;
  delta_iso = max(mmax-mmin);
  mmin = (mmin+mmax)/2-0.55*delta_iso;
  mmin(2) = 0;
  mmax = mmin+1.1*delta_iso;

  ///////////////////////////////////////////
  // Generate the view(s).
  ///////////////////////////////////////////
  xset('window', 2000);
  xname('2D Stick Figure (HuMAnS)');
  set('figure_style', 'old');
  xtape('clear', 2000);
  xset('wdim', 600, 600);
  xsetech(frect = [-1, 1, -1, 1]);
	
  xset('pixmap', 1);
  //f = gcf();
  //f.pixmap = 'on';
  for k = 1:size(Q, 2),
    tag = Tags(Q(:, k));
    tag1 = Tags(Q1(:, k));
    //clear_pixmap();
    xset("wwpc");
    ///////////////////////////////////////////
    // Draw the center of mass of the system,
    ///////////////////////////////////////////
   // param3d1([tag($,3)-0.02;tag($,3)+0.02;tag($,3)+0.02;tag($,3)-0.02;tag($,3)-0.02],...
	  //   [tag($,1)-0.02;tag($,1)+0.02;tag($,1)-0.02;tag($,1)+0.02;tag($,1)-0.02],...
	   //  list([tag($,2);tag($,2);tag($,2);tag($,2);tag($,2)],21), POV(1), POV(2),...
	    // "Z@X@Y", [1,4], [mmin(3),mmax(3),mmin(1),mmax(1),mmin(2),mmax(2)]);

    if exists('LAMBDA', 'local')|(exists('LAMBDA', 'all')&~exists('Q', 'local')) then
      contact = Contact(Q(:, k));
      contact = [contact(1:$/3), contact($/3+1:2*$/3), contact(2*$/3+1:$)];
      lambda = [LAMBDA(1:$/3, k), LAMBDA($/3+1:2*$/3, k), LAMBDA(2*$/3+1:$, k)];
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
		   POV(1), POV(2), "Z@X@Y", [0,0]);
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
		 POV(1), POV(2), "Z@X@Y", [0,0]);
      end;
    end;
    ///////////////////////////////////////////
    // and finally draw the 3D stick figure of the system.
    ///////////////////////////////////////////
    //param3d1([tag(StickFigure(:,1),1)';tag(StickFigure(:,2),1)'],...
	   //  [tag(StickFigure(:,1),2)';tag(StickFigure(:,2),2)'],...
	    // [tag(StickFigure(:,1),3)';tag(StickFigure(:,2),3)'],...
	     //POV(1), POV(2), "X@Z@Y", [0,0]);



		x = [tag(1,1), tag(3,1), tag(5,1)];
		z = [tag(1,3), tag(3,3), tag(5,3)];
		
		plot2d(x, z, rect = [-1, -1, 1, 1]);

		x1 = [tag1(1,1), tag1(3,1), tag1(5,1)];
		z1= [tag1(1,3), tag1(3,3), tag1(5,3)];
		
		plot2d(x1, z1, rect = [-1, -1, 1, 1], style=3);




    if exists('Message') then
      [a, b] = xgetech();
      xstring(0, b($), Message);
    end;
    //show_pixmap();
    xset("wshow");
  end;
  xset('pixmap', 0);

endfunction
