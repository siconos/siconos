// Copyright (C) INRIA 1999-2008
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
// @file ActuationModel/NoDynamics/TaskFunctionControl/Trajectory.scilab
// @author Florence Billet
//
// Affiliation(s): INRIA, team BIPOP
//
// Email(s): Florence.Billet@inria.fr
// 
// @brief Compute the position, velocity and acceleration desired at a given time t
// 


// Description:
// 


//%
// Compute the position, velocity and acceleration desired at a given time t
//  @param t (float) time
//
//  @return position (float vector, size = NDOF) desired position in task space
//  @return velocity (float vector, size = NDOF) desired velocity in task space
//  @return acceleration (float vector, size = NDOF) desired acceleration
//  in task space
//  @return contacts (int vector, variable size <= NCONTACTS ) number of the
//  lines of contact vector which correspond to effectives contacts at t.
//
function [position, velocity, acceleration, contacts] = Trajectory(t)

if(CVERSION) then
  if getTrajectoryName() == 'RX90Circle' then
    taille = [6, 1];
  end;
  if getTrajectoryName() == 'PA10Infinity' then
    taille = [7, 1];
  end;
 [position_c, velocity_c, acceleration_c, contacts_c] = fort("trajectory", t, 1, "d", "out", taille, 2, "d", taille, 3, "d", taille, 4, "d", [1,1], 5, "i");
end;

if(~CVERSION | DEBUG) then
  if getTrajectoryName() == 'RX90Circle' then

    a = %pi/3;
    r = 0.2;
    position = zeros(6, 1);
    velocity = zeros(6, 1);
    acceleration = zeros(6, 1);   

    position(1) = (0.45 + 0.07)*cos(%pi/2 - %pi/3) + ...
	0.38*cos(%pi/2 -  %pi/3 + %pi/6);
    position(2) = r*cos(a*t) + 0.42 + (0.45 + 0.07)*sin(%pi/2 - %pi/3) + ...
	0.38*sin(%pi/2 -  %pi/3 + %pi/6) - r;
    position(3) = r*sin(a*t);  
    position(4) = 0.1*cos(%pi/2 - %pi/3);
    position(5) = 0.15*cos(%pi/2 - %pi/3); 
    
    velocity(1) = 0;
    velocity(2) = - r*a*sin(a*t);
    velocity(3) = r*a*cos(a*t);
    
    acceleration(1) = 0;
    acceleration(2) = - r*a^2*cos(a*t);
    acceleration(3) = - r*a^2*sin(a*t);
    
    contacts = 0; 
    
  elseif getTrajectoryName() == 'PA10Infinity' then 
    
    a = %pi/3;
    ry = 0.4;
    rz = 0.3;
    
    position = zeros(7, 1);
    velocity = zeros(7, 1);
    acceleration = zeros(7, 1);

    position(1) = 0.45*cos(17*%pi/32 - %pi/2) + 0.48*sin(%pi/6);
    position(2) =  ry*sin(a*t)*cos(a*t)/(1 + (cos(a*t))^2) + ...
	0.48*cos(%pi/6) - 0.45*sin(17*%pi/32 - %pi/2) + 0.335;
    position(3) = rz*sin(a*t)/(1 + (cos(a*t))^2);
    position(6) = 0.158*sin(%pi/2 - (%pi/4 + %pi/6));
    
    velocity(1) = 0;
    velocity(2) = ry*(cos(a*t)^2*a/(1 + cos(a*t)^2) - ...
		      sin(a*t)^2*a/(1 + cos(a*t)^2) + ...
		      2*sin(a*t)^2*cos(a*t)^2*a/(1 + cos(a*t)^2)^2);
    
    velocity(3) = rz*(cos(a*t)*a/(1 + cos(a*t)^2) + ...
		      2*sin(a*t)^2*cos(a*t)*a/(1 + cos(a*t)^2)^2);
    
    
    acceleration(1) = 0;
    acceleration(2) = ry*(-4*cos(a*t)*a^2*sin(a*t)/(1 + cos(a*t)^2) + ...
			  6*cos(a*t)^3*a^2*sin(a*t)/(1 + cos(a*t)^2)^2 - ...
			  6*sin(a*t)^3*a^2*cos(a*t)/(1 + cos(a*t)^2)^2 + ...
			  8*sin(a*t)^3*cos(a*t)^3*a^2/(1 + cos(a*t)^2)^3);
    acceleration(3) = rz*(-sin(a*t)*a^2/(1 + cos(a*t)^2) + ...
			  6*cos(a*t)^2*a^2*sin(a*t)/(1 + cos(a*t)^2)^2 + ...
			  8*sin(a*t)^3*cos(a*t)^2*a^2/(1 + cos(a*t)^2)^3 - ...
			  2*sin(a*t)^3*a^2/(1 + cos(a*t)^2)^2);
    
    contacts = 0;  
    
  end;
end;

if(CVERSION) then
 if(DEBUG) then
  if(or(abs(position_c-position)>EPS)) then
   printf("max(abs(position_c-position)) = %g\n", max(abs(position_c-position)));
  end;
  if(or(abs(velocity_c-velocity)>EPS)) then
   printf("max(abs(velocity_c-velocity)) = %g\n", max(abs(velocity_c-velocity)));
  end;
  if(or(abs(acceleration_c-acceleration)>EPS)) then
   printf("max(abs(acceleration_c-acceleration)) = %g\n", max(abs(acceleration_c-acceleration)));
  end;
  if(contacts_c ~= contacts) then
   printf("contacts_c different de contacts\n");
  end;
 end;
position = position_c;
velocity = velocity_c;
acceleration = acceleration_c;
contacts = contacts_c;
end;
  
endfunction
