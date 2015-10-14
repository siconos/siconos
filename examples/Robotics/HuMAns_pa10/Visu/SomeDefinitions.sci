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
// @file LagrangianModel/Bip/SomeDefinitions.scilab
// @author Pierre-Brice Wieber
//
// Affiliation(s): INRIA, team BIPOP
//
// Email(s): Pierre-Brice.Wieber@inria.fr
// 
// @brief Define the Lagrangian Model Name, the contacts names, the contact solids and
// the Stickfigure and the name of H-Anim Joints used in VRML Visualization
// 
 
// Description:
// 
// Modifications:
// $Log: SomeDefinitions.sci,v $
// Revision 1.5  2005/07/25 16:15:58  billet
// Adding of VRML Visualization in LagrangianModel repertory in order to have a less complex tree structure
//
// Revision 1.4  2005/05/04 13:13:15  billet
// Correction and changes for doxygen documentation
//
// Revision 1.3  2005/05/03 11:45:18  billet
// Comments for doxygen documentation
//
// Revision 1.2  2005/04/07 11:41:30  wieber
// Added the variables ContactSolids (for use in Visu) and LagrangianModelName (for later use).
//
// Revision 3.0.0.1  2005/02/08 13:05:30  rpissard
// version start HuMAnS
//
// 


//Definitions for Lagrangian Model

LagrangianModelName = 'Pa10';



StickFigure = [1, 3;...
	       3, 5];



