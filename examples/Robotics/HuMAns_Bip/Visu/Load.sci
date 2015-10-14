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
// 
// Author(s): Pierre-Brice Wieber
// Affiliation(s): INRIA, team BIPOP
// Email(s): Pierre-Brice.Wieber@inria.fr
// 
// Description:
// 
// Modifications:
// $Log$
// Revision 1.6  2005/07/28 01:44:08  wieber
// Fixed some tabulations, aka, Florence should stop using nedit.
//
// Revision 1.5  2005/07/26 13:08:22  billet
// Modification of VRML Visualization files structure for more simplicity
//
// Revision 1.4  2005/07/25 16:15:58  billet
// Adding of VRML Visualization in LagrangianModel repertory in order to have a less complex tree structure
//
// Revision 1.3  2005/03/23 21:11:17  rpissard
// prefix libraries with lib and SCIDIR use for Makefiles
//
// Revision 1.2  2005/03/12 15:31:41  rpissard
// Unix Makefile tuning
//
// Revision 1.1.1.1  2005/02/08 13:05:34  rpissard
// version start HuMAnS
//
// 

exec('KickStart.sci');

[LIBPATH, LIBEXT] = LibTools();

exec('SomeDefinitions.sci');

idlib=link(LIBPATH+'/libLagrangianModel'+LIBEXT);

addinter(idlib, 'LagrangianGateway',...
	 ["Tags",...
	  "CreateVRML"]);


exec('Visu.sci');

