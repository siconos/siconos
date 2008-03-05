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
// $Log: KickStart.sci,v $
// Revision 1.6  2005/07/01 15:38:32  billet
// New version 1.0.3 in the banner.
//
// Revision 1.5  2005/05/19 08:59:47  billet
// Change of Version number
//
// Revision 1.4  2005/04/19 12:15:03  wieber
// New version 1.0.1 in the banner.
//
// Revision 1.3  2005/03/17 17:10:48  rpissard
// Doxygen building
//
// Revision 1.2  2005/03/12 15:31:28  rpissard
// Unix Makefile tuning
//
// Revision 3.0.0.1  2005/02/08 13:04:26  rpissard
// version start HuMAnS
//
// 


printf('        -------------------------------------------\n');
printf('                        HuMAnS-1.0.3\n');
printf('           A Scilab toolbox for Humanoid Motion\n');
printf('                  Analysis and Simulation\n');
printf('               Copyright (C) INRIA 1999-2005\n');
printf('        -------------------------------------------\n');

lines(0);
stacksize(1e7);

function todo = LoadModule(path)
  todo = "oldpath = getcwd();"+...
	 "chdir("""+path+""");"+...
	 "exec(""Load.sci"");"+...
	 "chdir(oldpath);";
endfunction;

function [LIBPATH, LIBEXT] = LibTools()
  if MSDOS then
    LIBPATH = 'libWindows';
    LIBEXT = '.dll';
  else
    os = unix_g('uname -s');
    if os == 'Linux' then
      LIBPATH = 'libLinux';
      LIBEXT = '.so';
    elseif os == 'Darwin' then
      LIBPATH = 'libDarwin';
      LIBEXT = '.so';
    else
      printf('Operating System not supported\n');
      pause;
    end;
  end;
endfunction;
