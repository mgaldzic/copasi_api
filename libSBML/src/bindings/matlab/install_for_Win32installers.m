%  Filename: install_Win32.m
%  
%  installs the MATLAB binding on a Windows platform

%
%  Description : File to install SBMLToolbox
%  Author(s)   : SBML Team <sbml-team@caltech.edu>
%  Organization: University of Hertfordshire STRC
%  Created     : 2003-10-01
%  Revision    : $Id: install_for_Win32installers.m 11627 2010-08-02 16:24:47Z sarahkeating $
%  $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/bindings/matlab/install_for_Win32installers.m $
%
%  Copyright 2003 California Institute of Technology, the Japan Science
%  and Technology Corporation, and the University of Hertfordshire
%
%  This library is free software; you can redistribute it and/or modify it
%  under the terms of the GNU Lesser General Public License as published
%  by the Free Software Foundation; either version 2.1 of the License, or
%  any later version.
%
%  This library is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY, WITHOUT EVEN THE IMPLIED WARRANTY OF
%  MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  The software and
%  documentation provided hereunder is on an "as is" basis, and the
%  California Institute of Technology, the Japan Science and Technology
%  Corporation, and the University of Hertfordshire have no obligations to
%  provide maintenance, support, updates, enhancements or modifications.  In
%  no event shall the California Institute of Technology, the Japan Science
%  and Technology Corporation or the University of Hertfordshire be liable
%  to any party for direct, indirect, special, incidental or consequential
%  damages, including lost profits, arising out of the use of this software
%  and its documentation, even if the California Institute of Technology
%  and/or Japan Science and Technology Corporation and/or University of
%  Hertfordshire have been advised of the possibility of such damage.  See
%  the GNU Lesser General Public License for more details.
%
%  You should have received a copy of the GNU Lesser General Public License
%  along with this library; if not, write to the Free Software Foundation,
%  Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
%
%  The original code contained here was initially developed by:
%
%      Sarah Keating
%      Science and Technology Research Centre
%      University of Hertfordshire
%      Hatfield, AL10 9AB
%      United Kingdom
%
%      http://www.sbml.org
%      mailto:sbml-team@caltech.edu
%
%  Contributor(s):
%
%

if (strcmp(isoctave(), '0'))
  matlab = 1;
else
  matlab = 0;
end;
% add the current directory to the Matlab search
% path and save
addpath(pwd);

% path2rc is deprecated by version 7.0.4 
% replaced by savepath
% but savepath doesnt exist in version 6.5.1 or lower

if (matlab)
  v = version;
  v_num = str2num(v(1));

  if (v_num < 7)
    saved = path2rc;
  else
    saved = savepath;
  end;
else
  saved = savepath;
end;

if (saved ~= 0)
    error('Directory NOT added to the path');
end;

% try the executable
% if it doesnt work the library files are not on the system path and need
% to be placed there
try
    M = TranslateSBML('test.xml');
catch
    % determine the matlabroot for windows executable
    % this directory is saved to the environmental variable PATH
    Path_to_libs = matlabroot;

    if (matlab)
      Path_to_libs = strcat(Path_to_libs, '\bin\win32');
      % determine the location of the library files
      lib{1} = '..\..\win32\lib\libsbml.lib';
      lib{2} = '..\..\win32\bin\libsbml.dll';
      lib{3} = '..\..\win32\lib\libxml2.lib';
      lib{4} = '..\..\win32\bin\libxml2.dll';
      lib{5} = '..\..\win32\lib\iconv.lib';
      lib{6} = '..\..\win32\bin\iconv.dll';

      for i = 1:6
        copyfile(lib{i}, Path_to_libs);
      end;
    else
      Path_to_libs = strcat(Path_to_libs, '\bin');

      % determine the location of the library files
      lib{1} = '..\..\win32\lib\libsbml.lib';
      lib{2} = '..\..\win32\bin\libsbml.dll';

      for i = 1:2
        copyfile(lib{i}, Path_to_libs);
      end;
    end;
end;

try
  M = TranslateSBML('test.xml');
  disp('Installation successful');
catch
  disp('Installation failed.');
end;
