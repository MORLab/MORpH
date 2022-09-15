function thirdPartyCheck(software)
% THIRDPARTYCHECK - Checks if a third-party software is installed
%
% Syntax:
%   THIRDPARTYCHECK('Manopt')
%
% Input Arguments:
%       *Required Input Arguments:*
%       - software:	String of third-party software
%                   {'Manopt','M.E.S.S.','GRANSO','CVX','YALMIP','SADPA','SAMDP','LINORM','HINORM'}
%
%-----------------------------------------------------------------------
% This file is part of
%
% <a href="https://github.com/MORLab/MORpH">MORpH</a> - a MATLAB toolbox to store, analyze,
% interconnect and reduce large-scale port-Hamiltonian models
%
% Authors:      Tim Moser
% E-Mail:       <a href="mailto:morlab.rt@ed.tum.de">morlab.rt@ed.tum.de</a>
% Website:      <a href="https://www.epc.ed.tum.de/en/rt/home">www.epc.ed.tum.de/rt</a>
% Copyright :   © 2022 Chair of Automatic Control, TUM (see LICENSE.md)
%-----------------------------------------------------------------------

switch software
    case 'Manopt'
        filename = 'importmanopt.m';
        iLink = '<a href="https://www.manopt.org/downloads.html">https://www.manopt.org/downloads.html</a>';
        dLink = '<a href="https://github.com/NicolasBoumal/manopt">https://github.com/NicolasBoumal/manopt</a>';
    case 'M.E.S.S.'
        filename = 'operatormanager.m';
        iLink =  '<a href="https://www.mpi-magdeburg.mpg.de/projects/mess">https://www.mpi-magdeburg.mpg.de/projects/mess</a>';
        dLink =  '<a href="https://github.com/mpimd-csc/mmess">https://github.com/mpimd-csc/mmess</a>';
    case 'GRANSO'
        filename = 'granso.m';
        iLink =  '<a href="http://www.timmitchell.com/software/GRANSO/">http://www.timmitchell.com/software/GRANSO/</a>';
        dLink =  '<a href="https://gitlab.com/timmitchell/GRANSO/">https://gitlab.com/timmitchell/GRANSO/</a>';
    case 'CVX'
        filename = 'cvx_version.m';
        iLink =  '<a href="http://cvxr.com/cvx/download/">http://cvxr.com/cvx/download/</a>';
        dLink =  '<a href="https://github.com/cvxr/CVX">https://github.com/cvxr/CVX</a>';
    case 'YALMIP'
        filename = 'yalmip.m';
        iLink =  '<a href="https://yalmip.github.io/download/">https://yalmip.github.io/download/</a>';
        dLink =  '<a href="https://github.com/yalmip/YALMIP">https://github.com/yalmip/YALMIP</a>';
    case 'SADPA'
        filename = 'sadpa.m';
        iLink =  '<a href="https://sites.google.com/site/rommes/software">https://sites.google.com/site/rommes/software</a>';
        dLink =  [];
    case 'SAMDP'
        filename = 'samdp.m';
        iLink =  '<a href="https://sites.google.com/site/rommes/software">https://sites.google.com/site/rommes/software</a>';
        dLink =  [];
    case 'LINORM'
        filename = 'linorm_subsp.m';
        iLink =  '<a href="http://www.math.tu-berlin.de/index.php?id=186267&L=1">http://www.math.tu-berlin.de/index.php?id=186267&L=1</a>';
        dLink =  [];
    case 'HINORM'
        filename = 'hinorm.m';
        iLink =  '<a href="http://www.mpi-magdeburg.mpg.de/mpcsc/software/infnorm">http://www.mpi-magdeburg.mpg.de/mpcsc/software/infnorm</a>';
        dLink =  [];
    otherwise
        error(strcat("Unknown third-party software ",software,"."));
end

% Create error message
errMsg = strcat("With the specified configuration, this function requires the software ",software,".",{newline},...
    "You can download the latest release from ",{newline},iLink);
if ~isempty(dLink)
    errMsg = strcat(errMsg,",",{newline},"or clone the repository from",{newline},dLink,".",{newline});
else
    errMsg = strcat(errMsg,".",{newline});
end
errMsg = strcat(errMsg,"Please follow the installation instructions and add the software to your MATLAB path.");

% Throw error if filename is not found
assert( exist( filename, 'file' ) == 2, string(errMsg));

end