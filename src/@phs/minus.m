function errSys = minus(sys1, sys2)
% MINUS - Returns the system model (ss/sss object) of the system which
%               subtracts the output of sys2 from the output of sys1.
%
% Syntax:
%   errSys = MINUS(sys1, sys2)
%
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys1: System 1, "positive", phs object
%       - sys2: System 2, "negative", phs object
%
% Output Arguments:
%       - errSys:   error system, ss(s) object
%
%
%-----------------------------------------------------------------------
% This file is part of 
%
% <a href="https://github.com/MORLab/MORpH">MORpH</a> - a MATLAB toolbox to store, analyze,
% interconnect and reduce large-scale port-Hamiltonian models
%
% Authors:      Julius Durmann, Tim Moser
% E-Mail:       <a href="mailto:morlab.rt@ed.tum.de">morlab.rt@ed.tum.de</a>
% Website:      <a href="https://www.epc.ed.tum.de/en/rt/home">www.epc.ed.tum.de/rt</a> 
% Copyright :   © 2022 Chair of Automatic Control, TUM (see LICENSE.md)
%-----------------------------------------------------------------------

if isa(sys2,'sss')
    % Convert to sss
    sys1 = sss(sys1);
    % Create error system
    errSys = sys1 - sys2;
elseif isa(sys2,'ss') || isa(sys2,'phs')
    % Convert to ss
    sys1 = ss(sys1);
    sys2 = ss(sys2);
    % Create error system
    errSys = sys1 - sys2;
else
    error('phs:minus:wrongInput',...
        'phs objects may be subtracted from other phs, ss and sss objects.');
end