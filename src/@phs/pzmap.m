function [varargout] = pzmap(varargin)
% PZMAP - plots and returns the poles and zeros of the port-Hamiltonian
%         system sys
%
% Syntax:
%   PZMAP(sys)
%   PZMAP(sys1,sys2,...)
%   [p, z] = PZMAP(sys)
%
% Description:
%       phs/pzmap is a wrapper for DynamicSystem/pzmap
%
%       pzmap(sys) creates a plot of the poles and zeros of the PH system
%       sys
%
%       [p, z] = pzmap(sys) returns and plots the poles and zeros of the PH
%       system sys
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:		phs object
%
% Output Arguments:
%       - p:        poles of sys
%       - z:        zeros of sys
%
% Examples:
%       The following code creates a simple phs object and runs the
%       pzmap-function:
%
%       J = [0 -1; 1 0]; R = [0 0; 0 1]; Q = eye(2); G = [1; 0];
%       sys = phs(J, R, Q, G);
%       pzmap(sys);
%
% See Also:
%       phs, DynamicSystem/pzmap
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

% Transform phs (phsRed) objects to ss objects
for i = 1:length(varargin)
    if isa(varargin{i},'phs') || isa(varargin{i},'phsRed')
        varargin{i} = ss(varargin{i});
    end
end

% Call DynamicSystem/pzmap
varargout = cell(1,nargout);
[varargout{:}] = pzmap(varargin{:});

end