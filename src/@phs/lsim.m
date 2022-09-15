function [varargout] = lsim(varargin)
% LSIM - simulates the port-Hamiltonian system sys for given input u and
%        time t
%
% Syntax:
%   y = LSIM(sys, u, t)
%   ... (see DynamicSystem/lsim)
%
% Description:
%       phs/lsim is a wrapper for DynamicSystem/lsim!
%
%       y = lsim(sys, u, t) simulates the system for input u
%       and returns the system response y corresponding to the time
%       vector t
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:		phs object
%       - u:        Input vector
%       - t:        Time vector
%
% Output Arguments:
%       - y: 		System response
%
% Examples:
%       The following code creates a simple phs object and runs the
%       lsim-function for a sine input signal:
%
%       J = [0 -1; 1 0]; R = [0 0; 0 1]; Q = eye(2); G = [1; 0];
%       sys = phs(J, R, Q, G);
%       t = linspace(0, 100, 1000);
%       u = sin(2*pi*0.1*t);
%       lsim(sys, u, t);
%
% See Also:
%       phs, DynamicSystem/lsim
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

% Call DynamicSystem/lsim
varargout = cell(1,nargout);
[varargout{:}] = lsim(varargin{:});

end