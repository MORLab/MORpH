function [varargout] = step(varargin)
% STEP - Plots and returns the step response of the port-Hamiltonian system
%        sys
%
% Syntax:
%   STEP(sys)
%   STEP(sys,tfinal)
%   STEP(sys,t)
%   STEP(sys1,sys2,...)
%   [y, t] = STEP(sys)
%   [y, t] = STEP(sys, tfinal)
%   ... (see DynamicSystem/step)
%
% Description:
%       phs/step is a wrapper for DynamicSystem/step!
%
%       step(sys) plots the step-response of the phs-object sys.
%       In case of MIMO-systems only the first output stimulated by a step
%       in the first input is plotted.
%
%       [y, t] = step(sys, tfinal) plots and returns the step-response of
%       the system sys calculated to tfinal.
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:		phs object
%       *Optional Input Arguments:*
%       - tfinal:   end time of simulation
%       - t:        time vector
%       - *:        see DynamicSystem/step
%
% Output Arguments:
%       - y:        vector/matrix of step-response
%       - t:        vector of time corresponding to y
%       - *:        See DynamicSytem/step
%
% Examples:
%       The following code creates a simple phs object and runs the
%       step-function:
%
%       J = [0 -1; 1 0]; R = [0 0; 0 1]; Q = eye(2); G = [1; 0];
%       sys = phs(J, R, Q, G);
%       step(sys);
%
% See Also:
%       phs, DynamicSystem/step
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

% Transform phs (phsRed) objects to ss-objects
for i = 1:length(varargin)
    if isa(varargin{i},'phs') || isa(varargin{i},'phsRed')
        varargin{i} = ss(varargin{i});
    end
end

% Call DynamicSystem/step
varargout = cell(1,nargout);
[varargout{:}] = step(varargin{:});

end