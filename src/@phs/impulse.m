function [varargout] = impulse(varargin)
% IMPULSE - plots and returns the impulse response of the
%           port-Hamiltonian system sys
%
% Syntax:
%   IMPULSE(sys)
%   IMPULSE(sys, tfinal)
%   IMPULSE(sys, t)
%   [y, t] = IMPULSE(sys)
%   [y, t] = IMPULSE(sys, tfinal)
%   ... (see DynamicSystem/impulse)
%
% Description:
%       phs/impulse is a wrapper for DynamicSystem/impulse!
%
%       impulse(sys) plots the impulse response of the PH system sys.
%       In case of MIMO-systems, independent impulse commands are
%       applied to each input channel.
%
%       [y, t] = impulse(sys, tfinal) plots the impulse response of the PH
%       system sys and returns the vectors t and y of the system response.
%       It therefore transforms the phs-object into a ss-object and
%       then runs the function ss/impulse.
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:		phs object
%       *Optional Input Arguments:*
%       - tfinal:	end time of simulation
%       - t:        time vector
%
% Output Arguments:
%       - y: 		vector of system response (output-variable y)
%       - t: 		vector of time corresponding to vector y
%
% Examples:
%       The following code creates a simple phs object and runs the
%       impulse-function:
%
%       J = [0 -1; 1 0]; R = [0 0; 0 1]; Q = eye(2); G = [1; 0];
%       sys = phs(J, R, Q, G);
%       impulse(sys);
%
% See Also:
%       phs, ss/impulse
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

% Call DynamicSystem/impulse
varargout = cell(1,nargout);
[varargout{:}] = impulse(varargin{:});

end