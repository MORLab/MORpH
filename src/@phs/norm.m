function [n, varargout] = norm(varargin)
% NORM - Computes the H2 or H-infinity norm of the system sys
%
% Syntax:
%   n = NORM(sys)
%   n = NORM(sys, 2)
%   n = NORM(sys, inf)
%   n = NORM(sys, inf, TOL)
%   [GPEAK, FPEAK] = norm(sys, inf)
%   [...] = norm(..., Opts)
%
% Description:
%       phs/norm is a wrapper for DynamicSystem/norm!
%
%       n = norm(sys) returns the H2-norm of sys
%
%       n = norm(sys, 2) / n = norm(sys, inf)  returns the H2- or the
%       H-infinity norm respectively
%
%       n = norm(sys, inf, TOL) specifies a relative accuracy TOL for the
%       computed infinity norm (default: TOl = 1e-2)
%
%       [GPEAK,FPEAK] = norm(sys, inf) returns the peak value GPEAK at the
%       corresponding frequency FPEAK
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:      phs object
%       *Optional Input Arguments:*
%       - type: 	specify which norm should be computed
%       - TOL:		tolerance (default: 1e-2)
%
% Output Arguments:
%       - n:                norm
%       - [GPEAK,FPEAK]: 	see 'Description'
%
% Examples:
%       The following code creates a simple phs object and runs the
%       norm-function:
%
%       J = [0 -1; 1 0]; R = [0 0; 0 1]; Q = eye(2); G = [1; 0];
%       sys = phs(J, R, Q, G);
%       norm(sys);
%
% See Also:
%       DynamicSystem/norm, phs
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

% Call DynamicSystem/norm
varargout = cell(1,nargout);
[n, varargout{:}] = norm(varargin{:});

end