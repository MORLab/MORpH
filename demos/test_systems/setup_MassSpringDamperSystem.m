function sys = setup_MassSpringDamperSystem(n,m,k,c,type)
% SETUP_MASSSPRINGDAMPERSYSTEM - This function creates a port-Hamiltonian
% model for a mechanical system of connected masses, springs and dampers (see [1],[2])
%
% Syntax:
%   sys = SETUP_MASSSPRINGDAMPERSYSTEM(n,m,k,c,type)
%
% Input Arguments:
%   - n:        Dimension n of the system
%   - m, k, c:  These parameters change the mass, spring constant and
%               damping coefficient respecitvely, can be scalars or vectors
%               with length n/2
%   - type:     {"SISO", "MIMO"}: Request SISO or MIMO model.
%
% Output Arguments:
%   - sys:      port-Hamiltonian model (phs-object)
%
% See Also:
%       setup_LadderNetworkSystem
%
%   References:
%   [1] S. Gugercin, R. V. Polyuga, C. Beattie, and A. van der Schaft. Interpolation-based H2
%       model reduction for port-Hamiltonian systems. In Proceedings 48th IEEE Conference on Decision
%       and Control, pp. 5362–5369, 2009.
%   [2] R. V. Polyuga and A. van der Schaft. Effort- and flow-constraint reduction methods for structure
%       preserving model reduction of port-Hamiltonian systems. Systems Control Lett., 61(3):412–421, 2012.
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

assert(mod(n,2) == 0, ...
    'The function generates port-Hamiltonian representations for even-dimensional Ladder Networks')
narginchk(4,5)

% Transform scalar inputs to vectors
if isscalar(m)
    m = ones(n/2,1)*m;
end
if isscalar(k)
    k = ones(n/2,1)*k;
end
if isscalar(c)
    c = ones(n/2,1)*c;
end

% Transform to row vectors
if size(m,2) == 1
    m = m';
end
if size(k,2) == 1
    k = k';
end
if size(c,2) == 1
    c = c';
end

%% Q
Q = zeros(n);
s = size(Q);
% Diagonal elements
idx = sub2ind(s,1:2:s(1)-1,1:2:s(2)-1);
Q(idx) = Q(idx) + k(1:end);
idx = sub2ind(s,3:2:s(1)-1,3:2:s(2)-1);
Q(idx) = Q(idx) + k(1:end-1);
idx = sub2ind(s,2:2:s(1),2:2:s(2));
Q(idx) = Q(idx) + 1./m;
% Super-/sub-diagonal elements
idx = sub2ind(s,1:2:s(1)-3,3:2:s(2)-1);
Q(idx) =  - k(1:end-1);
idx = sub2ind(s,3:2:s(1)-1,1:2:s(2)-3);
Q(idx) = - k(1:end-1);

%% R
R = zeros(n);
s = size(R);
idx = sub2ind(s,2:2:s(1),2:2:s(2));
R(idx) = c;

%% J
J = zeros(n);
s = size(J);
idx = sub2ind(s,2:2:s(1),1:2:s(2)-1);
J(idx) = -1;
idx = sub2ind(s,1:2:s(1)-1,2:2:s(2));
J(idx) = 1;

%% G
if nargin == 4 || strcmp(type, 'SISO')
    % SISO
    G = zeros(n,1);
    G(2, 1) = 1;
elseif strcmp(type,'MIMO')
    % MIMO
    G = zeros(n,2);
    G(2, 1) = 1;
    G(4, 2) = 1;
else
    % SISO
    warning("Unknown keyword for type (must be ""SISO"" or ""MIMO""). Creating SISO model...")
    % SISO
    G = zeros(n,1);
    G(2, 1) = 1;
end

%% Create system
Opts.inputValidation = false;
Opts.verbose = false;   % Suppress input validation warning
sys = phs(J,R,Q,G,Opts);
sys.Opts.verbose = true;    % Turn verbose mode on again

end

