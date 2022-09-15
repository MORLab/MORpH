function sys = setup_LadderNetworkSystem(n,L,C,Res)
% SETUP_LADDERNETWORKSYSTEM - This function creates a port-Hamiltonian model for an electrical
%   system of connected coils, capacitors and resistors (see [1],[2])
%
% Syntax:
%   sys = SETUP_LADDERNETWORKSYSTEM(n,L,C,Res)
%
% Input Arguments:
%       *Required Input Arguments:*
%       - n:	dimension of the system
%       - L:    Inductance of all coils
%       - C:    Capacitance of all capacitors
%       - Res:  Resistance of all resistors
%
% Output Arguments:
%       - sys: 	port-Hamiltonian model (phs-object)
%
% See Also:
%       setup_MassSpringDamperSystem
%
%   References:
%       [1] S. Gugercin, R. Polyuga, C. Beattie, and A. van der Schaft, 
%           "Interpolation-based H2 model reduction for portHamiltonian systems,” 
%           in Proceedings of the 48h IEEE Conference on Decision and Control (CDC), 
%           2009, pp. 5362–5369. 
%       [2] Polyuga (2010), "Model Reduction of Port-Hamiltonian Systems", Dissertation, 
%           Figure 5.1
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

if mod(n,2) ~= 0
    error('This system must have an even order');
end

% J
J = diag(ones(n-1,1), 1) + diag(-1*ones(n-1,1), -1);

% R
if length(Res) == 1
    Res = repmat(Res,n/2+1,1);
end
diagR = zeros(n,1);
index_diagR = 2:2:n;
diagR(index_diagR) = Res(1:end-1);
diagR(end) = diagR(end) + Res(end);
R = diag(diagR);

% Q
diagQ = zeros(n,1);
index_diagQ_C = 1:2:n-1;
index_diagQ_L = 2:2:n;
diagQ(index_diagQ_C) = 1./C;
diagQ(index_diagQ_L) = 1./L;
Q = diag(diagQ);

% G
G = zeros(n,1);
G(1) = 1;

% Create system
sys = phs(J,R,Q,G);

end