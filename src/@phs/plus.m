function sys = plus(sys1,sys2)
% PLUS - Returns a system which consists of two original systems sys 1 and
%           sys2 in parallel connection
%
% Description:
%            ------> [ sys1 ] -------
%            |                      |
%            |                      v
%   u ------ +                      o -------> y
%            |                      ^
%            |                      |
%            ------> [ sys2 ] -------
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys1:	first system (phs object)
%       - sys2: second system (phs object)
%
% Output Arguments:
%       - sys:  sys1+sys2 (parallel) (phs object)
%
% See Also:
%      	phs
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

if isa(sys2,'phs')
    % Extract system matrices
    [J1, R1, Q1, G1, E1, P1, S1, N1] = getMatrices(sys1);
    [J2, R2, Q2, G2, E2, P2, S2, N2] = getMatrices(sys2);

    % Add systems
    J = blkdiag(J1,J2);
    R = blkdiag(R1,R2);
    Q = blkdiag(Q1,Q2);
    G = [G1;G2];
    E = blkdiag(E1,E2);
    P = [P1;P2];
    S = S1+S2;
    N = N1+N2;

    sys = phs(J,R,Q,G,E,P,S,N);

elseif isa(sys2,'ss')
    % Convert to ss
    sys1 = ss(sys1);
    sys2 = ss(sys2);

    sys = sys1 + sys2;

elseif isa(sys2,'sss')
    % Convert to sss
    sys1 = sss(sys1);
    sys2 = sss(sys2);

    sys = sys1 + sys2;

else
    error('phs:plus:wrongInput',...
        'phs objects may be added to other phs, ss and sss objects.');
end