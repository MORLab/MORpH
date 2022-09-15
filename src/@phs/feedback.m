function [fbSys] = feedback(sys_plant, sys_contr)
% FEEDBACK - returns a system which consists of two original systems sys_plant
%       and sys_contr in a negative feedback loop
%
% Description:
%   u ------ o ----> [ sys_plant ] -----+--------> y
%            ^ -                        |
%            |                          |
%             ------ [ sys_contr ] <----
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys_plant:	open loop system (phs-object)
%       - sys_contr: 	control system (phs-object)
%
% Output Arguments:
%       - fbSys:  closed loop system (phs-object)
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

% Create surrogate variables to simplify equations
J1 = sys_plant.J;   J2 = sys_contr.J;
R1 = sys_plant.R;   R2 = sys_contr.R;
Q1 = sys_plant.Q;   Q2 = sys_contr.Q;
G1 = sys_plant.G;   G2 = sys_contr.G;
E1 = sys_plant.E;   E2 = sys_contr.E;
P1 = sys_plant.P;   P2 = sys_contr.P;
S1 = sys_plant.S;   S2 = sys_contr.S;
N1 = sys_plant.N;   N2 = sys_contr.N;
D1 = S1 + N1;       D2 = S2 + N2;
m = size(G1,2);

% Check some conditions
assert(isequal(size(G1,2),size(G2,2)),'phs:feedback:inputDimensionsDontMatch',...
    'Input dimensions of sys_plant and sys_contr must be equal.');

assert(rank(eye(m)+D2*D1) == m,'phs:feedback:algebraicLoop',...
    'The matrix [I + (sys_contr.S+sys_contr.N)*(sys_plant.S+sys_plant.N)] must be invertible.');

% Compute general state-space matrices
K = eye(m) + D2*D1;
JR = [J1-R1 - (G1-P1)*(K\D2)*(G1+P1)'     , -(G1-P1)*(K\(G2+P2)'); ...
    (G2-P2)*(eye(m)-D1*(K\D2))*(G1+P1)' , J2-R2 - (G2-P2)*D1*(K\(G2+P2)')];
B = [(G1-P1)*(eye(m)-(K\D2)*D1)  ; (G2-P2)*(eye(m)-D1*(K\D2))*D1];
C = [(eye(m)-D1*(K\D2))*(G1+P1)' , -D1*(K\(G2+P2)')];
SN = (eye(m)-D1*(K\D2))*D1;

% Create pH system
E = blkdiag(E1,E2);
J = 0.5*(JR - JR');
R = -0.5*(JR + JR');
Q = blkdiag(Q1,Q2);
G = 0.5*(C'+B);
P = 0.5*(C'-B);
S = 0.5*(SN + SN');
N = 0.5*(SN - SN');

fbSys = phs(J,R,Q,G,E,P,S,N);

end