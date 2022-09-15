function [sysGyr] = gyrator(sys1, sys2, M)
% GYRATOR - returns a system which consists of two original systems sys1
%       and sys2 with external input/output pairs and internal input/output
%       pairs connected by a gyrator interconnection
%           u1_int = -M*y2_int
%           u2_int = M'*y1_int
%
% Description:
%
%         u1_ext ----> [          ] ----> y1_ext
%                      [   sys1   ]
%          u1_int ---- [          ] ------- y1_int
%               __|_______________________|
%              |  |_______________________
%              |                          |
%       u2_int [M']--- [          ] ----[-M] y2_int
%                      [   sys2   ]
%         u2_ext ----> [          ] ----> y2_ext
%
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys1:	system one (phs-object)
%       - sys2: system two (phs-object)
%
% Output Arguments:
%       - sysGyr:  interconnected system (phs-object)
%
% See Also:
%      	feedback, transformer
%
%-----------------------------------------------------------------------
% This file is part of
%
% <a href="https://github.com/MORLab/MORpH">MORpH</a> - a MATLAB toolbox to store, analyze,
% interconnect and reduce large-scale port-Hamiltonian models
%
% Authors:      Tim Moser
% E-Mail:       <a href="mailto:morlab.rt@ed.tum.de">morlab.rt@ed.tum.de</a>
% Website:      <a href="https://www.epc.ed.tum.de/en/rt/home">www.epc.ed.tum.de/rt</a>
% Copyright :   © 2022 Chair of Automatic Control, TUM (see LICENSE.md)
%-----------------------------------------------------------------------

% Get matrices
[J1, R1, Q1, G1, E1, P1, S1, N1] = getMatrices(sys1);
[J2, R2, Q2, G2, E2, P2, S2, N2] = getMatrices(sys2);

% Parse coupling matrix
nInt = size(M,1);

assert(nInt <= min(size(G1,2),size(G2,2)) && ...
    nInt == size(M,2) && ...
    isreal(M),...
    'phs:gyrator:wrongInput',...
    ['The coupling matrix M must be real, square and have dimensions ' ...
    'smaller than or equal to the number of inputs for both systems.']);

% Split inputs and ouputs
G1i = G1(:,1:nInt); G1e = G1(:,nInt+1:end);
P1i = P1(:,1:nInt); P1e = P1(:,nInt+1:end);
G2i = G2(:,1:nInt); G2e = G2(:,nInt+1:end);
P2i = P2(:,1:nInt); P2e = P2(:,nInt+1:end);

D1 = S1+N1; D2 = S2+N2;
D1_11 = D1(1:nInt,1:nInt);
D1_12 = D1(1:nInt,nInt+1:end);
D1_21 = D1(nInt+1:end,1:nInt);
D1_22 = D1(nInt+1:end,nInt+1:end);

D2_11 = D2(1:nInt,1:nInt);
D2_12 = D2(1:nInt,nInt+1:end);
D2_21 = D2(nInt+1:end,1:nInt);
D2_22 = D2(nInt+1:end,nInt+1:end);

% Check for algebraic loops
assert(rank(eye(nInt)+M*D2_11*M'*D1_11) == nInt,'phs:gyrator:algebraicLoop',...
    "The matrix [I + M*D2(1,1)*M'*D1(1,1)] must be invertible.");
assert(rank(eye(nInt)+M'*D1_11*M*D2_11) == nInt,'phs:gyrator:algebraicLoop',...
    "The matrix [I + M'*D1(1,1)*M*D2(1,1)] must be invertible.");

Xi1 = -(eye(nInt)+M*D2_11*M'*D1_11)\M;
Xi2 = (eye(nInt)+M'*D1_11*M*D2_11)\M';

% Construct interconnected model
E = blkdiag(E1,E2);
Q = blkdiag(Q1,Q2);
JR = [J1-R1 + (G1i-P1i)*Xi1*D2_11*M'*(G1i+P1i)' , (G1i-P1i)*Xi1*(G2i+P2i)';...
    (G2i-P2i)*Xi2*(G1i+P1i)'                  , J2-R2 - (G2i-P2i)*Xi2*D1_11*M*(G2i+P2i)'];
B = [(G1e-P1e)+(G1i-P1i)*Xi1*D2_11*M'*D1_12 , (G1i-P1i)*Xi1*D2_12; ...
    (G2i-P2i)*Xi2*D1_12                    , (G2e-P2e)-(G2i-P2i)*Xi2*D1_11*M*D2_12];
C = [(G1e+P1e)'+D1_21*Xi1*D2_11*M'*(G1i+P1i)' , D1_21*Xi1*(G2i+P2i)';...
    D2_21*Xi2*(G1i+P1i)'                     , (G2e+P2e)'-D2_21*Xi2*D1_11*M*(G2i+P2i)'];
D = [D1_21*Xi1*D2_11*M'*D1_12 + D1_22 , D1_21*Xi1*D2_12;...
    D2_21*Xi2*D1_12                  , D2_22 - D2_21*Xi2*D1_11*M*D2_12];

% Decompose in symmetric and skew-symmetric parts
J = 0.5*(JR - JR');
R = -0.5*(JR + JR');
G = 0.5*(C' + B);
P = 0.5*(C' - B);
S = 0.5*(D + D');
N = 0.5*(D - D');

sysGyr = phs(J,R,Q,G,E,P,S,N);

end
