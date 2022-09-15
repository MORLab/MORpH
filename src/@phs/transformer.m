function [sysTr] = transformer(sys1, sys2, M)
% TRANSFORMER - Returns a system which consists of two original systems sys1
%       and sys2 with external input/output pairs and internal input/output
%       pairs connected by a transformer interconnection
%           u1_int = -M*u2_int
%           y2_int = M'*y1_int
%
% Description:
%
%         u1_ext ----> [          ] ----> y1_ext
%                      [   sys1   ]
%          u1_int ---- [          ] ---- y1_int
%                 |                    |
%        u2_int [-M]-- [          ] --[M'] y2_int
%                      [   sys2   ] 
%         u2_ext ----> [          ] ----> y2_ext
% 
%                                   
% Input Arguments:
%       *Required Input Arguments:*
%       - sys1:	system one (phs object)
%       - sys2: system two (phs object)
% 
% Output Arguments:
%       - sysTr:  interconnected system (phs object)
%
% See Also:
%      	gyrator, feedback
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

% Get matrices
[J1, R1, Q1, G1, E1, P1, S1, N1] = getMatrices(sys1);
[J2, R2, Q2, G2, E2, P2, S2, N2] = getMatrices(sys2);

% Check coupling dimensions
nInt = size(M,1);

assert(nInt <= min(size(G1,2),size(G2,2)) && ...
    nInt == size(M,2) && ...
    isreal(M),...
    'phs:transformer:wrongInput',...
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

% Construct interconnected model
E = blkdiag(E1,E2,zeros(nInt));
Q = blkdiag(Q1,Q2,zeros(nInt));
JR = [blkdiag(J1-R1,J2-R2), [-(G1i-P1i)*M;G2i-P2i];...
     [M'*(G1i+P1i)',-(G2i+P2i)',-M'*D1_11*M-D2_11]];
B = [blkdiag(G1e-P1e,G2e-P2e); M'*D1_12, -D2_12];
C = [blkdiag((G1e+P1e)',(G2e+P2e)'), [-D1_21*M; D2_21]];
D = blkdiag(D1_22,D2_22);

% Decompose in symmetric and skew-symmetric parts
J = 0.5*(JR - JR');
R = -0.5*(JR + JR');
G = 0.5*(C' + B);
P = 0.5*(C' - B);
S = 0.5*(D + D');
N = 0.5*(D - D');

sysTr = phs(J,R,Q,G,E,P,S,N);

end
