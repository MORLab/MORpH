function [sys] = setup_RandomStaircasePHDAE(dimsX,nInp)
% SETUP_RANDOMSTAIRCASEPHDAE - Creates a random pHDAE in staircase form
%
% Syntax:
%   sys = SETUP_RANDOMSTAIRCASEPHDAE(dimsX, nInp)
%
% Input Arguments:
%       *Required Input Arguments:*
%       - dimsX:	1x4 vector [n1,n2,n3,n4] with dimensions of partitioned state
%                   vector x = [x1;x2;x3;x4] in R^{n1+n2+n3+n4}
%       - nInp:     number of inputs and outputs
%
% Output Arguments:
%       - sys: 	phDAE in staircase form (phs-object)
%
% See Also:
%       staircaseCheck, staircaseDims, toStaircase
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

% Input parsing
if any(size(dimsX)~=[1,4]) || any(floor(dimsX)~=dimsX) || any(dimsX<0) || dimsX(1) ~= dimsX(4)
    error('MORpH:setup_RandomStaircasePHDAE:wrongInput', 'Please enter an 1x4 vector with non-negative integers and n1 == n4 for dimsX.');
end
if any(size(nInp)~=[1,1]) || floor(nInp)~=nInp || nInp<0
    error('MORpH:setup_RandomStaircasePHDAE:wrongInput', 'Number of inputs must be a non-negative integer.');
end

n = sum(dimsX);
n1 = dimsX(1); n2 = dimsX(2); n3 = dimsX(3); n4 = dimsX(4);

% Q
Q = eye(n);

% E
eta1 = -1+2*rand(n1); eta2 = -1+2*rand(n2);
E = blkdiag(eta1'*eta1,eta2'*eta2,zeros(n3+n4));

% J
eta = -1+2*rand(n); J = eta-eta';
J(n-n4+1:end,n1+1:end) = zeros(n4,n2+n3+n4);
J(n1+1:end,n-n4+1:end) = zeros(n2+n3+n4,n4);

% G
G = -1+2*rand(n,nInp);

% N
eta = -1+2*rand(nInp,nInp); N = eta-eta';

% W
eta = -1+2*rand(floor(1+n*rand),n1+n2+n3+nInp);
W = eta'*eta;

% R, P, S
R = W(1:(n1+n2+n3),1:(n1+n2+n3)); R = blkdiag(R,zeros(n4));
P = W(1:(n1+n2+n3),end-nInp+1:end); P = [P;zeros(n4,nInp)];
S = W(end-nInp+1:end,end-nInp+1:end);

% Create system
sys = phs(J,R,Q,G,E,P,S,N);

end

