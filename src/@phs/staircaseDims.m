function [dims] = staircaseDims(sys)
% STAIRCASEDIMS - Computes the partitioning of the state vector x = [x1;x2;x3;x4] in R^(n1+n2+n3+n4) 
%                 of a pHDAE in staircase form
%
% Syntax:
%   dims = staircaseDims(sys)
%
% Description:
%
% The function computes the partitioning x = [x1;x2;x3;x4] in R^(n1+n2+n3+n4)               
% of the pHDAE system
%     E*dx/dt = (J-R)*x(t)  + (G-P)*u(t),
%           y = (G+P)'*x(t) + (S+N)*u(t),
%
% in staircase form with
%
%      [E11  0  0 0]      [J11 J12 J13 J14]      [R11 R12 R13  0 ]      [G1]      [P1]
%  E = [ 0  E22 0 0], J = [J21 J22 J23  0 ], R = [R21 R22 R23  0 ], G = [G2], P = [P2]
%      [ 0   0  0 0]      [J31 J32 J33  0 ]      [R31 R32 R33  0 ]      [G3]      [P3]
%      [ 0   0  0 0]      [J41  0   0   0 ]      [ 0   0   0   0 ]      [G4]      [0 ]
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:      phs object 
%
% Output Arguments:
%       - dims: vector with [n1,n2,n3,n4]
%
% See Also:
%       staircaseCheck, toStaircase
%
% References:
%       [1] F. Achleitner, A. Arnold, and V. Mehrmann. Hypocoercivity and controllability in linear
%           semi-dissipative ODEs and DAEs. ZAMM Z. Angew. Math. Mech., 2021.
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

%% Input parsing
narginchk(1,1);
% Check phs input type
if ~isa(sys,'phs')
    error('phs:staircaseCheck:wrongInput', 'Model is not an object of the phs-class.');
end

%% Compute n1+n2
zRows_E = all(sys.E==0,2);
if all(zRows_E)
    n12 = 0;
else
    n12 = find(~zRows_E,1,'last');
end

%% Compute n3
JR = sys.J-sys.R;
JRpart = JR(n12+1:end,n12+1:end);
zRows_JRpart = all(JRpart==0,2);
if all(zRows_JRpart)
    n3 = 0;
else
    n3 = find(~zRows_JRpart,1,'last');
end

%% Compute n1, n2, n4
n = sys.dim;
n4 = n - n3 - n12;
n1 = n4;
n2 = n12 - n1;

dims = [n1,n2,n3,n4];
