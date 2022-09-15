function [sysP, Pol_1] = decomposePHDAE(sys, varargin)
% DECOMPOSEPHDAE - Decomposes a PHDAE into a proper PHODE system and an
% improper polynomial part
%
% Syntax:
%   [sysP, Pol_1] = DECOMPOSEPHDAE(sys)
%   [sysP, Pol_1] = DECOMPOSEPHDAE(sys, Opts)
%
% Description:
%
% The transfer function of a general pHDAE system
%     E*dx/dt = (J-R)*Q*x(t)  + (G-P)*u(t),
%           y = (G+P)'*x(t) + (S+N)*u(t),
%
% can be decomposed into
%     G(s) = Gp(s) + Pol_1*s
%
% where Gp(s) is the transfer function of a pHODE system. This
% decomposition is computed by this function.
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:      pHDAE (phs object)
%       *Optional Input Arguments:*
%       - Opts:  structure with execution parameters
%           - .toStaircase.*: Other options that will be passed to this function;
%                     Please refer to its documentation (doc toStaircase)
%           - .phs.*: Other options that will be passed on to the class 'phs';
%                     Please refer to its documentation (doc phs).
%
% Output Arguments:
%       - sysP:  proper part of system (phs object) with transfer function Gp(s)
%       - Pol_1: Coefficient of the linear polynomial part of G(s)
%
% See Also:
%       toStaircase, decomposePHDAE
%
% References:
%       [1] F. Achleitner, A. Arnold, and V. Mehrmann. “Hypocoercivity and 
%           controllability in linear semi-dissipative Hamiltonian ordinary 
%           differential equations and differential-algebraic equations." 
%           In: ZAMM. Z. Angew. Math. Mech. (2021).
%       [2] T. Moser et al. Structure-Preserving Model Order Reduction for Index 
%           Two Port-Hamiltonian Descriptor Systems. 
%           arXiv Preprint arXiv:2206.03942. 2022. url: https://arxiv.org/abs/2206.03942
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

%% Parse Inputs
narginchk(1,2);
[sys, Opts] = parseInputs(sys, varargin{:});

%% Transform system to staircase form
if ~sys.hasStaircase
    [sys, dims] = toStaircase(sys,Opts.toStaircase);
else
    dims = staircaseDims(sys);
end

%% Bring system to almost Kronecker canonical form
[J, R, ~, G, E, P, S, N] = getMatrices(sys);
A = J-R;

n1 = dims(1); n2 = dims(2); n3 = dims(3); n4 = dims(4); n = sum(dims);
n1_start = 1; n1_end = n1;
n2_start = n1+1; n2_end = n1+n2;
n3_start = n1+n2+1; n3_end = n1+n2+n3;
n4_start = n-n4+1; n4_end = n;

if n1 == 0
    La = eye(n2);
    Za = eye(n2);

    if n3 > 0 % Index 1
        A23 = A(n2_start:n2_end,n3_start:n3_end);
        A32 = A(n3_start:n3_end,n2_start:n2_end);
        A33 = A(n3_start:n3_end,n3_start:n3_end);

        La = [La , -A23/A33 ; zeros(n3,n2) , A33\eye(n3)];
        Za = [Za , zeros(n2,n3); -A33\A32 , eye(n3)];
    end

else % Index 2

    if n3 > 0
        A11 = A(n1_start:n1_end,n1_start:n1_end);
        A12 = A(n1_start:n1_end,n2_start:n2_end);
        A21 = A(n2_start:n2_end,n1_start:n1_end);
        A13 = A(n1_start:n1_end,n3_start:n3_end);
        A31 = A(n3_start:n3_end,n1_start:n1_end);
        A23 = A(n2_start:n2_end,n3_start:n3_end);
        A32 = A(n3_start:n3_end,n2_start:n2_end);
        A33 = A(n3_start:n3_end,n3_start:n3_end);
        A14 = A(n1_start:n1_end,n4_start:n4_end);
        A41 = A(n4_start:n4_end,n1_start:n1_end);

        La = [[eye(n1) , zeros(n1,n2) , -A13/A33 , (-A11+(A13/A33)*A31)/A41] ; ...
            [zeros(n2,n1) , eye(n2) , -A23/A33 , (-A21 + (A23/A33)*A31)/A41]; ...
            [zeros(n3,n1+n2) , A33\eye(n3) , (-A33\A31)/A41]; ...
            [zeros(n4,n1+n2+n3) , -A41\eye(n1)]];

        Za = [[eye(n1) , zeros(n1,n2+n3+n4)] ; ...
            [zeros(n2,n1) , eye(n2) , zeros(n2,n3+n4)]; ...
            [zeros(n3,n1) , -A33\A32 , eye(n3),zeros(n3,n4)]; ...
            [zeros(n4,n1) , A14\(-A12+(A13/A33)*A32) , zeros(n4,n3),A14\eye(n4)]];

    else % n3 == 0
        A11 = A(n1_start:n1_end,n1_start:n1_end);
        A12 = A(n1_start:n1_end,n2_start:n2_end);
        A21 = A(n2_start:n2_end,n1_start:n1_end);
        A14 = A(n1_start:n1_end,n4_start:n4_end);
        A41 = A(n4_start:n4_end,n1_start:n1_end);

        La = [[eye(n1) , zeros(n1,n2) , (-A11)/A41] ; ...
            [zeros(n2,n1) , eye(n2) , (-A21)/A41]; ...
            [zeros(n4,n1+n2) , -A41\eye(n1)]];

        Za = [[eye(n1) , zeros(n1,n2+n4)] ; ...
            [zeros(n2,n1) , eye(n2) , zeros(n2,n4)]; ...
            [zeros(n4,n1) , A14\(-A12) , A14\eye(n4)]];
    end
end

% Define canonical system (not pH)
E_c = La*E*Za;
A_c = La*A*Za;
B_c = La*(G-P);
C_c = (G+P)'*Za;

%% Extract proper system with Gp(s)
Ap = A_c(n2_start:n2_end,n2_start:n2_end);
Bp = B_c(n2_start:n2_end,:);
Cp = C_c(:,n2_start:n2_end);
Ep = E_c(n2_start:n2_end,n2_start:n2_end);
Dp = S + N - C_c(:,n4_start:n4_end)*B_c(n1_start:n1_end,:) + ...
    C_c(:,n1_start:n1_end)*B_c(n4_start:n4_end,:)-C_c(:,n3_start:n3_end)*B_c(n3_start:n3_end,:);

Jp = 0.5*(Ap-Ap');
Rp = -0.5*(Ap+Ap');
Gp = 0.5*(Cp'+Bp);
Pp = 0.5*(Cp'-Bp);
Sp = 0.5*(Dp+Dp');
Np = 0.5*(Dp-Dp');

sysP = phs(Jp,Rp,eye(n2),Gp,Ep,Pp,Sp,Np,Opts.phs);

%% Extract polynomial coefficient Pol_1 of linear improper part Pol_1*s
Pol_1 = full(C_c(:,n4_start:n4_end)*E_c(n1_start:n1_end,n1_start:n1_end)*B_c(n4_start:n4_end,:));
if min(eig(Pol_1)) < 0
    warning('MORpH:decomposePHDAE:decompositionError', 'Decomposition did not yield a positive semidefinite polynomial part.');
end

end

function [sys, Opts] = parseInputs(sys, varargin)
% Check phs input type
if ~isa(sys,'phs')
    error('MORpH:decomposePHDAE:wrongInput', 'Model is not an object of the phs-class.');
end

% Opts
if ~isempty(varargin) && isstruct(varargin{end})
    Opts = varargin{end};
else
    Opts = struct();
end

% Option parsing
OptsAdmissible.toStaircase = struct;
OptsAdmissible.phs = struct;
Opts = phsMOR_parseOpts(Opts,OptsAdmissible);

end