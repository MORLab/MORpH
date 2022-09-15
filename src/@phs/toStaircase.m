function [sysStair, dims] = toStaircase(sys, varargin)
% TOSTAIRCASE - Transforms a pHDAE into staircase form
%
% Syntax:
%   sysStair = TOSTAIRCASE(sys)
%   sysStair = TOSTAIRCASE(sys, Opts)
%
% Description:
%
% The function transforms a general pHDAE system
%     E*dx/dt = (J-R)*Q*x(t)  + (G-P)*u(t),
%           y = (G+P)'*x(t) + (S+N)*u(t),
%
% into staircase form with partitioned state vector x = [x1;x2;x3;x4] in
% R^(n1+n2+n3+n4) and such that
%
% (i)   The system satisfies the pH structural constraints.
% (ii)  The system matrices have the following zero-patterns:
%      [E11  0  0 0]      [J11 J12 J13 J14]      [R11 R12 R13  0 ]      [G1]      [P1]
%  E = [ 0  E22 0 0], J = [J21 J22 J23  0 ], R = [R21 R22 R23  0 ], G = [G2], P = [P2],
%      [ 0   0  0 0]      [J31 J32 J33  0 ]      [R31 R32 R33  0 ]      [G3]      [P3]
%      [ 0   0  0 0]      [J41  0   0   0 ]      [ 0   0   0   0 ]      [G4]      [0 ]
%
% (iii) rank(E11) = n1, rank(E22) = n2
% (iv)  rank(J14) = n1 = n4
% (v)   rank(J33-R33) = n3
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:      phs object
%       *Optional Input Arguments:*
%       - Opts:
%             - rankTol:    Tolerance for rank decisions
%                           [{1e-12} / positive double]
%             - .scaling.*: Options passed to function 'scaling' 
%                           Please refer to the function's documentation
%                           for more information.
%             - .phs.*:     Options passed to class 'phs'   
%                           Please refer to the function's documentation
%                           for more information.
%
% Output Arguments:
%       - sysStair: pHDAE in staircase form
%       - dims: Splitting of the state vector [n1,n2,n3,n4]
%
% See Also:
%       staircaseDims, staircaseCheck, scaling
%
% References:
%       [1]  F. Achleitner, A. Arnold, and V. Mehrmann. Hypocoercivity and controllability in linear
%            semi-dissipative ODEs and DAEs. ZAMM Z. Angew. Math. Mech., In Press, 2021.
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

%% Scale system (if necessary)
sys = scaling(sys,Opts.scaling);

E = sys.E;
J = sys.J;
R = sys.R;
G = sys.G;
P = sys.P;
S = sys.S;
N = sys.N;

n = size(E,1);

%% Transform E to semi-explicit form
zRows_E = all(sys.E==0,2);
if all(zRows_E)
    n12 = 0;
else
    n12 = find(~zRows_E,1,'last');

    if isequal(E-blkdiag(E(1:n12,1:n12),zeros(n-n12)), zeros(size(E))) ...
            && svds(E(1:n12,1:n12), 1, 'smallest') > Opts.rankTol
        % E already has the required form
    else
        % Eigendecomposition of E
        [Ue,De] = eig(full(E));
        [~,ind]=sort(diag(De),'descend');
        De = De(ind,ind);
        Ue = Ue(:,ind);
        n12 = rank(De,Opts.rankTol);
        T1 = Ue;

        % Transform system
        E = T1'*E*T1;
        J = T1'*J*T1;
        R = T1'*R*T1;
        G = T1'*G;
        P = T1'*P;
    end
end

%% Transform J(3:4,3:4)-R(3:4,3:4) to semi-explicit form
if n-n12 > 0
    % Analyze zero pattern of J-R
    JR34 = J(n12+1:end, n12+1:end) - R(n12+1:end, n12+1:end);
    zRows_JRpart = all(JR34==0,2);
    if all(zRows_JRpart)
        n3 = 0;
    else
        n3 = find(~zRows_JRpart,1,'last');

        if isequal(JR34-blkdiag(JR34(1:n3,1:n3),zeros(n-n12-n3)), zeros(size(JR34))) ...
                && svds(JR34(1:n3,1:n3), 1, 'smallest') > Opts.rankTol
            % JR34 already has the required form
        else
            % Transform system
            [~,D2,V2] = svd(full(JR34));
            n3 = rank(D2,Opts.rankTol);

            % New system (pH)
            T2 = blkdiag(eye(n12),V2);
            E = T2'*E*T2;
            J = T2'*J*T2;
            R = T2'*R*T2;
            G = T2'*G;
            P = T2'*P;
        end
    end
else
    n3 = 0;
end

%% Decomposition J(4,1:2) = [J41 , 0]
if (n-n12-n3) > 0

    % Get remaining dimensions
    n4 = n-n12-n3;
    n1 = n4;
    n2 = n12-n1;

    Jpart = J(n12+n3+1:end,1:n12);

    if isequal(Jpart-[Jpart(:,1:n1),zeros(n4,n2)], zeros(size(Jpart)))
        % JR34 already has the required form
    else
        % Transform system
        [Uj,~,Vj] = svd(full(Jpart));

        % New system (pH)
        T3 = blkdiag(Vj,eye(n3),Uj);
        E = T3'*E*T3;
        J = T3'*J*T3;
        R = T3'*R*T3;
        G = T3'*G;
        P = T3'*P;

    end
else
    n1 = 0;
    n2 = n12;
    n4 = 0;
end

%% Block-diagonalization of E(1:2,1:2)
if n1 > 0
    if ~isequal(E(n1+1:n12,1:n1), zeros(n2,n1))
        L12 = [eye(n1), -E(n1+1:n12,1:n1)'/E(n1+1:n12,n1+1:n12); zeros(n2,n1) , eye(n2)];
        Le = blkdiag(L12, eye(n3+n4));

        % New system (pH)
        T4 = Le';
        E = T4'*E*T4;
        J = T4'*J*T4;
        R = T4'*R*T4;
        G = T4'*G;
        P = T4'*P;
    end
end

%% Enforce zero patterns
if n12 > 0
    E = blkdiag(E(1:n1,1:n1),E(n1+1:n12,n1+1:n12),zeros(n3+n4));
end

if n1 > 0
    J(n-n4+1:end, n1+1:end) = zeros(n4, n2+n3+n4);
    J(n1+1:end, n-n4+1:end) = zeros(n2+n3+n4, n4);
    R = blkdiag(R(1:n12+n3, 1:n12+n3), zeros(n4));
    P = [P(1:(n12+n3),:); zeros(n4,size(P,2))];
end

%% Final check if system fulfills staircase constraints
try
    sysStair = phs(J,R,eye(size(E)),G,E,P,S,N,Opts.phs);
    dims = [n1,n2,n3,n4];
    staircaseCheck(sysStair,Opts.rankTol);
catch
    error('MORpH:toStaircase:finalSystemNotStaircase', ...
        ['Transformation to staircase form failed. \n' ...
        'Consider changing the tolerance for rank decisions.']);
end

end

function [sys, Opts] = parseInputs(sys, varargin)
% Check phs input type
if ~isa(sys,'phs')
    error('phs:staircaseCheck:wrongInput', 'Model is not an object of the phs-class.');
end

% Opts
if ~isempty(varargin) && isstruct(varargin{end})
    Opts = varargin{end};
else
    Opts = struct();
end

% Option parsing
OptsAdmissible.rankTol = 1e-12;
OptsAdmissible.scaling = struct;
OptsAdmissible.phs = struct;
Opts = phsMOR_parseOpts(Opts,OptsAdmissible);
end
