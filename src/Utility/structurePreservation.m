function [J_red, R_red, Q_red, G_red, E_red, P_red, S_red, N_red, V, W] = structurePreservation(sys, V, Opts)
% STRUCTUREPRESERVATION - Returns reduced system matrices of sys with V by
%       application of PH structure preserving compuations
%
% Syntax:
%   [J_red, R_red, Q_red, G_red, E_red, P_red, S_red, N_red, V, W]
%                                   = STRUCTUREPRESERVATION(sys, V, Opts)
%
% Description:
%       This function uses the Krylov-subspace matrix V (possibly obtained
%       by the arnoldi algorithm) to compute reduced system matrices J_red,
%       R_red, etc. of a reduced port-Hamiltonian system. The structure
%       preservation mechanism can be set by Opts.structurePreservation
%       (see Input Arguments below).
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:  phs object which is to be reduced
%       - V:    Krylov subspace matrix, should be orthogonal
%       - Opts: structure with execution parameters
%           - .structurePreservation: determines how structure is
%                   preserved by the algorithm:
%                   These options are available:
%                   > 'scaling': use scaled energy coordinates and W = V; [1]
%                       Note: Algorithm expects system where Q = eye(N);
%                   > 'specialInverse': W = Q*V*inv(V'*Q*V); [3]
%                   > 'Cholesky': V'*Q*V = L'*L, V <- V*inv(L), W <- QV; [2]
%                   > 'Cholesky+': V'*Q'*E*V = L'*L, V <- V*inv(L), W <- QV
%                   > 'QV':       V <- V, W <- QV [4]
%                   > false: Does not perform structure preservation, only
%                            V is returned
%                   [{'Cholesky'} / 'scaling' / 'specialInverse' / 'QV' / ...
%                    'Cholesky+']
%           - .reorthog: {'orth', true, false}, default: 'mgs'
%                   Determines if the Krylov-subspace matrix V will be
%                   orthogonalized before further computations
%                   > 'orth' uses MATLAB function orth
%                   > 'qr' uses MATLAB function qr
%                   > false does not use orthogonalization
%                   > true uses default ('qr')
%                   [{true} / 'qr' / 'orth' / false]
% Output Arguments:
%       - J_red, R_red, Q_red, E_red, P_red, S_red, N_red:
%                   Reduced system matrices
%       - V, W:     Projection matrices (V may differ from input!)
%
% Examples:
%       See usage in arnoldiPH.
%
% See Also:
%       arnoldiPH
%
% References:
%       [1] Polyuga (2010), "Model Reduction of Port-Hamiltonian Systems",
%           Dissertation, 2010.
%       [2] S. Gugercin et al. “Structure-Preserving Tangential Interpolation 
%           for Model Reduction of Port-Hamiltonian Systems." 
%           In: Automatica 48.9 (2012), pp. 1963–1974.
%       [3] T. Wolf, B. Lohmann, R. Eid, and P. Kotyczka. Passivity and 
%           structure preserving order reduction of linear port-Hamiltonian 
%           systems using Krylov subspaces. Eur. J. Control, 16(4):401–406, 2010.
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

OptsAdmissible = struct();
OptsAdmissible.reorthog = {true,'orth','qr',false};
OptsAdmissible.structurePreservation = {'Cholesky','scaling','specialInverse','Cholesky+','QV',false};
Opts = phsMOR_parseOpts(Opts, OptsAdmissible);

% Orthogonalization of V
old_dim = size(V,2);
if isfield(Opts,'reorthog')
    switch Opts.reorthog
        case false

        case true
            [V,~] = qr(V,0);
        case 'orth'
            V = orth(V);
        case 'qr'
            [V,~] = qr(V,0);
        otherwise
            error('MORpH:structurePreservation:unknownOption',...
                ['Input ''', Opts.reorthog, ''' is not supported for option reorthog.'])
    end
end
if size(V,2) < old_dim
    warning('MORpH:structurePreservation:dimensionsLost',...
        ['Some dimensions were lost during reorthogonalization!\n'...
        'This may happen because the reduced order exceeds the minimal'...
        ' order of the system.\n'...
        'Dropped from ', num2str(old_dim), ' to ', num2str(size(V,2)),...
        ' dimensions'])
end

% Structure preservation strategies
sV = size(V);
J = sys.J; R = sys.R; Q = sys.Q; G = sys.G; E = sys.E; P = sys.P; S = sys.S; N = sys.N;
switch Opts.structurePreservation
    case 'scaling'
        assert(isequal(Q,eye(sys.dim)))
        W = V;

        J_red = W'*J*V;
        R_red = W'*R*V;
        Q_red = eye(sV(2), sV(2));  % identity matrix due to scaling
        if sys.isImplicit     % avoid having small numerical inaccuracies leading to descriptor systems!
            E_red = W'*E*V;
        else
            E_red = eye(length(J_red));
        end
        G_red = W'*G;
        P_red = W'*P;
        S_red = S;
        N_red = N;
    case 'specialInverse'
        W = Q*V/(V'*Q*V);

        J_red = W'*J*W;
        R_red = W'*R*W;
        Q_red = V'*Q*V;

        if sys.isImplicit     % avoid having small numerical inaccuracies leading to descriptor systems!
            E_red = W'*E*V;
        else
            E_red = eye(length(J_red)); % (W'*V = I by definition of W)
        end
        G_red = W'*G;
        P_red = W'*P;
        S_red = sys.S;
        N_red = sys.N;
    case 'Cholesky'
        L = chol(V'*Q*V);
        V = V/L;
        W = Q*V;

        J_red = W'*J*Q*V;
        R_red = W'*R*Q*V;
        Q_red = eye(size(J_red));
        if sys.isImplicit     % avoid having small numerical inaccuracies leading to descriptor systems!
            E_red = W'*E*V;
        else
            E_red = eye(length(J_red)); % (W'*V = I by definition of W)
        end
        G_red = W'*G;
        P_red = W'*P;
        S_red = sys.S;
        N_red = sys.N;
    case 'Cholesky+'
        L = chol(V'*Q'*E*V);
        V = V/L;
        W = Q*V;

        J_red = W'*J*Q*V;
        R_red = W'*R*Q*V;
        E_red = eye(length(J_red)); % (W'*E*V = I by definition of W and V)
        Q_red = eye(length(J_red));
        G_red = W'*G;
        P_red = W'*P;
        S_red = sys.S;
        N_red = sys.N;
    case 'QV'
        W = Q*V;

        J_red = W'*J*Q*V;
        R_red = W'*R*Q*V;
        E_red = W'*E*V;
        Q_red = eye(length(J_red));
        G_red = W'*G;
        P_red = W'*P;
        S_red = sys.S;
        N_red = sys.N;
    case false
        % Used for modelFctPH
        W = [];
        J_red = [];
        R_red = [];
        E_red = [];
        Q_red = [];
        G_red = [];
        P_red = [];
        S_red = [];
        N_red = [];
    otherwise
        error('Structure Preservation method unknown')
end

end
