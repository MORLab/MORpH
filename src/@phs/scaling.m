function sys = scaling(sys, varargin)
% SCALING - Applies a state-space transformation to a port-Hamiltonian
%   system such that its energy matrix Q equals the identity matrix.
%
% Syntax:
%   sys_scaled = SCALING(sys)
%
% Description:
%       Find new state coordinates z, x = S*z, such that S'*Q'S = I.
%       See the references for more details.
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:  phs object
%       *Optional Input Arguments:*
%           - Opts:
%               - .rankTol: tolerance for rank decision and splitting of Q
%                           [{1e-16} / positive double]                    
%
% Output Arguments:
%       - sys: phs object where Q = I
%
% Examples:
%       sys = setup_MassSpringDamperSystem(20,2,1,1)
%       isequal(sys.Q, eye(sys.dim))
%       sys_scaled = scaling(sys)
%       isequal(sys_scaled.Q, eye(sys.dim))
%       bode(sys,sys_scaled)
%
% See Also:
%       phs, makeExplicit
%
% References:
%       [1] C. Mehl, V. Mehrmann, and M. Wojtylak. Distance problems for dissipative Hamiltonian
%           systems and related matrix polynomials. Linear Algebra Appl., pages 335–366, 2021.
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

%% Input parsing
narginchk(1,2)

if ~isempty(varargin) && isstruct(varargin{1})
    Opts = varargin{1};
else
    Opts = struct();
end

OptsAdmissible.rankTol = 1e-16;
Opts = phsMOR_parseOpts(Opts,OptsAdmissible);

% Suppress warnings
warning('off', 'phs:phs:changedProperty')

%% Transform system
if ~isequal(sys.Q, eye(sys.dim))
    if svds(sys.Q, 1, 'smallest') < Opts.rankTol

        % Compute SVD of Q
        [U,D,V] = svd(full(sys.Q));
        rankQ = rank(D,Opts.rankTol);

        % Transform system
        sys.E = U'*sys.E*V;
        sys.J = U'*sys.J*U;
        sys.R = U'*sys.R*U;
        sys.Q = D;
        sys.G = U'*sys.G;
        sys.P = U'*sys.P;

        % Extract upper part
        sys.E = sys.E(1:rankQ,1:rankQ);
        sys.J = sys.J(1:rankQ,1:rankQ);
        sys.R = sys.R(1:rankQ,1:rankQ);
        sys.Q = sys.Q(1:rankQ,1:rankQ);
        sys.G = sys.G(1:rankQ,:);
        sys.P = sys.P(1:rankQ,:);

    end

    % Q is invertible
    sys.E = sys.Q'*sys.E;
    sys.J = sys.Q'*sys.J*sys.Q;
    sys.R = sys.Q'*sys.R*sys.Q;
    sys.G = sys.Q'*sys.G;
    sys.P = sys.Q'*sys.P;
    sys.Q = speye(size(sys.E));
end

warning('on', 'phs:phs:changedProperty')

warning('phs:scaling:changedProperty',...
    ['The system was transformed to scaled energy coordinates.\n' ...
    'Consider running phs.inputValidation(sys) to validate the new representation.']);