function sys = makeExplicit(sys, varargin)
% MAKEEXPLICIT - Makes Descriptor system sys explicit by applying a
%       state-space transformation such that E=I
%
% Syntax:
%   sysExp = MAKEEXPLICIT(sys)
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:  phs object
%
% Output Arguments:
%       - sys:   phs object of explicit PH system with E=I
%       *Optional Input Arguments:*
%       - Opts:  structure with execution parameters
%           - keepScaled: Keeps a scaled system with Q=I scaled 
%                         [{false} /true]
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

% Input parsing
narginchk(1,2)

if ~isempty(varargin) && isstruct(varargin{1})
    Opts = varargin{1};
else
    Opts = struct();
end

if sys.isDAE
    error("phs:makeExplicit:E_not_invertible",...
        "Transformation to explicit form is not feasible for pHDAEs.");
end

OptsAdmissible.keepScaled = {false,true};
Opts = phsMOR_parseOpts(Opts,OptsAdmissible);

if Opts.keepScaled && ~isequal(sys.Q, eye(sys.dim))
    error("phs:makeExplicit:Q_not_identity",...
        "Please first use scaling() to bring the system to scaled energy coordinates. ");
end

% Suppress warnings
warning('off', 'phs:phs:changedProperty')

if sys.isImplicit
    if Opts.keepScaled
        % Compute (sparse) Cholesky factor of E
        [Le,flag,P] = chol(sparse(sys.E));
        assert( flag == 0, "Cholesky decomposition of E failed");
        Le = Le*P';
        % Transform system
        sys.J = Le'\sys.J/Le;
        sys.R = Le'\sys.R/Le;
        sys.G = Le'\sys.G;
        sys.P = Le'\sys.P;
        sys.E = speye(sys.dim);

    else
        sys.Q = sys.Q/sys.E;
        sys.E = speye(sys.dim);
    end
end

warning('on', 'phs:phs:changedProperty')

warning('phs:makeExplicit:changedProperty',...
    ['The system was transformed to explicit form.\n' ...
    'Consider running phs.inputValidation(sys) to validate the new representation.']);