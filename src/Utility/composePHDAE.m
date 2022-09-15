function [sysC] = composePHDAE(sys, Pol_1, varargin)
% COMPOSEPHDAE - Adds an improper polynomial part to the transfer function
% of a pHDAE sys
%
% Syntax:
%   sysC = COMPOSEPHDAE(sys, Pol_1)
%   sysC = COMPOSEPHDAE(sys, Pol_1, Opts)
%
% Description:
%
% Composition of a transfer function
%     Gc(s) = G(s) + Pol_1*s
%
% where G(s) is the transfer function of a pHDAE system sys, into a pHDAE
% system
%     E*dx/dt = (J-R)*Q*x(t)  + (G-P)*u(t),
%           y = (G+P)'*x(t) + (S+N)*u(t),
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:      phs object
%       - Pol_1:    square matrix coefficient of linear polynomial part; positive
%                   semidefinite and with same dimensions as sys.S and sys.N
%       *Optional Input Arguments:*
%       - Opts:  structure with execution parameters
%           - .phs.*: Other options that will be passed on to the class 'phs';
%                     Please refer to its documentation (doc phs).
%
% Output Arguments:
%       - sysC:  pHDAE with integrated improper part
%
% See Also:
%       decomposePHDAE
%
% References:
%       [1] T. Moser et al. Structure-Preserving Model Order Reduction for Index 
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
narginchk(2,3);
[sys, Pol_1, Opts] = parseInputs(sys, Pol_1, varargin{:});

%% Compose system
m = size(Pol_1,1);
if ~isequal(Pol_1,zeros(m,m))
    % Factorization of Pol_1
    [Us,Ds] = eig(full(Pol_1));
    Ld = Us*sqrt(Ds);

    % System composition
    Ec = blkdiag(sys.E,eye(m),zeros(m));
    Jc = blkdiag(sys.J,[zeros(m),-eye(m);eye(m),zeros(m)]);
    Rc = blkdiag(sys.R,zeros(2*m));
    Qc = blkdiag(sys.Q,eye(2*m));
    Gc = [sys.G;zeros(m);Ld'];
    Pc = [sys.P;zeros(2*m,m)];
    Sc = sys.S;
    Nc = sys.N;

    if isa(sys,'phsRed')
        sysC = phsRed(Jc,Rc,Qc,Gc,Ec,Pc,Sc,Nc,Opts.phs);
        sysC.method = sys.method;
        sysC.info = sys.info;
        sysC.parameters = sys.parameters;
    else
        sysC = phs(Jc,Rc,Qc,Gc,Ec,Pc,Sc,Nc,Opts.phs);
    end
else
    sysC = sys;
end

end

function [sys, Pol_1, Opts] = parseInputs(sys, Pol_1, varargin)
% Check phs input type
if ~isa(sys,'phs')
    error('MORpH:composePHDAE:wrongInput', 'Model is not an object of the phs-class.');
end
if ~isnumeric(Pol_1) || ~isequal(size(Pol_1),size(sys.S))
    error('MORpH:composePHDAE:wrongInput', 'Pol_1 has incorrect dimensions or datatype.');
elseif min(eig(full(Pol_1))) < 0
    error('MORpH:composePHDAE:wrongInput', 'Pol_1 is not positive semidefinite.');
end
% Opts
if ~isempty(varargin) && isstruct(varargin{end})
    Opts = varargin{end};
else
    Opts = struct();
end

% Option parsing
OptsAdmissible.phs = struct;
Opts = phsMOR_parseOpts(Opts,OptsAdmissible);

end