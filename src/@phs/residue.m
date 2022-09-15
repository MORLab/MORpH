function [r, p, d] = residue(sys)
% RESIDUE - Computes the residues r belonging to the poles p and the
%           feedthrough term d of the system sys
%
% Syntax:
%   [r, p, d] = RESIDUE(sys)
%   [r, p, d] = sys.RESIDUE
%
% Description:
%       RESIDUE determines the poles and residues, i.e. the asymptotic
%       points and the coefficients of the partial fraction expansion of
%       the system's transfer function:
%                r1       r2             rn
%       G(s) = ------ + ------ + ... + ------ + d
%              s - p1   s - p2         s - pn
%
%       RESIDUE transforms the system sys to diagonal form using eig. It
%       then uses the transformed input and output matrices to determine
%       the residues.
%       Note that higher-order poles pi which would normally contribute to
%       polynomial functions in the denominator are split into multiple
%       occurences of (s-pi) due to how eig handles multiple eigenvalues.
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys
%
% Output Arguments:
%       - r: 	Cell vector of residues (scalars or matrices)
%       - p:    Vector of poles
%       - d:    Feedthrough term (scalar or matrix)
%
% Examples:
%       The code below creates a simple mass-spring-damper system and
%       computes the poles, residue, and feedthrough term.
%
%       sys = setup_MassSpringDamperSystem(10,2,1,1)
%       [r, p, d] = sys.residue
%
% See Also:
%       phs, eig
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

% Bring system to diagonal form
[V, D] = sys.eig;

if rank(V, 1e-7) < length(V)
    warning('phs:residue:MatrixNotDiagonalizable',...
        ['The dynamic matrix (J-R)*Q of the system is not diagonalizable. ' ...
        'The pole-residue decomposition of non-diagonalizable systems is currently not supported.']);
end

B = (sys.E*V)\(sys.G - sys.P);
C = (sys.G + sys.P)'*sys.Q*V;

% Compute residues from projected B and C matrices
r = cell(sys.dim, 1);
for k = 1:sys.dim
    r{k} = C(:,k)*B(k,:);
end

% Poles
p = diag(D);

% Feedthrough
d = sys.N + sys.S;

end