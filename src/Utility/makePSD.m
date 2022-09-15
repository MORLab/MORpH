function Ap = makePSD(A,varargin)
% MAKEPSD - Makes a symmetric matrix positive semidefinite
%
% Syntax:
%   Ap = MAKEPSD(A)
%   Ap = MAKEPSD(A,1e-20)
%
% Input Arguments:
%       *Required Input Arguments:*
%       - A:	real matrix
%       *Optional Input Arguments:*
%       - delta:	minimum eigenvalue of Ap; can be used to enforce positive
%                   definiteness
%
% Output Arguments:
%       - Ap: 	positive (semi-)definite matrix
%
% See Also:
%       enforcePHStructure
%
% References:
%       [1] N. J. Higham, "Computing a nearest symmetric positive semidefinite matrix",
%           Linear Algebra and its Applications, vol. 103 (1988), pp. 103-118.
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

narginchk(1,2);

if nargin == 1
    % Use eps as default
    delta = eps;
else
    delta = varargin{1};
end

% Enforce symmetry
A = (A+A')/2;

% Eigendecomposition
[Z,ev] = eig(A);

% Projection onto positive (semi-)definite manifold
Ap = Z*diag(max(diag(ev),delta))*Z';