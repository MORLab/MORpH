function [varargout] = eig(sys, varargin)
% EIG - computes eigenvalues and eigenvectors of the
% generalized eigenvalue problem (J-R)*Q*V = E*V*D for a PH system sys
%
% Syntax:
%   lambda = sys.EIG
%   lambda = EIG(sys)
%   lambda = sys.EIG(varargin)
%   lambda = EIG(sys, varargin)
%   [V, D] = sys.EIG
%   [V, D] = EIG(sys, varargin)
%   [V, D, W] = sys.EIG
%   [V, D, W] = EIG(sys, varargin)
%
% Description:
%       The function returns eig(full((sys.J-sys.R)*sys.Q)), i.e. the
%       eigenvalues of the systems dynamic matrix.
%       For descriptor systems, eig(full((sys.J-sys.R)*sys.Q), full(sys.E))
%       is called, i.e. the solutions lambda of the generalized eigenvalue
%       problem
%           (J-R)*Q * v = lambda * E * v
%       are computed.
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:  phs-object
%       *Optional Input Arguments:*
%       - see optional input arguments of eig
%
% Output Arguments:
%       - lambda:   Eigenvalues / Poles
%       - V, W:     Matrices of right and left eigenvectors
%       - D:        Matrix of (generalized) eigenvalues
%
% Examples:
%       This example creates a simple PH system and calls eig.
%
%       J = [0 -1; 1 0]; R = [0 0; 0 1]; Q = eye(2); G = [1; 0];
%       sys = phs(J, R, Q, G);
%       lambda = sys.eig
%
% See Also:
%       eig, phs
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

varargout = cell(nargout, 1);
if ~sys.isImplicit
    [varargout{:}] = eig(full((sys.J-sys.R)*sys.Q), varargin{:});
else
    [varargout{:}] = eig(full((sys.J-sys.R)*sys.Q), full(sys.E), varargin{:});
end

end