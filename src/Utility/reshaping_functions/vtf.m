function M = vtf(v, n, m)
% VTF - Creates a nxm matrix from a vector with length n*m
%
% Syntax:
%   M = VTF(v, n, m)
%
% Input Arguments:
%       *Required Input Arguments:*
%       - v:	vector with length n*m
%       - n:    number of rows in M
%       - m:    number of columns in M
%
% Output Arguments:
%       - M: 	nxm matrix
%
% See Also:
%       ftv
%
% References:
%       [1] P. Schwerdtner and M. Voigt. Structure Preserving Model Order Reduction 
%           by Parameter Optimization. arXiv Preprint arXiv:2011.07567. 2020. 
%           url: https://arxiv.org/abs/2011.07567.
%
%-----------------------------------------------------------------------
% This file is part of
%
% <a href="https://github.com/MORLab/MORpH">MORpH</a> - a MATLAB toolbox to store, analyze,
% interconnect and reduce large-scale port-Hamiltonian models
%
% Authors:      Maximilian Bonauer
% E-Mail:       <a href="mailto:morlab.rt@ed.tum.de">morlab.rt@ed.tum.de</a>
% Website:      <a href="https://www.epc.ed.tum.de/en/rt/home">www.epc.ed.tum.de/rt</a>
% Copyright :   © 2022 Chair of Automatic Control, TUM (see LICENSE.md)
%-----------------------------------------------------------------------

% Make full nxm matrix from n*m vector
M = reshape(v, [n,m]);

end