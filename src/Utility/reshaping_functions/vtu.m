function M = vtu(v)
% VTU - Creates an upper triangular nxn matrix from a vector with
% length n(n+1)/2
%
% Syntax:
%   M = VTU(v)
%
% Input Arguments:
%       *Required Input Arguments:*
%       - v:	vector with length n(n+1)/2
%
% Output Arguments:
%       - M: 	nxn matrix
%
% See Also:
%       utv
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

% Make Upper triangular matrix from vector
n = round((sqrt(8 * length(v) + 1) - 1) / 2);
M = zeros(n);
M(tril(ones(n),0)==1) = v;
M = M.';

end