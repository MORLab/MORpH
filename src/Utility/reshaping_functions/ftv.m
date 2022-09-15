function v = ftv(M)
% FTV - Creates a vector with length n*m from an nxm matrix
%
% Syntax:
%   v = FTV(M)
%
% Input Arguments:
%       *Required Input Arguments:*
%       - M:	nxm matrix
% 
% Output Arguments:
%       - v: 	vector with length n*m
%
% See Also:
%       vtf
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

% Make n*m vector from full nxm matrix
    v = reshape(M, [numel(M), 1]);

end