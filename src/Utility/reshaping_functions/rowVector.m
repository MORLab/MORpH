function vec = rowVector(vec)
% ROWVECTOR - Returns any vector as row vector, throws error for
%             matrices
%
% Syntax:
%   vec = ROWVECTOR(vec)
%
% Input Arguments:
%       *Required Input Arguments:*
%       - vec:	vector, row or column
%
% Output Arguments:
%       - vec:  vector, row
%
%-----------------------------------------------------------------------
% This file is part of
%
% <a href="https://github.com/MORLab/MORpH">MORpH</a> - a MATLAB toolbox to store, analyze,
% interconnect and reduce large-scale port-Hamiltonian models
%
% Authors:      Julius Durmann
% E-Mail:       <a href="mailto:morlab.rt@ed.tum.de">morlab.rt@ed.tum.de</a>
% Website:      <a href="https://www.epc.ed.tum.de/en/rt/home">www.epc.ed.tum.de/rt</a>
% Copyright :   © 2022 Chair of Automatic Control, TUM (see LICENSE.md)
%-----------------------------------------------------------------------

if isrow(vec)
    return
end
if iscolumn(vec)
    vec = vec';
    return
end
error("Input is not a vector!");

end