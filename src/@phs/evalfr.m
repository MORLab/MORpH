function fresp = evalfr(sys, s)
% EVALFR - Evaluates the transfer function of sys at a single complex
% variable s
%
% Syntax:
%   fresp = EVALFR(sys, s)
%   fresp = sys.EVALFR(s)
%
% Description:
%       evalfr calls sys.transferFunction to obtain the transfer function.
%       It then evaluates the function at s and returns the result.
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:  phs object
%       - s:   (complex) scalar variable
%
% Output Arguments:
%       - fresp: transfer function value at s (sparse matrix for MIMO systems)
%
%
% See Also:
%       phs, phs/transferFunction
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

if ~isscalar(s)
    error('phs:evalfr:frequencyMustBeScalar', ...
        ['evalfr can only compute the response at a single frequency.\n' ...
        'Use freqresp to compute the system response at multiple frequencies.']);
end
G = sys.transferFunction;
fresp = G(s);

end