function G = transferFunction(sys)
% TRANSFERFUNCTION - Returns a transfer function handle for the phs object
%                   sys
%
% Syntax:
%   G = sys.TRANSFERFUNCTION
%   G = TRANSFERFUNCTION(sys)
%
% Description:
%       The transfer function G(s) of the port-Hamiltonian system sys is
%       computed as follows:
%       G = @(s) (sys.G + sys.P)'*sys.Q
%           * (s*sys.E - (sys.J - sys.R)*sys.Q)\(sys.G - sys.P)
%           + (sys.S + sys.N);
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys: phs object
%
% Output Arguments:
%       - G: function handle
%
% See Also:
%       phs, phs/evalfr
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

% Transfer function
G = @(s) (sys.G + sys.P)'*sys.Q * ((s*sys.E - (sys.J - sys.R)*sys.Q)\(sys.G - sys.P)) + (sys.S + sys.N);

end