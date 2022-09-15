function [redSys, varargout] = pHMOR_DAE_Wrapper(fctHandle, varargin)
% PHMOR_DAE_WRAPPER - Wrapper to apply pH-MOR methods suitable for pHODEs to
% pHDAEs. The wrapper is based on a decomposition of the original transfer
% function to find a pHODE representation of the proper part to which pHODE
% methods may be applied and subsequently, (possibly) improper parts are
% attached again.
%
% Syntax:
%   [redSys,V] = PHMOR_DAE_WRAPPER(@irkaPH,sys,s0,b,Opts)
%   [redSys,hsv]  = PHMOR_DAE_WRAPPER(@ecm,sys,redOrder,Opts)
%
% Description:
%
% The transfer function of the pHDAE system is first decomposed into
%     G(s) = Gp(s) + Pol_1*s
% where Gp(s) is the transfer function of a pHODE system. After applying
% the method with 'fctHandle' to this pHODE part, the improper part is
% attached back to the model.
%
% Input Arguments:
%       *Required Input Arguments:*
%       - fctHandle: function handle of the pH-MOR method
%       *Optional Input Arguments:*
%       - additional inputs to @fctHandle
%       - Opts:
%           - .(func2str(fctHandle)): Opts struct forwarded to MOR method
%           - .decomposePHDAE.*:      Other options that will be passed on to the function;
%                                     Please refer to its documentation (doc decomposePHDAE).
%           - .composePHDAE.*:        Other options that will be passed on to the function;
%                                     Please refer to its documentation (doc composePHDAE).
%
% Output Arguments:
%       - redSys: reduced pHDAE system
%       - all possible output combinations of @fctHandle
%
% See Also:
%       decomposePHDAE, composePHDAE
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

% Parse inputs
[fctHandle, sys, varargin, Opts] = parseInputs(fctHandle, varargin{:});

% Compute pHODE representation for proper part
[sysP, Pol_1] = decomposePHDAE(sys, Opts.decomposePHDAE);

% Reduce proper part
varargout = cell(nargout, 1);
[varargout{:}] = fctHandle(sysP, varargin{:},Opts.(func2str(fctHandle)));
redSysP = varargout{1};
varargout(1) = [];

% Attach improper part
redSys = composePHDAE(redSysP, Pol_1, Opts.composePHDAE);

end

function [fctHandle, sys, varargin, Opts] = parseInputs(fctHandle, varargin)
% Check phs input type
if ~isa(fctHandle,'function_handle')
    error('MORpH:pHMOR_DAE_Wrapper:wrongInput', 'First argument is not a function_handle.');
end

% phs
sys = varargin{1};
varargin(1) = [];

% Opts
if ~isempty(varargin) && isstruct(varargin{end})
    Opts = varargin{end};
    varargin(end) = [];
else
    Opts = struct();
end
% Option parsing
OptsAdmissible.(func2str(fctHandle)) = struct;
OptsAdmissible.decomposePHDAE = struct;
OptsAdmissible.composePHDAE = struct;
Opts = phsMOR_parseOpts(Opts,OptsAdmissible);

end