function [varargout] = bode(varargin)
% BODE - This function plots the bode-plot of the port-Hamiltonian system sys
%
% Syntax:
%   BODE(sys)
%   BODE(sys, omega)
%   BODE(sys1, sys2, ..., sysN)
%   BODE(sys1, sys2, ..., sysN, omega)
%   BODE(sys1,'-r',sys2,'--k');
%   BODE(sys1,'-r',sys2,'--k',omega);
%   [mag, phase, omega] = BODE(sys)
%   [mag, phase, omega] = BODE(sys, omega)
%   ... (see DynamicSystem/bode)
%   [...] = BODE(..., Opts)
%
% Description:
%       phs/bode is a wrapper for DynamicSystem/bode!
%
%       bode(sys) plots the bode diagram of the port-Hamiltonian system
%       sys. It therefore transforms the phs-object into a ss-object and
%       then runs the function DynamicSystem/bode.
%
%       This function uses phs/freqresp to transform any phs-object into
%       frd-objects by default. For simple systems, this can speed up
%       computation. See also Opts.useFrd in the Input Arguments section.
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:          phs-object, system of interest
%       *Optional Input Arguments:*
%       - omega:        Can be a cell array with min and max value for
%                       frequency OR
%                       a vector of frequencies (e.g. created using
%                       logspace function)
%       - sys2, ...:    Other phs-objects whose bode plot should be plottet
%                       simultaneously (only works if no output is
%                       requested)
%       - <color specifiers>: see Syntax, provides line style options for
%                           plotting
%       - Opts:         Structure with execution parameters
%           - .useFrd:  If true, all phs systems from the input will be
%                       transformed to frd-objects with frequencies
%                       determined by phs/freqresp. If false, all phs
%                       systems will be transformed to ss-objects.
%                       [{true} / false]
%
% Output Arguments:
%       *Optional output arguments:*
%       - mag:          magnitude values, multidimensional array of size
%                       (number of system inputs, number of system outputs,
%                       length of frequency vector omega)
%       - phase:        phase values (in degrees), multidimensional array
%                       of size (number of system inputs, number of system
%                       outputs, length of frequency vector omega)
%       - omega:        vector of frequencies
%       - *:            See DynamicSystem/bode
%
%
% Examples:
%       The following code creates a simple phs object and runs the
%       bode-function:
%
%       J = [0 -1; 1 0]; R = [0 0; 0 1]; Q = eye(2); G = [1; 0];
%       sys = phs(J, R, Q, G);
%       bode(sys);
%
% See Also:
%   phs, DynamicSystem/bode
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

% Parse inputs
Opts = struct();
if isa(varargin{end}, 'struct')
    Opts = varargin{end};
    varargin = varargin(1:end-1);
end
OptsAdmissible.useFrd = true;
Opts = phsMOR_parseOpts(Opts, OptsAdmissible);

omega = [];
if isa(varargin{end}, 'double') || isa(varargin{end}, 'cell')
    omega = varargin{end};
end

% Transform phs (phsRed) objects to frd- or ss-objects
for i = 1:length(varargin)
    if isa(varargin{i},'phs') || isa(varargin{i},'phsRed')
        if Opts.useFrd
            if isempty(omega)
                [H, W] = varargin{i}.freqresp();
            else
                [H, W] = varargin{i}.freqresp(omega);
            end
            varargin{i} = frd(H, W);
        else
            varargin{i} = ss(varargin{i});
        end
    end
end

% Call DynamicSystem/bode
varargout = cell(1,nargout);
[varargout{:}] = bode(varargin{:});

end