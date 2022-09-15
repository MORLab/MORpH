function bodemag(varargin)
% BODEMAG - plots only the magnitude of the bode plot of the
%           port-Hamiltonian system sys
%
% Syntax:
%   BODEMAG(sys)
%   ... (see DynamicSystem/bodemag)
%   [...] = BODEMAG(..., Opts)
%
% Description:
%       phs/bodemag is a wrapper for DynamicSystem/bodemag!
%
%       bodemags(sys) plots the magnitude part of the bode diagram of the
%       PH-system sys
%
%       This function uses phs/freqresp to transform any phs-object into
%       frd-objects by default. For simple systems, this can speed up
%       computation. See also Opts.useFrd in the Input Arguments section.
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:       phs-object
%       *Optional Input Arguments:*
%       - Opts:     Structure with execution parameters
%           - .useFrd:  If true, all phs systems from the input will be
%                       transformed to frd-objects with frequencies
%                       determined by phs/freqresp. If false, all phs
%                       systems will be transformed to ss-objects.
%                       [{true} / false]
%
% Examples:
%       The following code creates a simple phs object and runs the
%       bodemag-function:
%
%       J = [0 -1; 1 0]; R = [0 0; 0 1]; Q = eye(2); G = [1; 0];
%       sys = phs(J, R, Q, G);
%       bodemag(sys);
%
% See Also:
%       phs, DynamicSystem/bodemag
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

% Call DynamicSystem/bodemag
bodemag(varargin{:});

end