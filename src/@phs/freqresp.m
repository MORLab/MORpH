function [H, W] = freqresp(sys, varargin)
% FREQRESP - Computes the frequency response H of the PH system sys at
%           frequencies W (will be determined based on the system dynamics
%           if not provided).
%
% Syntax:
%   [H, W] = FREQRESP(sys)
%   H = FREQRESP(sys,W)
%   H = FREQRESP(sys,W,units)
%
% Description:
%       See description of freqresp.
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:	phs-object
%       *Optional Input Arguments:*
%       - W:        Frequency grid (radians/second)
%       - units:    Specify frequency units of W
%                   ['rad/TimeUnit', 'cycles/TimeUnit', 'rad/s', 'Hz',
%                    'kHz', 'MHz', 'GHz', 'rpm']
%
% Output Arguments:
%       - H: 	Frequency response values belonging to W; complex
%       - W: 	Frequency grid (3-Dimensional array: p x p x length(w))
%                   where p is the number of in-/outputs of sys
%
% See Also:
%       DynamicSystem/freqresp, phs
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

% call DynamicSystem/freqresp
[H, W] = freqresp(ss(sys),varargin{:});
