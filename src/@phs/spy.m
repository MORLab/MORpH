function figureHandle = spy(sys,figureHandle)
% SPY - Creates a figure with the spy-plots of all system matrices
%
% Syntax:
%   sys.SPY
%   figureHandle = sys.SPY
%   figureHandle = SPY(sys)
%   sys.SPY(figureHandle)
%   SPY(sys,figureHandle)
%
% Description:
%       spy calls spy for all system matrices and collects them in a single
%       figure. You may provide the figure handle to be plotted in by
%       yourself.
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:  phs-object
%       *Optional Input Arguments:*
%       - figureHandle: Figure handle to be plotted in
%
% Output Arguments:
%       - figureHandle: handle to the created Figure
%
% See Also:
%       spy, phs
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

if nargin == 1
    figureHandle = figure();
elseif nargin == 2
    set(0,'currentfigure', figureHandle)
end

matOrder = {'E','J','G','N','Q','R','P','S'};
for i = 1:8
    subplot(2,4,i)
    spy(sys.(matOrder{i}));
    title(matOrder{i})
end

end