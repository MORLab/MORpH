function shifts = initShifts(sys, strategy, number, Opts)
% INITSHIFTS - Create initial shifts for Krylov subpace reduction
%
% Syntax:
%   shifts = INITSHIFTS(sys, strategy, number)
%   shifts = INITSHIFTS(sys, strategy, number, Opts)
%
% Description:
%       initShifts creates a column vector of <number> initial shifts for
%       the PH system <sys>. The provided string <strategy> lets you choose
%       between different strategies, including:
%
%       'zeros':    A vector of zeros
%       'linear':   Linearly spaced, real-valued shifts from
%                   10^(Opts.range(1)) to 10^Opts.range(2)
%       'logarithmic':  Logarithmically spaced, real-valued shifts from
%                   10^(Opts.range(1)) to 10^Opts.range(2)
%       'eig_large':    Shifts based on the largest eigenvalues (absolute
%                   value) of the system
%       'eig_small':    Shifts based on the smallest eigenvalues (absolute
%                   value) of the system
%       'diag':     Shifts based on the small and large diagonal entries of
%                   the system's dynamic matrix (J-R)*Q
%       'eig_circle':   Shifts sampled randomly in the right half of the
%                   complex plane between two circles with radius
%                   |smallest eig| and |largest eig|
%
%       You can pass an optional structure Opts to set the range values for
%       the above methods (these values will also be used as fallback
%       values if eigs returns nan).
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:  phs object
%       - strategy: ['zeros' / 'linear' / 'logarithmic' / 'eig_large' /...
%                    'eig_small' / 'diag' / 'eig_circle']
%       - number:   Number of shifts to be created
%
%       *Optional Input Arguments:*
%       - Opts:  structure with execution parameters
%           - range:    vector [min, max] of range/fallback values
%                       [{[-10,10]} / 1-by-2 vector]
%
% Output Arguments:
%       - shifts:   Column vector of <number> shifts
%
% Examples:
%   This example creates a simple mass-spring-damper system of order 100
%   and creates 20 shifts.
%
%       sys = setup_MassSpringDamperSystem(100, 2, 1, 1)
%       shifts = initShifts(sys, 'eig_circle', 20)
%       scatter(real(shifts), imag(shifts))
%       grid on
%
% See Also:
%       arnoldiPH, irkaPH
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

if nargin < 4
    Opts = struct();
end

OptsAdmissible.range = [-10, 10];
Opts = phsMOR_parseOpts(Opts, OptsAdmissible);

switch strategy
    case 'zeros'
        shifts = zeros(number, 1);
    case 'linear'
        shifts = linspace(10^(Opts.range(1)), 10^Opts.range(2), number);
    case 'logarithmic'
        shifts = logspace(Opts.range(1), Opts.range(2), number)';
    case 'eig_large'
        shifts = - sys.eigs(number, 'largestabs');
        shifts(isnan(shifts)) = 0;
    case 'eig_small'
        shifts = - sys.eigs(number, 'smallestabs');
        shifts(isnan(shifts)) = 0;
    case 'diag'
        % Take smallest and largest diagonal entries of dynamic
        % matrices. This choice is motivated by the Gerschgorin
        % circles.
        diags = sort(diag(full((sys.J - sys.R)*sys.Q)));
        shifts = - [diags(1:ceil(number/2)); diags(end-ceil(number/2):end)];
        if length(shifts) > number
            shifts = shifts(1:number);
        end
    case 'eig_circle'
        % Take smallest and largest eigenvalue (absolute values) to
        % draw two half circles in the right complex half plane.
        % Randomly distribute values between these half circles.
        outer_radius = abs(sys.eigs(1, 'largestabs'));
        inner_radius = abs(sys.eigs(1, 'smallestabs'));
        if isnan(outer_radius) || isnan(inner_radius)
            inner_radius = 0;
            outer_radius = Opts.range(2);
            warning("MORpH:initShifts:eigs_not_converged", "Eigs did not converge! Falling back to default values...")
        end
        mag = inner_radius + (outer_radius - inner_radius)*rand(floor(number/2), 1);
        phase = -pi/2 + pi*rand(floor(number/2), 1);   % Angles in the right half plane
        shifts = [mag .* cos(phase) + 1i*mag .* sin(phase); ...
            mag .* cos(phase) - 1i*mag .* sin(phase)];
        if mod(number, 2) ~= 0
            shifts = [shifts; inner_radius + (outer_radius - inner_radius)*rand(1)];
        end
    otherwise
        msg = strcat("The provided strategy " + strategy + " for shift initialization is not implemented");
        error("MORpH:initShifts:strategyNotImplemented", msg);
end

end