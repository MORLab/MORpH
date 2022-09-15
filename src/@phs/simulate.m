function [tout, Y, X] = simulate(sys, u, times, x0, Opts)
% SIMULATE - Simulates the phs system sys according to a given input over time
%
% Syntax:
%   [tout, Y, X] = sys.simulate(u, times, x0)
%   [tout, Y, X] = sys.simulate(u, times, x0; Opts)
%   [tout, Y, X] = simulate(sys, u, times, x0)
%   [tout, Y, X] = simulate(sys, u, times, x0, Opts)
%
% Description:
%       This function solves the differential equation of the system
%           x_dot(t) = (J - R)*Q*x(t) + (B - P)*u(t)
%       by using a differential equation solver such as ode45.
%
%       Since the solver determines the required time steps itself, the
%       values of u will be linearly interpolated if necessary.
%
%       You can use a different solver by specifying the option Opts.solver
%       to be any function handle with the same signature as ode45.
%
%       If you wish to obtain the values at exactly the specified time
%       steps, set Opts.matchTime to true. This will then interpolate the
%       results of the ode solver with spline.
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys: phs object
%       - u: Input data as vector. Dimensions must be: (N, p) where N is
%               the number of time steps and p is the input dimension of
%               the system
%       - times: Vector of time steps
%       - x0: The initial state of the system (column vector)
%       *Optional Input Arguments:*
%       - Opts:  structure with execution parameters
%           - .solver: 	    Function handle to any ode solver which has the
%                           same signature as ode45
%                           [{@ode45} / function handle]
%           - .matchTime: 	Set to true if Y and X should be interpolated
%                           to match the time steps from times
%                           [{false} / true]
%
% Output Arguments:
%       - tout: 	Vector of time steps corresponding to Y and X,
%                   dimensions: (Nout)
%       - Y: 		System output at time steps tout,
%                   dimensions: (Nout, p)
%       - X:        System states at time steps tout,
%                   dimensions: (Nout, n)
%
% Examples:
%       The code below creates a simple mass-spring-damper system and
%       simulates it for a sinusoidal input.
%
%       sys = setup_MassSpringDamperSystem(n, 2, 1, 1, 'SISO');
%       Opts.solver = @ode45;
%       Opts.matchTime = false;
%
%       t = linspace(0, 100, 1000)';
%       u = sin(t);
%       x0 = [1; zeros(n-1, 1)];
%
%       [tout, Y] = sys.simulate(u, t, x0, Opts);
%       Y_lsim = sys.lsim(u, t, x0);
%
%       figure()
%       title("Sine response")
%       plot(tout, Y);
%
% See Also:
%       phs, ode45, spline
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

%% Input parsing
narginchk(4, 5);
if nargin < 5
    Opts = struct();
end
OptsAdmissible.solver = @ode45;
OptsAdmissible.matchTime = false;
Opts = phsMOR_parseOpts(Opts, OptsAdmissible);
solver = Opts.solver;
if size(x0, 2) ~= 1
    x0 = x0';
end

%% Solve differential equations
if sys.isImplicit
    fprintf("Trying to transform system to explicit representation...\n")
    if sys.isDAE
        error("phs:simulate:DAEnotSupported", "'simulate' does not support DAE systems");
    end
    sys = sys.makeExplicit;
    fprintf("Internal system representation is now explicit.\n")
end

f = @(t, x) (sys.J - sys.R)*sys.Q*x + (sys.G-sys.P)*interpolate_u(t);
[tout, X] = solver(f, times, x0);
Y = ((sys.G + sys.P)'*sys.Q*X')';

%% Match time if necessary
if Opts.matchTime
    Y = spline(tout, Y', times)';
    X = spline(tout, X', times)';
    tout = times;
end

%% Auxiliary
    function u_i = interpolate_u(t)
        % Find index where time series first exceeds t
        idx = find(times >= t, 1);
        if times(idx) == t
            u_i = u(idx, :);
        else
            % Linear interpolation
            u_i = (t - times(idx-1))/(times(idx) - times(idx-1)) * (u(idx, :) - u(idx-1, :)) + u(idx-1, :);
        end
        u_i = u_i';
    end

end

