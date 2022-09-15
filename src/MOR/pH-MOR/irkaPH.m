function [redSys, V, s0, b, W, nLU] = irkaPH(sys, redOrder, varargin)
% IRKAPH - calculates the reduced order system redSys by application of the
%          Iterative Rational Krylov Algorithm (IRKA algorithm) for PH
%          systems [1]
%
% Syntax:
%   redSys = IRKAPH(sys, redOrder)
%   redSys = IRKAPH(sys, redOrder, Opts)
%   redSys = IRKAPH(sys, s0)
%   redSys = IRKAPH(sys, s0, Opts)
%   redSys = IRKAPH(sys, s0, b)
%   redSys = IRKAPH(sys, s0, b, Opts)
%   [redSys, V, s0, b, W, nLU] = IRKAPH(sys, redOrder, ...)
%
% Description:
%       redSys = irkaPH(sys, redOrder) returns a reduced PH system of order
%       redOrder.
%
%       redSys = irkaPH(sys, s0) returns a reduced PH system
%       where s0 is used as initial guess for the interpolation
%       points.
%
%       redSys = irkaPH(sys, s0, b) returns a reduced PH system
%       where s0 and b are used as initial guess for the interpolation
%       points and tangent directions.
%
%       Argument Opts can be used for further adjustments (see below).
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:		phs object
%       - redOrder: desired order of the reduced system OR
%                   vector of initial shifts
%       *Optional Input Arguments:*
%       - Opts:     Structure with execution parameters
%           - .tol:         Tolerance for the exit condition. If the
%                           stopping criterion defined by Opts.stopCrit
%                           falls below this value, iteration will stop.
%                           [{1e-3} / positive double]
%           - .maxIter: 	Defines the maximum amount of iterations
%                           [{500} / positive integer]
%			-.stopCrit:     Defines the stopping criterion
%                           Available options:
%                           > 's0': (Relative) shift difference between
%                               two iterations
%                           > 'sysr': (Relative) reduced system difference
%                               (norm) between two iterations
%                           > 's0+tanDir': 's0' combined with difference in
%                               tangent directions (angle) between two
%                               iterations
%                           > 'combAll': Criteria 's0' AND 'sysr' must be
%                               smaller than tolerance
%                           > 'combAny': Criteria 's0' OR 'sysr' must be
%                               smaller than tolerance
%                           [{'s0'} / 'sysr' / 's0+tanDir' / 'combAll' /...
%                           'combAny]
%           - .degTol:      Defines the angle tolerance for stopping
%                           criterion 's0+tanDir',
%                           [{5} / positive double]
%           - .interactive: Enables the interactive mode which allows to
%                           stop and plot the result of every iteration
%                           [{false} / true]
%           - .initShifts:  Sets method for initial shift selection
%                           [{'eig_circle'} / 'zeros' / 'linear' / ...
%                           'logarithmic' / 'eig_large', 'eig_small', 'diag']
%                           See 'initShifts' for detailed information.
%           - .verbose:     Turn warnings and messages on by setting to
%                           true
%                           [{false} / true]
%           - .summary:     Enable/disable final summary of the algorithm
%                           [{true} / false]
%           - .arnoldiPH.*  Other options that will be passed on to the
%                           used method 'arnoldiPH';
%                           Plase refer to documentation of the respective
%                           algorithm.
%
% Output Arguments:
%       - redSys: 	reduced model (phsRed object)
%       - V:        Matrix that is used for reduction
%       - s0:       Shifts of the last iteration
%       - b:        Tangent directions of the last iteration
%       - W:        Matrix that is used for reduction
%       - nLU:      Total amount of LU decompositions
%
% Examples:
%       see demo_irkaPH
%
% See Also:
%       phs, arnoldiPH, initShifts, demo_irkaPH
%
% References:
%       [1] S. Gugercin et al. “Structure-Preserving Tangential Interpolation 
%           for Model Reduction of Port-Hamiltonian Systems." 
%           In: Automatica 48.9 (2012), pp. 1963–1974.
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

narginchk(2,5);

[s0, b, Opts] = parseInputs(sys, redOrder, varargin{:});

% Prepare command line outputs
stopCritLabel = {Opts.stopCrit};
if strcmp(stopCritLabel,'combAll')||strcmp(stopCritLabel,'combAny')
    stopCritLabel = {'s0','sysr'};
end
if strcmp(stopCritLabel,'s0+tanDir')
    stopCritLabel = {'s0','tanDir'};
end

% Suppress warnings
warning('off', 'MORpH:phs:noInputValidation') % Input validation in the end!


%% IRKAPH
% Initialization
iterations_to_prompt = 0;   % iterations until prompt for interactive mode; not required in other modes
nLU = 0;                    % accumulated number of LU decompositions

switch Opts.stopCrit
    case 'combAny'
        nStopVal = 2;
    case 'combAll'
        nStopVal = 2;
    case 's0+tanDir'
        nStopVal = 2;
    otherwise
        nStopVal = 1;
end
stopCrit_evolution = zeros(Opts.maxIter, nStopVal); % keep track of stopping criterion evolution

iteration = 0;      % iteration counter
stop = false;       % initialize stop
redSys = [];        % initialize empty reduced order system

% if interactive mode: plot original system
if Opts.interactive
    [stop, iterations_to_prompt] = interactive_prompt(sys, redSys, iteration, iterations_to_prompt);
end

% ----------------------------- IRKA ITERATION -----------------------------
while (~stop && (iteration < Opts.maxIter || Opts.interactive))

    iteration = iteration + 1;
    if Opts.verbose
        fprintf("Iteration %i\n", iteration);
    end

    % Determine next reduced system
    redSys_old = redSys;
    [redSys, V, W, nLUi] = arnoldiPH(sys, s0, b, Opts.arnoldiPH);
    nLU = nLU + nLUi;

    % Determine mirrored eigenvalues (and left eigenvectors if MIMO)
    s0_old = s0;
    b_old = b;
    if ~sys.isMIMO
        % SISO system
        s0 = - redSys.eig();
        b = ones(1, length(s0));
    else
        % MIMO system
        [~, L, Y] = redSys.eig();
        s0 = -diag(L)';
        b = (redSys.G-redSys.P)'*Y;
    end

    % Check if algorithm converged
    [stop,stopCrit] = stoppingCriterion(s0,s0_old,redSys,redSys_old,Opts,b,b_old);
    stopCrit_evolution(iteration,:) = stopCrit;

    % Output to console
    if Opts.verbose
        for i=1:length(stopCritLabel)
            fprintf(" Stopping criterion (%s): %.3e\n",stopCritLabel{i},stopCrit(i));
        end
    end

    % Interactive mode --> plot and ask user
    if Opts.interactive
        [stop, iterations_to_prompt] = interactive_prompt(sys, redSys, iteration, iterations_to_prompt);
    end

end
% ------------------------- END OF IRKA ITERATION --------------------------

% Check phs structure of system if not already done
warning('on', 'MORpH:phs:noInputValidation')
if ~redSys.Opts.inputValidation % ~Opts.irkaPH.arnoldiPH.phs.inputValidation
    redSys.Opts.inputValidation = true;
    phs.inputValidation(redSys);
end

% Reset s0 to s0_old for correct output
% (shifts that were used for creating the reduced model)
s0 = s0_old;
b = b_old;

% Change information of the reduced system to irka-information
redSys.method = @irkaPH;
arnoldiPHParameters = redSys.parameters;
arnoldiPHInfo = redSys.info;
redSys.parameters = Opts;
redSys.parameters.arnoldiPH = arnoldiPHParameters;
redSys.info = struct();
redSys.info.arnoldiPH = arnoldiPHInfo;
redSys.info.iter = iteration;
redSys.info.stopCritEvolution = stopCrit_evolution(1:iteration, :);
redSys.info.nLU = nLU;

% Display final results
if Opts.summary
    if stop
        converged = 'yes';
    else
        converged = 'no';
    end

    fprintf("=====================================\n");
    fprintf(strcat(" IRKAPH stopped after ", num2str(iteration), " iterations\n"));
    fprintf(strcat(" Converged: ", converged, "\n"));
    for i=1:length(stopCritLabel)
        fprintf(" Stopping criterion (%s): %.3e\n",stopCritLabel{i},stopCrit(i));
    end
    fprintf("=====================================\n");
end
end


%% ======================== AUXILIARY FUNCTIONS ===========================

function [stop,stopCrit] = stoppingCriterion(s0,s0_old,redSys,redSys_old,Opts,varargin)
% Determine if stopping creterion (defined in Opts) is satisfied
% Syntax:
%   [stop,stopCrit] = stoppingCriterion(s0,s0_old,redSys,redSys_old,Opts)
%   [stop,stopCrit] = stoppingCriterion(s0,s0_old,redSys,redSys_old,Opts,Rt,Rt_old)

switch Opts.stopCrit
    case 's0'
        if length(s0) == length(s0_old)
            stopCrit = norm((s0 - s0_old)./s0_old, 1)/redSys.dim;
            stop = stopCrit <= Opts.tol;
        else
            % Dimensions are not equal, probably because some
            % dimensions were not linearly independent (see
            % structurePreservation)
            stopCrit = nan;
            stop = false;
        end
    case 'sysr'
        if isempty(redSys_old) || isempty(redSys)
            % One of the systems not available, e.g. in first irka iteration
            stop = false;
            stopCrit = nan;
        else
            stopCrit = norm(redSys-redSys_old)/norm(redSys);
            stop = stopCrit <= Opts.tol;
        end
    case 's0+tanDir'
        Rt = varargin{1};
        Rt_old = varargin{2};
        if length(s0) == length(s0_old) && length(Rt) == length(Rt_old)
            % Shift convergence
            stopCrit(1) = norm((s0 - s0_old)./s0_old, 1)/redSys.dim;
            stop(1) = stopCrit(1) <= Opts.tol;

            % Tangential directions
            angleRt = zeros(1,redSys.dim); % Initialization
            for iDir = 1:redSys.dim
                angleRt(iDir) = abs(rad2deg(subspace(Rt_old(:,iDir),Rt(:,iDir))));
            end
            stopCrit(2) = max(angleRt);
            stop(2) = stopCrit(2) <= Opts.degTol;

            % Tangential direction convergence
            stop = all(stop);
        else
            % Dimensions are not equal, probably because some
            % dimensions were not linearly independent (see
            % structurePreservation)
            stopCrit = [nan, nan];
            stop = false;
        end
    case 'combAll'
        Opts.stopCrit = 's0';
        [stop(1),stopCrit(1)] = stoppingCriterion(s0,s0_old,redSys,redSys_old,Opts,varargin);
        Opts.stopCrit = 'sysr';
        [stop(2),stopCrit(2)] = stoppingCriterion(s0,s0_old,redSys,redSys_old,Opts,varargin);
        stop = all(stop);
    case 'combAny'
        Opts.stopCrit = 's0';
        [stop(1),stopCrit(1)] = stoppingCriterion(s0,s0_old,redSys,redSys_old,Opts,varargin);
        Opts.stopCrit = 'sysr';
        [stop(2),stopCrit(2)] = stoppingCriterion(s0,s0_old,redSys,redSys_old,Opts,varargin);
        stop = any(stop);
    otherwise
        error('MORpH:irkaPH:unknownStopCrit','Stopping criterion not implemented!');
end
end

function [stop, iterations_to_prompt] = interactive_prompt(sys, redSys, iteration, iterations_to_prompt)
stop = false;
if iteration == 0
    figure();
    hold on
    bode(sys);
    set(findall(gcf,'type','line'),'linewidth',4)
    disp('Press any key to continue')
    pause
else
    hold on
    iterations_to_prompt = iterations_to_prompt - 1;
    if iterations_to_prompt <= 0
        set(findall(gcf,'type','line'),'LineStyle','--')
        bode(redSys);
        inp = input(['Press <Enter> to continue and \n<Q> to stop iteration '...
            'if you are satisfied with the result.\nYou may also '...
            'input the number of iterations until the next prompt.    '],...
            's');
        if isempty(inp)
            inp = 'Enter';
            iterations_to_prompt = 1;
        else
            if ~isnan(str2double(inp))
                iterations_to_prompt = str2double(inp);
            else
                iterations_to_prompt = 1;
            end
        end
        if min(inp == 'Q') || min(inp == 'q')
            hold off
            stop = true; % stop IRKA iteration
        end
    end
    hold off
end
end

function [s0, b, Opts] = parseInputs(sys, redOrder, varargin)
% Check phs input type
if ~isa(sys,'phs')
    error('MORpH:irkaPH:wrongInput', 'Original model is not an object of the phs-class');
end

% Catch incompatible DAE systems
if sys.isDAE
    error('irkaPH currently only supports pHODE systems.')
end

% Opts
if ~isempty(varargin) && isstruct(varargin{end})
    Opts = varargin{end};
    varargin(end) = [];
else
    Opts = struct();
end

% Option parsing
OptsAdmissible.tol = 1e-3;
OptsAdmissible.degTol = 5;  % for stopping criterion 's0+tanDir', tolerance in degree
OptsAdmissible.stopCrit = {'s0','sysr','s0+tanDir','combAll','combAny'};
OptsAdmissible.maxIter = 500;
OptsAdmissible.interactive = {false,true};
OptsAdmissible.verbose = {false,true};
OptsAdmissible.summary = {true,false};
OptsAdmissible.arnoldiPH.structurePreservation = 'specialInverse'; % only default - correct values will be checked in arnoldiPH
OptsAdmissible.arnoldiPH.phs.inputValidation = {false, true}; % default: only perform phs validation in the end
OptsAdmissible.initShifts = {'eig_circle', 'zeros', 'linear', 'logarithmic', 'eig_large', 'eig_small', 'diag'};

Opts = phsMOR_parseOpts(Opts,OptsAdmissible);

% s0
if isscalar(redOrder)
    % Only reduced order provided -> initialize shifts
    s0 = initShifts(sys, Opts.initShifts, redOrder);
elseif isvector(redOrder)
    % Shifts provided
    s0 = redOrder;
    redOrder = length(redOrder);
    Opts.initShifts = 'custom';
else
    error('MORpH:irkaPH:badInputPattern', ...
        ['You may either provide a reduced order (scalar) or a vector '...
        'of initial shifts as second input to irkaPH.']);
end

if redOrder > sys.dim
    error('MORpH:irkaPH:wrongInput', ...
        'Reduced Order cannot exceed original system order')
end

% b
switch length(varargin)
    case 0
        % b is initialized as one-vectors
        b = ones(size(sys.G,2), redOrder);
    case 1
        b = varargin{1};
    otherwise
        error('MORpH:irkaPH:badInputPattern',...
            'The number of inputs is not correct.')
end

Opts.startShifts = s0;
Opts.startTangent = b;

end