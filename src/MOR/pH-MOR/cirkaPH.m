function [sysr, V, W, s0_old, Rt_old, sysm, s0mTot, nLU] = cirkaPH(sys, redOrder, varargin)
% CIRKA - Confined Iterative Rational Krylov Algorithm for port-Hamiltonian systems
%
% Syntax:
%       sysr                 = CIRKAPH(sys, redOrder)
%       sysr                 = CIRKAPH(sys, s0)
%       sysr                 = CIRKAPH(sys, s0, Rt)
%       sysr                 = CIRKAPH(sys, s0, ..., Opts)
%       [sysr, V, W, s0, Rt, sysm, s0mTot, nLU] = CIRKAPH(sys, ... )
%
% Description:
%       This function executes an adapted version of the Confined Iterative Rational Krylov
%       Algorithm (CIRKA) as proposed by Castagnotto et al. in [1,2].
%
%       The algorithm is based on constructing a model function, i.e. a
%       surrogate representing the full oder model locally about some
%       frequencies, and running IRKAPH with respect to the surrogate model.
%       The model function is updated until convergence.
%
%
% Input Arguments:
%       *Required Input Arguments:*
%       -sys:			full order model (phs object)
%       -s0:			vector of initial shifts
%       -redOrder:      (alternatively) reduced order
%       -Rt:            initial right tangential directions for MIMO
%
%       *Optional Input Arguments:*
%       -Opts:			structure with execution parameters
%           -.s0m:      initial shifts for surrogate;
%                       [{[s0,s0]} / vector ]
%           -.Rtm:      initial right tangential directions for surrogate;
%                       [{[Rt,Rt]} / matrix ]
%           -.maxIter:  maximum number of iterations;
%						[{20} / positive integer]
%           -.tol:		convergence tolerance;
%						[{1e-3} / positive double]
%           -.degTol:	convergence tolerance for subspace angles (deg);
%						[{5} / positive double]
%           -.stopCrit: convergence criterion for CIRKA (see also: irkaPH)
%                       ['s0' / 's0+tanDir' / 'sysr' / 'sysm' / {'combAny'} / 'combAll' ]
%           -.verbose:	show text output during iterations;
%						[{false} / true]
%           -.summary:  Show summary in command window after computation.
%                       [{true} / false]
%           -.plot:     plot results;
%                       [{false} / true]
%           -.suppressWarn: suppress warnings;
%                       [{false} / true]
%           -.clearInit: reset the model function after first iteration;
%                       [{true}, false]
%           -.modelFctPH.*  Other options that will be passed on to the
%                       used method (modelFctPH);
%                       Please refer to its documentation (doc modelFctPH).
%           -.irkaPH.*  Other options that will be passed on to the
%                       used method (irkaPH);
%                       Please refer to its documentation (doc irkaPH).
%
% Output Arguments:
%       -sysr:              reduced order model (phsRed object)
%       -V,W:               resulting projection matrices (V = Vm*Virka, W = Wm*Wirka)
%       -s0:                final choice of shifts
%       -Rt:                matrices of right tangential directions
%       -sysm:              resulting model function
%       -s0mTot:            shifts for model function
%       -nLU:               number of (high-dimensional) LU decompositions
%
% Examples:
%       sys = setup_MassSpringDamperSystem(100,2,1,1,'SISO')
%       sysr = cirkaPH(sys,20)
%
%       //Note: The computational advantage of the model function framework
%       is given especially for truly large scale systems, where the
%       solution of a sparse LSE becomes much more expensive than of a
%       small dense LSE. You can compare the advantage of CIRKA vs IRKA by
%       comparing the nLU output for the same reduction task.
%
% See Also:
%       irkaPH, modelFctPH, demo_cirkaPH
%
% References:
%       [1] A. Castagnotto, H. K. F. Panzer, and B. Lohmann, “Fast H2-
%           optimal model order reduction exploiting the local nature
%           of Krylov-subspace methods,” in 2016 European Control
%           Conference (ECC), 2016, pp. 1958–1969
%       [2] A. Castagnotto and B. Lohmann, “A new framework for
%           H2-optimal model reduction,” Mathematical and Computer
%           Modelling of Dynamical Systems, vol. 24, no. 3, pp. 236–257, 2018
%       [3] T. Moser, J. Durmann, and B. Lohmann. "SurrogateBased H2 Model 
%           Reduction of Port-Hamiltonian Systems." In: 2021 Eur. Control 
%           Conf. (ECC). 2021, pp. 2058–2065.
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
[s0_old, Rt_old, s0m, Rtm, Opts] = parseInputs(sys, redOrder, varargin{:});
startShifts = s0_old;
startTangent = Rt_old;

% Show warning if structure preservation methods do not coincide
if ~strcmp(Opts.modelFctPH.structurePreservation, Opts.irkaPH.arnoldiPH.structurePreservation)
    warning('MORpH:cirkaPH:differentStructurePreservation',...
        ['Opts.modelFctPH.structurePreservation is different to Opts.irkaPH.arnoldiPH.structurePreservation.\n',...
        'This is suboptimal. \nNote that the (optionally) returned matrices V and W are not correct any more!'])
end

% Suppress warnings
if Opts.suppressWarn
    warning('off','Control:analysis:NormInfinite3');
    warning('off','Control:analysis:InaccurateResponse');
end

%% CirkaPH algorithm
% Initialization
kIter = 0;                      % Outer iteration counter
kIrka = zeros(1,Opts.maxIter);  % Inner iteration counter
nLU = 0;                        % Number of high order nLU decompositions
switch Opts.stopCrit
    case 'combAny'
        nStopVal = 3;
    case 'combAll'
        nStopVal = 3;
    case 's0+tanDir'
        nStopVal = 2;
    otherwise
        nStopVal = 1;
end
stopVal = zeros(Opts.maxIter, nStopVal); % keep track of stopping Critierion evolution

% Initialize reduced models with zero systems (asymptotically stable for norm computations)
sysrOld = phs(zeros(length(s0_old)),1e-12*eye(length(s0_old)),eye(length(s0_old)),zeros(length(s0_old),size(sys.G,2)));
sysmOld = phs(zeros(length(Opts.s0m)),1e-12*eye(length(Opts.s0m)),eye(length(Opts.s0m)),zeros(length(Opts.s0m),size(sys.G,2)));

% Generate initial model function / intermediate model sysm
[sysm, s0mTot, RtmTot, Vm, Wm,nLUi] = modelFctPH(sys,s0m,Rtm,Opts.modelFctPH);
nLU = nLU + nLUi;

%% CIRKA iteration
if Opts.verbose, fprintf('Starting model function MOR...\n'); end
if Opts.plot, [~,~,sysFrd] = bode(sys); end

stop = false;
while ~stop && kIter < Opts.maxIter
    kIter = kIter + 1;
    if Opts.verbose, fprintf(sprintf('modelFctMor: k=%i\n',kIter));end

    % Update model function / intermediate model (sys --> sysm)
    if kIter > 1
        if kIter == 2 && Opts.clearInit
            % Reset the model function after the first step (optional)
            s0m = [rowVector(s0_old),rowVector(s0m(1:length(s0m)-length(s0_old)))];
            Rtm = [Rt_old,Rtm(:,1:length(s0m)-length(s0_old))];
            [sysm, s0mTot, RtmTot, Vm, Wm,nLUk]  = modelFctPH(sys,s0m,Rtm,Opts.modelFctPH);
        else
            [sysm, s0mTot, RtmTot, Vm, Wm, nLUk] = modelFctPH(sys,s0_old,Rt_old,s0mTot,RtmTot,Vm,Wm,Opts.modelFctPH);
        end
        nLU = nLU + nLUk; % Update count of high order nLU decompositions
    end

    % Reduce model with irkaPH (sysm --> sysr)
    [sysr, V_irka, s0, Rt, W_irka] = irkaPH(sysm, s0_old, Rt_old, Opts.irkaPH);
    kIrka(kIter) = sysr.info.iter;

    if Opts.plot
        figure; bodemag(sysFrd,sss(sys),ss(sysm),ss(sysr))
        legend('FOM','ModelFct','ROM');
        title(sprintf('kIter=%i, nModel=%i',kIter,sysm.dim));
        pause
    end

    % Check convergence
    [stop, stopVal(kIter,:)] = verifyStopCrit(s0_old, s0, sysr, sysrOld, ...
        Rt_old, Rt, sysm, sysmOld, Opts.stopCrit, Opts);

    % Detect full order
    if length(s0mTot)> sys.dim
        stop = true;
    end

    % Detect stagnation
    if ~stop && kIter > 1 && max(abs(stopVal(kIter-1,:)-stopVal(kIter,:))) < Opts.tol*1e-3
        stop = true;
    end

    % Update parameters for next iteration
    s0_old = s0;
    Rt_old = Rt;
    sysmOld = sysm;
    sysrOld = sysr;

    % Output to console
    if Opts.verbose
        fprintf(1,'\tkIrka: %03i\n',        kIrka(kIter));
        fprintf(1,'\tstopVal (%s): %s\n',   Opts.stopCrit,sprintf('%3.2e\t',stopVal(kIter,:)));
        fprintf(1,'\tModelFct size: %i \n', length(s0mTot));
    end
end % cirkaPH iteration

%% Set outputs
% Check phs structure of system if not already done
if ~sysr.Opts.inputValidation % ~Opts.irkaPH.arnoldiPH.phs.inputValidation
    sysr.Opts.inputValidation = true;
    phs.inputValidation(sysr);
end

% Overwrite irkaPH information
sysr.method = @cirkaPH;

parametersIrkaPH = sysr.parameters;
sysr.parameters = Opts;
sysr.parameters.irkaPH = parametersIrkaPH;

infoIrkaPH = sysr.info;
sysr.info = struct();
sysr.info.irkaPH = infoIrkaPH;

% Set cirkaPH parameters and info
sysr.parameters.startShifts = startShifts;      % store initial shifts
sysr.parameters.startTangent = startTangent;    % store initial tangent directions
sysr.info.s0 = s0_old;                  % store final shifts
sysr.info.Rt = Rt_old;                  % store final tangent directions
sysr.info.kIrka = kIrka;
sysr.info.originalOrder = sys.dim;
sysr.info.modelFctOrder = sysm.dim;
sysr.info.finalStopCrit = stopVal(kIter,:);
sysr.info.stopCritEvolution = stopVal(1:kIter,:);
sysr.info.relH2err = norm(sysm-sysr)/norm(sysm);
sysr.info.kIrka = kIrka(1:kIter);
sysr.info.nHighDimLU = nLU;
sysr.info.s0m = s0m;
sysr.info.Rtm = Rtm;

% Sort parameter and info fields
% (try/catch because orderfields does not work if struct names are not
% identical and complete)
try
    sysr.parameters = orderfields(sysr.parameters,...
        {'initShifts','startShifts','startTangent','s0m','Rtm','stopCrit',...
        'tol','degTol','maxIter', ...
        'modelTol','clearInit','stableModelFct',...
        'irkaPH','modelFctPH','verbose','summary','plot','suppressWarn'});
catch
    sysr.parameters = orderfields(sysr.parameters);
end
try
    sysr.info = orderfields(sysr.info,...
        {'s0','Rt','s0m','Rtm','finalStopCrit','stopCritEvolution',...
        'kIrka','nHighDimLU','originalOrder','modelFctOrder',...
        'relH2err','irkaPH'});
catch
    sysr.info = orderfields(sysr.info);
end

% Output matrices
V = Vm*V_irka;
W = Wm*W_irka;

% Display warnings and text output
if kIter >= Opts.maxIter && ~Opts.suppressWarn
    warning('sssMOR:cirka:maxiter','modelFctMor did not converge within maxiter');
end
if Opts.summary
    if stop     % stop = true indicates that algorithm converged
        fprintf('CIRKA step %03u - Convergence (%s): %s \n', ...
            kIter, Opts.stopCrit, sprintf('% 3.1e', stopVal(kIter,:)));
    else
        fprintf('CIRKA step %03u - not converged (%s): %s \n', ...
            kIter, Opts.stopCrit, sprintf('% 3.1e', stopVal(kIter,:)));
    end
end

% Reactivate warnings
if Opts.suppressWarn
    warning('on','Control:analysis:NormInfinite3');
    warning('on','Control:analysis:InaccurateResponse');
end

end

%% ===================== AUXILIARY FUNCTIONS ======================
function [s0, Rt, s0m, Rtm, Opts] = parseInputs(sys, redOrder, varargin)
narginchk(2,4);

% Check phs input type
if ~isa(sys,'phs')
    error('MORpH:cirkaPH:wrongInput', 'Original model is not an object of the phs-class');
end

% Check phs input type
if sys.isDAE
    error('MORpH:cirkaPH:wrongInput', 'CIRKA-PH currently only works for pHODE systems');
end

% Check if Opts is provided
if ~isempty(varargin) && isstruct(varargin{end})
    Opts = varargin{end};
    varargin = varargin(1:end-1);
else
    Opts = struct();
end

% Admissible value sets
OptsAdmissible.stopCrit= {'combAny','s0','s0+tanDir','sysr','sysm','combAll'}; %stopping criterion for CIRKA
OptsAdmissible.verbose         = {false,true};     % verbose mode cirkaPH
OptsAdmissible.summary         = {true, false};    % summary at end of computation
OptsAdmissible.plot            = {false,true};     % display text and plots
OptsAdmissible.suppressWarn    = {false,true};     % suppress warnings
OptsAdmissible.clearInit       = {true,false};     % reset the model fct after initialization?
OptsAdmissible.stableModelFct  = {true,false};     % make sysm stable
OptsAdmissible.irkaPH          = struct();         % make sure that irkaPH opts exists
OptsAdmissible.modelFctPH      = struct();         % make sure that modelFctPH opts exists
OptsAdmissible.initShifts = {'zeros', 'eig_circle', 'linear', 'logarithmic', 'eig_large', 'eig_small', 'diag'};
% Subfunctions
OptsAdmissible.modelFctPH.updateModel = {'new','all'}; % shifts used for the model function update
% Default values
OptsAdmissible.maxIter = 20;        % maximum number of CIRKA iterations
OptsAdmissible.tol     = 1e-3;      % tolerance for stopping criterion
OptsAdmissible.degTol  = 5;     	% tolerance for angles between directions
OptsAdmissible.modelTol= OptsAdmissible.tol;  %shift tolerance for model function
OptsAdmissible.modelFctPH.structurePreservation = 'specialInverse'; % Structure preservation model function
OptsAdmissible.irkaPH.arnoldiPH.structurePreservation = 'specialInverse';   % Structure preservation reduced model
OptsAdmissible.irkaPH.arnoldiPH.enforcePH = false;
OptsAdmissible.irkaPH.arnoldiPH.phs.inputValidation = false;    % Input validation will be performed in cirkaPH at the end
OptsAdmissible.irkaPH.tol = 1e-4;   % tolerance of irkaPH
% Suppress warnings from low-level functions to reduce console output
OptsAdmissible.irkaPH.summary                       = false;
OptsAdmissible.irkaPH.arnoldiPH.verbose             = false;
OptsAdmissible.irkaPH.arnoldiPH.phs.inputValidation = false;
OptsAdmissible.irkaPH.arnoldiPH.phs.verbose         = false;
OptsAdmissible.modelFctPH.phs.verbose               = false;

Opts = phsMOR_parseOpts(Opts,OptsAdmissible);

% s0
if isscalar(redOrder)
    % Only reduced order provided -> initialize shifts
    s0 = initShifts(sys, Opts.initShifts, redOrder);
elseif isvector(redOrder)
    % Shifts provided
    s0 = redOrder;
else
    error('MORpH:irkaPH:badInputPattern', ...
        ['You may either provide a reduced order (scalar) or a vector '...
        'of initial shifts as second input to irkaPH.']);
end

% Rt
switch length(varargin)
    case 0
        Rt = ones(size(sys.G,2),length(s0));
    case 1
        Rt = varargin{1};
end

if ~isfield(Opts, 's0m')
    Opts.s0m     = [rowVector(s0), rowVector(s0)];   % default surrogate shifts
end
if ~isfield(Opts, 'Rtm')
    Opts.Rtm     = [Rt,Rt]; 	% default right tangential directions for surrogate
end

s0m = Opts.s0m;
Rtm = Opts.Rtm;
end % function parseInputs

function [stop,stopVal] = verifyStopCrit(s0, s0new, sysr, sysrOld, Rt, Rtnew, sysm, sysmOld, stopCrit, Opts)
% Compute stopping criterion for CIRKAPH
s0 = rowVector(s0);
s0new = rowVector(s0new);
switch stopCrit
    case 's0'   % Shift convergence
        if length(s0) == length(s0new)
            stopVal = norm((s0 - s0new)./s0, 1)/sysr.dim;
            stop = stopVal <= Opts.tol;
        else
            % Dimensions are not equal, probably because some
            % dimensions were not linearly independent (see
            % structurePreservation)
            stopVal = nan;
            stop = false;
        end
    case 's0+tanDir' % Shift convergence + tangent direction convergence
        if length(s0) == length(s0new) && length(Rt) == length(Rtnew)
            % Shift convergence
            stopVal(1) = norm((s0 - s0new)./s0new, 1)/sysr.dim;
            stop(1) = stopVal(1) <= Opts.tol;

            % Tangential directions
            angleRt = zeros(1,sysr.dim); % Initialization
            for iDir = 1:sysr.dim
                angleRt(iDir) = abs(rad2deg(subspace(Rt(:,iDir),Rtnew(:,iDir))));
            end
            stopVal(2) = max(angleRt);
            stop(2) = stopVal(2) <= Opts.degTol;

            % Tangential direction convergence
            stop = all(stop);
        else
            % Dimensions are not equal, probably because some
            % dimensions were not linearly independent (see
            % structurePreservation)
            stopVal = [nan, nan];
            stop = false;
        end

    case 'sysr' % Reduced system convergence
        stopVal = norm(sysr-sysrOld)/norm(sysr);
        stop = (stopVal <= Opts.tol);

    case 'sysm' % Model system convergence
        stopVal=norm(sysm-sysmOld)/norm(sysm);
        stop = (stopVal <= Opts.tol);

    case 'combAll'
        [stop(1),stopVal(1)] = verifyStopCrit(s0, s0new, sysr, sysrOld, Rt, Rtnew, sysm, sysmOld, 's0', Opts);
        [stop(2),stopVal(2)] = verifyStopCrit(s0, s0new, sysr, sysrOld, Rt, Rtnew, sysm, sysmOld, 'sysr', Opts);
        [stop(3),stopVal(3)] = verifyStopCrit(s0, s0new, sysr, sysrOld, Rt, Rtnew, sysm, sysmOld, 'sysm', Opts);
        stop = all(stop);

    case 'combAny'
        [stop(1),stopVal(1)] = verifyStopCrit(s0, s0new, sysr, sysrOld, Rt, Rtnew, sysm, sysmOld, 's0', Opts);
        [stop(2),stopVal(2)] = verifyStopCrit(s0, s0new, sysr, sysrOld, Rt, Rtnew, sysm, sysmOld, 'sysr', Opts);
        [stop(3),stopVal(3)] = verifyStopCrit(s0, s0new, sysr, sysrOld, Rt, Rtnew, sysm, sysmOld, 'sysm', Opts);
        stop = any(stop);

    otherwise
        error('The stopping criterion selected is incorrect or not implemented');
end

end % function verifyStopCrit
