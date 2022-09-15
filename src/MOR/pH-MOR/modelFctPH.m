function [sysm,s0mTot,RtmTot,V,W,nLU] = modelFctPH(sys,s0m,varargin)
% MODELFCPH - Computes or updates the model function of a phs object
%
% Syntax:
%       sysm = MODELFCTPH(sys,s0m)
%       [sysm, s0mTot, V, W]        = MODELFCTPH(sys,s0m,Opts)
%       [sysm, ...]                 = MODELFCTPH(sys,s0m,Rtm,Opts)
%       [sysm, s0mTot, V, W]        = MODELFCTPH(sys,s0m,s0mTot,V,W,Opts)
%       [sysm, ...]                 = MODELFCTPH(sys,s0m,Rtm,s0mTot,RtmTot,V,W,Opts)
%       [sysm, ...]                 = MODELFCTPH(sys,...,Opts)
%       [sysm, s0mTot, V, W, nLU]   = MODELFCTPH(sys,s0m,...)
%
% Description:
%       This function generates a surrogate model, called "model function" sysm,
%       of a high-dimensional model sys by interpolation at the complex
%       frequencies defined in the vector s0m.
%
%       If additional inputs s0mTot, V, W are passed, then the model function
%       is updated combining the already existing information, defined by
%       s0mTot, to s0m. V and W are updated respectively.
%
%       The optional structure Opts allows the definition of additional
%       execution parameters.
%
%       This function is used in the model function based MOR approach
%       introduced in [1] and [2].
%
%
% Input Arguments:
%       *Required Input Arguments:*
%       -sys:			full oder model (phs object)
%       -s0m:			vector of shifts for model function (update)
%       -Opts:			structure with execution parameters
%            -.updateModel: chooses which shifts of s0m are included in update;
%                           [{'new'} / 'all']
%            -.modelTol:    tolerance for identifying new shifts;
%                           [{1e-3} / positive float]
%            -.degTol  :    convergence tolerance for subspace angles (deg);
%                           [{5} / positive float]
%            -.structurePreservation: See structurePreservation.m
%            -.plot    :    generate analysis plots;
%                           [{false} / true]
%           - .phs.*        Other options that will be passed on to the class 'phs';
%                           Please refer to documentation (phs).
%           - .arnoldiPH.*  Other options that will be passed on to the
%                           used method 'arnoldiPH';
%                           Plase refer to documentation of the respective
%                           algorithm.
%
%       *Optional Input Arguments:*
%       -s0mTot:        vector of all shifts already used
%       -Rtm:           matrix of tangent directions for model function
%                       (update)
%       -RtmTot:        matrix with all tangent directions already used as
%                       columns
%       -V,W:           corresponding projection matrices V,W
%
%
% Output Arguments:
%       -sysm:          (updated) model function;
%       -s0mTot:        cumulated vector of approximation frequencies
%       -RtmTot:        cumulated matrix of approximation tangent
%                       directions
%       -V,W:           (updated) projection matrices
%       -nLU:           number of LU decompositions
%
% See also:
%       cirkaPH, structurePreservation, arnoldiPH, demo_cirkaPH
%
% References:
%       [1] H. K. F. Panzer, “Model order reduction by Krylov subspace
%           methods with global error bounds and automatic choice of
%           parameters," Dissertation, Technical University of Munich,
%           Munich, 2014.
%       [2] A. Castagnotto, H. K. F. Panzer, and B. Lohmann, “Fast H2-
%           optimal model order reduction exploiting the local nature
%           of Krylov-subspace methods," in 2016 European Control
%           Conference (ECC), 2016, pp. 1958–1969.
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
[s0mTot, Rtm, RtmTot, V, W, Opts] = parseInputs(sys,s0m,varargin{:});

% New model function or update?
if isempty(s0mTot)
    % New model function
    s0mTot = rowVector(s0m);
    RtmTot = Rtm;
    V = []; W = [];
else
    % Model function update
    [s0m,Rtm] = updateModelFctShifts(s0mTot,s0m,RtmTot,Rtm,Opts);
    s0mTot = [rowVector(s0mTot), rowVector(s0m)];
    RtmTot = [RtmTot, Rtm];
end

%%  Compute the model function
% Check model function is not larger than original
if length(s0mTot)<size(sys.J,1)
    %   Update model function
    [sysm,V,W,nLU] = updateModelFct(sys,rowVector(s0m),Rtm,V,W,s0mTot,Opts);
else
    warning('MORpH:modelFct:sizeLimit',...
        ['Model function is already as big as the original.',...
        ' Returning the original model.']);
    sysm    = sys;
    V       = speye(sys.dim);
    W       = V;
    nLU     = 0;
end

end

% =================== AUXILIARY FUNCTIONS =====================
function [sysm,V,W,nLU] = updateModelFct(sys,s0,Rt,V,W,s0mTot,Opts)
% Compute the updated model function

warning('off', 'MORpH:phs:noInputValidation') % Input validation in the end!
if isempty(V)
    % First run: use Arnoldi
    [~, V, ~, nLU] = arnoldiPH(sys, s0, Rt, Opts.arnoldiPH);
else
    % Update
    if isempty(s0)
        warning('MORpH:modelFct:noUpdate',...
            's0 is empty, so no update will be performed. Try a lower tolerance.');
        %Leave V,W as they are
        nLU = 0;

    else
        if size(Rt,1) == 1
            %SISO
            [~, Vnew, ~, nLU] = arnoldiPH(sys, s0, Rt, Opts.arnoldiPH);
            V = [V,Vnew];
            [V,~] = qr(V,0);
        else
            % MIMO
            [~, Vnew, ~, nLU] = arnoldiPH(sys, s0, Rt, Opts.arnoldiPH);
            V = [V,Vnew];
            [V,~] = qr(V,0);
        end
    end
end
% Create model function from W and V --> phsRed object
[J_red, R_red, Q_red, G_red, E_red, P_red, S_red, N_red, V, W] = structurePreservation(sys,V,Opts);
warning('on', 'MORpH:phs:noInputValidation')
sysm = phsRed(J_red, R_red, Q_red, G_red, E_red, P_red, S_red, N_red, Opts.phs);
sysm.method = @modelFctPH;
sysm.parameters = Opts;
sysm.info.originalOrder = sys.dim;
sysm.info.s0mTot = s0mTot;
end % function updateModelFct

function [s0m,Rtm] = updateModelFctShifts(s0mTot,s0new,RtmTot,RtmNew,Opts)
% Distinguish between new and old shifts
switch Opts.updateModel
    case 'all'
        s0m = s0new;
        Rtm = RtmNew;

        % give robustness warning from MIMO
        if size(Rtm,1)> 1
            warning('MORpH:modelFct:updateAllMimo',...
                'The update option ''all''for MIMO models is not robust enough to cover higher multiplicities');
        end

    case 'new'
        % a) find shifts alredy used
        [idx,idxTot] = isMember(s0new,s0mTot,Opts.modelTol);
        % add new shifts and respective tangential directions
        s0m = s0new(~idx);
        Rtm = RtmNew(:,~idx);

        % b) find old shifts with different tangential directions
        if any(idx) && (size(Rtm,1) > 1) % only for MIMO
            idxOld   = find(idx);
            idxTot   = idxTot(idxOld);
            for iO = 1:length(idxOld)
                angR = abs(rad2deg(subspace(RtmNew(:,idxOld(iO)),RtmTot(:,idxTot(iO)))));
                if angR > Opts.degTol % new tangential direction
                    if Opts.verbose
                        fprintf(2,'A new tangential direction to old shift found\n');
                    end
                    s0m = [s0m, s0new(idxOld(iO))];
                    Rtm = [Rtm, RtmNew(:,idxOld(iO))];
                end
            end
        end

    otherwise
        error('selected model function update is not valid');
end
end % function updateModelFctShifts

function [idx, idxOld] = isMember(sNew,sOld,tol)
% ISMEMBER - Determine which elements are within a radius
% Returns boolean vector idx which indicates if shift at index is in sOld
% and vector of positions in old shift vector idxOld
idx     = zeros(1,length(sNew));
idxOld  = zeros(1,length(sNew));
for iS = 1:length(sNew)
    dist    = abs(sOld-sNew(iS));
    [mD,iD] = min(dist);
    if mD <= abs(sOld(iD))*tol
        idx(iS)     = true;
        idxOld(iS)  = iD;
    end
end
end

function [s0mTot, Rtm, RtmTot, V, W, Opts] = parseInputs(sys,s0m,varargin)
narginchk(2,8);

% Opts
if ~isempty(varargin) && isstruct(varargin{end})
    Opts        = varargin{end};
    varargin    = varargin(1:end-1);
else
    Opts = struct();
end

OptsAdmissible.updateModel = {'new','all'};
OptsAdmissible.modelTol = 1e-3;
OptsAdmissible.degTol = 5; %deg
OptsAdmissible.plot = {false,true};
OptsAdmissible.phs = struct();
OptsAdmissible.verbose = true;
OptsAdmissible.arnoldiPH = struct();
OptsAdmissible.arnoldiPH.structurePreservation = {false};
OptsAdmissible.arnoldiPH.phs.inputValidation = {false};

Opts = phsMOR_parseOpts(Opts,OptsAdmissible);

% Remaining parameters in varargin
switch length(varargin)
    case 0
        Rtm = ones(size(sys.G,2),length(s0m));
        s0mTot = [];
        RtmTot = [];
        V = [];
        W = [];
    case 1
        Rtm = varargin{1};
        s0mTot = [];
        RtmTot = [];
        V = [];
        W = [];
    case 3
        s0mTot  = varargin{1};
        V       = varargin{2};
        W       = varargin{3};
        RtmTot = [];
        Rtm = ones(size(sys.G,2),length(s0m));
    case 5
        Rtm     = varargin{1};
        s0mTot  = varargin{2};
        RtmTot  = varargin{3};
        V       = varargin{4};
        W       = varargin{5};
    otherwise
        error('MORpH:modelFct:nargin','Wrong number of input arguments')
end

end % function parseInputs
