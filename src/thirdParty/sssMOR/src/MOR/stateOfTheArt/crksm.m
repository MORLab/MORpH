function [sysr,V,W,Z,data] = crksm(varargin)
% CRKSM - Cumulative Rational Krylov Subspace Method for cumulative reduction and/or approximately solving Lyapunov equations 
%
% Syntax:
%       [sysr,V,W,S,data]             = CRKSM(sys, s0_inp)
%       [sysr,V,W,S,data]             = CRKSM(sys, s0_inp, Rt) 
%       [sysr,V,W,R,data]             = CRKSM(sys, [], s0_out) 
%       [sysr,V,W,R,data]             = CRKSM(sys, [], s0_out, [], Lt)
%       [sysr,V,W,S,data]             = CRKSM(sys, s0_inp, s0_out)  
%       [sysr,V,W,S,data]             = CRKSM(sys, s0_inp, s0_out, Rt, Lt)
%       [sysr,V,W,R,data]             = CRKSM(sys.', s0_inp, s0_out)  
%       [sysr,V,W,R,data]             = CRKSM(sys.', s0_inp, s0_out, Rt, Lt) 
%       [sysr,...,data]               = CRKSM(sys,...,Opts_rksm)
%
% Description:
%       This function performs the Cumulative Rational Krylov Subspace Method  
%       (CRKSM) as proposed by Druskin and Simoncini in [1] and [2].
%
%       On the one hand, |crksm| can be used for cumulative model order reduction,
%       where the projection matrices V and W are computed iteratively and 
%       cumulated step by step using rational Krylov subspaces. Thus, the
%       order of the reduced model is ad hoc chosen by the algorithm until 
%       the tolerance for the MOR stopping criterion is reached.
%
%       On the other hand, |crksm| can also be used for approximately
%       solving linear Lyapunov equations:
%
%              $$    A X E^T + E X A^T + B B^T = 0 \quad   (1)    $$
%
%              $$    A^T Y E + E^T Y A + C^T C = 0 \quad   (2)    $$
%
%       To this end, the cumulatively computed projection matrices are used
%       to project (1) and (2) to obtain _reduced_ Lyapunov equations. These 
%       reduced Lyapunov equations are then solved using _direct_ techniques (MATLAB
%       built-in lyapchol using Hammarling's method). This way, the Cholesky 
%       factorization X=S*S' of the solution of (1) and the low rank factor 
%       Y = R*R' of the dual (generalized) equation (2) are computed. 
%
%       Auf symmetry eingehen, 'both' / 'control' / 'observe' eingehen
%
%
%       This function computes input (and output) Krylov subspaces
%       corresponding to the shifts s0_inp and s0_out. For MIMO systems,
%       tangential directions Rt and Lt can be defined. Otherwise, block
%       Krylov subspaces will be computed.
%
%       Auf MIMO tangential, MIMO block eingehen. Auf adaptiveShifts
%       eingehen
%
% Input Arguments:
%		*Required Input Arguments:*
%       -sys:                   sss-object containing LTI system 
%       -A/B/C/D/E:             system matrices
%       -s0_inp:                initial expansion points for input Krylov 
%                               subspace (must have at least two entries)
%
%		*Optional Input Arguments:*
%       -s0_out:                initial expansion points for output Krylov subspace
%       -Rt/Lt:                 initial right/left tangential directions (MIMO case)
%       -Opts:                  a structure containing following options
%       *Option Settings for the whole crksm function*
%           -.purpose:          purpose of using CRKSM [{'lyapunov'} / 'MOR']
%                -'lyapunov':       use CRKSM for approximately solving Lyapunov equations
%                -'MOR':            use CRKSM for cumulative and adaptive model order reduction
%           -.maxiter:          specify maximal number of iterations
%                               [{200} / positive integer]
%           -.lse:              use LU or hessenberg decomposition (refer to solveLse)
%                               [{'sparse'} / 'full' / 'hess' / 'iterative' / 'gauss']
%           -.orth:             orthogonalization of new projection direction (refer to arnoldi) 
%                               [{'2mgs'} / 'dgks' / 'mgs' / false]
%       *Option Settings for Lyapunov Equations Purposes*
%           -.equation:         specify Lyapunov equation to be solved inthe two-sided case
%                               [{'control'} / 'observe']
%           -.stopCritLyap:     specify stopping criteria [{'residualLyap'} / 'normChol']
%                -'residualLyap':   compute residual of Lyapunov equation
%                -'normChol':       compare the norm of the last two Cholesky factors
%           -.crksmNorm:        specify norm [{2} / 'fro']
%                -2:                use 2-norm (Euclidian Norm)
%                -'fro':            use Frobenius Norm
%           -.lowrank:          compute the low-rank factor of the final solution 
%                               [{0} / 1]
%           -.restolLyap:       tolerance for -.stopCrit = 'residualLyap'
%                               [{1e-8} / positive float]
%           -.rctol:            tolerance for -.stopCrit = 'normChol'
%                               [{1e-12} / positive float]
%       *Option Settings for MOR Purposes*
%           -.real:             keep the projection matrices real
%                               [{true} / false]
%           -.restolMOR:        tolerance for MOR stopping criterion comparing the last two reduced models sysr
%                               [{1e-3} / positive float]
%       *Option Settings for the choice of shifts*
%           -.shifts:           choose shifts and how they should be used [{'dynamical'} / 'fixedCyclic']
%                -'fixedCyclic':    same initial shifts defined by the user are cyclically used for the whole function
%                -'dynamical':      new shifts are obtained dynamically during the programme 
%           -.getShiftsStrategy:choose strategy for dynamically getting shifts with |getShifts| 
%                -'adaptive':       a new shift is computed online for every single iteration using Druskin's method [1,2]
%                -'eigs':           new shifts are computed online using the mirrored images of the reduced eigenvalues (Ritz values)
%           -.shiftTol:         tolerance for the choice of new shifts
%                               [{0.1} / positive float]
%
% Output Arguments:
%       -sysr:                  reduced system
%       -V,W:                   projection matrices spanning Krylov subspaces
%       -Z (S/R):               Cholesky factor X=S*S' of Lyapunov equation A*X*E'+E*X*A'+B*B'=0 or
%                               Cholesky factor Y=R*R' of Lyapunov equation A'*Y*E+E'*Y*A+C'*C=0
%       -data:                  struct containing the following data
%           -.Norm:             norm values of the iterations
%           -.Shifts_Input:     used input shifts s0_inp during the whole programme
%           -.Shifts_Output:    used output shifts s0_out during the whole programme
%
% Examples:
%       MOR_purpose: This code computes ....
%
%> sys = sss('fom');
%> Opts.purpose = 'MOR';
%> sysr = crksm(sys, -eigs(sys,8).');
%> bode(sys,'-',sysr,'--r');
%
%       Lyapunov_purpose: This code computes ....
%
% See Also: 
%       rk, arnoldi, solveLse, sss/lyapchol, mess_lradi, getShifts,
%       initializeShifts
%
% References:
%       * *[1] Druskin, Simoncini (2011)*, Adaptive Rational Krylov Subspaces
%       for large-scale dynamical systems
%       * *[2] Druskin, Simoncini, Zaslavsky (2014)*, Adaptive Tangential
%       Interpolation in Rational Krylov Subspaces for MIMO Dynamical Systems
%       * *[3] Kuerschner (2016)*, Efficient Low-Rank Solution of Large-Scale Matrix Equations
%       * *[4] Wolf (2014)*, H2 Pseudo-Optimal Moder Order Reduction
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sssMOR">sssMOR</a>, a Sparse State-Space, Model Order 
% Reduction and System Analysis Toolbox developed at the Chair of 
% Automatic Control, Technische Universitaet Muenchen. For updates 
% and further information please visit <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Paul Heidenreich, Maria Cruz Varona
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  30 Nov 2017
% Copyright (c) 2016-2017 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%% Create Def-struct containing default values
% Note: there may be no point named Def.reuseLU, otherwise one gets a conflict with solveLse/lyapchol/bilyapchol

% general default option settings (for the whole function)
Def.purpose              = 'lyapunov';        % [{'lyapunov'} / 'MOR']
Def.maxiter              = 200;               % default number of iterations
%**solveLse options
Def.lse                  = 'sparse';          % [{'sparse'} / 'full' / 'hess' / 'iterative' / 'gauss']
%**arnoldi options
Def.orth                 = '2mgs';            % [{'2mgs'} / 'dgks' / 'mgs' / false]

% default option settings for Lyapunov equation purposes
Def.equation             = 'control';         % [{'control'} / 'observe']
Def.stopCrit             = 'residualLyap';    % [{'residualLyap'} / 'normChol']
Def.crksmNorm            = 2;                 % [{2} / 'fro']
Def.lowrank              = 0;                 % [{0} / 1] 
Def.restolLyap           = 1e-8;              % default tolerance for -.stopCrit = 'residualLyap'
Def.rctol                = 1e-12;             % default tolerance for -.stopCrit = 'normChol' 

% default option settings for MOR
Def.real                 = true;              % [{true} / false]
Def.restolMOR            = 1e-4;              % default tolerance

% default option settings for the choice of shifts
Def.shifts               = 'dynamical';       % choose usage of shifts [{'dynamical'} / 'fixedCyclic']
Def.getShiftsStrategy    = 'adaptive';        % choose GETSHIFTS strategy for updating shifts: [ {'adaptive'} / 'eigs' ]
Def.shiftTol             = 0.1;               % default value for new shifts

% define data struct
data = struct('Norm',[],'Shifts_Input',[],'Shifts_Output',[]);

%% Parsing of Inputs

if isa(varargin{end},'struct') 
    Opts = varargin{end};
    varargin = varargin(1:end-1);
end

% create the options structure
if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end

% input of sys-object
if isa(varargin{1},'ss') || isa(varargin{1},'sss') || isa(varargin{1},'ssRed')

    % read in sys object
    sys = varargin{1};

    % use hess if sys is ssRed object
    if isa(sys,'ssRed')
        Opts.lse='hess'; 
        if isempty(sys.E), sys.E = eye(size(sys.A)); end %ssRed robust compatibility
    end

    % check usage and inputs
    if length(varargin) == 2                            % usage: CRKSM(sys, s0_inp)
        s0_inp = varargin{2};                s0_out  = [];        
        Rt     = [];                         Lt      = [];
        input  = 1;                          pointer = @blockV;
        Opts.hermite = false;        
    elseif length(varargin) == 3       
        if size(varargin{3},1) == size(sys.B,2) && sys.isSiso == 0 &&...
           (size(varargin{3},2) == size(varargin{2},2) || size(varargin{3},2) == size(varargin{2},1)) 
                                                        % usage: CRKSM(sys, s0_inp, Rt)
            s0_inp = varargin{2};            s0_out  = [];           
            Rt     = varargin{3};            Lt      = [];
            input  = 2;                      pointer = @tangentialV;
            Opts.hermite = false;            
        elseif isempty(varargin{2})                     % usage: CRKSM(sys, [], s0_out)
            s0_inp = [];                     s0_out  = varargin{3};  
            Rt     = [];                     Lt      = [];
            input  = 3;                      pointer = @blockW;
            Opts.hermite = false;            
        else                                            % usage: CRKSM(sys, s0_inp, s0_out)
            s0_inp = varargin{2};            s0_out  = varargin{3};      
            Rt     = [];                     Lt      = [];           
            input  = 4; 
            if size(sys.B,2) ~= size(sys.C,1)
                error('Block Krylov for m~=p is not supported in crksm');
            end
            if all(s0_inp == s0_out)
                % in hermite case, m must be equal to p (m = p)!!
                pointerV = @blockV;         Opts.hermite = true;
            else
                pointerV = @blockV;         pointerW = @blockW;     Opts.hermite = false;
            end           
        end 
    elseif length(varargin) == 5
        if isempty(varargin{2})                         % usage: CRKSM(sys, [], s0_out, [], Lt)
            s0_inp = [];                     s0_out  = varargin{3};  
            Rt     = [];                     Lt      = varargin{5};
            input  = 5;                      pointer = @tangentialW;
            Opts.hermite = false;            
        else                                            % usage: CRKSM(sys, s0_inp, s0_out, Rt, Lt)
            s0_inp = varargin{2};            s0_out  = varargin{3};  
            Rt     = varargin{4};            Lt      = varargin{5};   
            input  = 6;                      
            if all(s0_inp == s0_out)
                pointerV = @tangentialV;         Opts.hermite = true;
            else
                pointerV = @tangentialV;         pointerW = @tangentialW;     Opts.hermite = false;
            end
        end  
    else
        error('Input not compatible with current crksm implementation');
    end   
elseif length(varargin) > 1
    % here, possible implementation of matrix input (A,B,C,D,E)
end

% check shifts and tangential directions
if  (~exist('s0_inp', 'var') || isempty(s0_inp)) && ...
    (~exist('s0_out', 'var') || isempty(s0_out))
    error('sssMOR:rk:NoExpansionPoints','No expansion points assigned.');
end

% check if extended (EKSM) or rational Krylov (RKSM) case
if (input == 1 && s0_inp(1,1) == 0 && s0_inp(1,2) == inf) || ...
   (input == 3 && s0_out(1,1) == 0 && s0_out(1,2) == inf) || ...
   (input == 4 && s0_inp(1,1) == 0 && s0_inp(1,2) == inf && s0_out(1,1) == 0 && s0_out(1,2) == inf)
    Opts.shifts = 'fixedCyclic';
    fprintf('Extended Krylov Subspace Method (EKSM / K-PIK) will be used, shifts are used in a fixed cyclic way \n');  
else
    % sort expansion points & tangential directions
    s0_inp = shiftVec(s0_inp);
    s0_inp = single(s0_inp); % cast operation is important for roboustness of ismemeber command
    s0old = s0_inp;
    if Opts.real
            s0_inp = cplxpair(single(s0_inp));  %make sure shifts can be paired 
            s0_inp = double(s0_inp);
    else
        s0_inp = sort(s0_inp);
    end

    if ~isempty(Rt)
        if size(Rt,2) ~= length(s0_inp),         error('Inconsistent size of Rt');end
        if size(Rt,1) ~= size(sys.B,2),          error('Inconsistent size of Lt');end
        [~,cplxSorting] = ismember(s0_inp,s0old); 
        Rt = Rt(:,cplxSorting);
        % check if tangential directions are complex conjugated  
        if mod(size(find(imag(s0_inp)),2),2) ~= 0 && mod(size(find(imag(Rt(1,:))),2),2) ~= 0 &&...
           sum(sum(imag(Rt),2)) ~= 0
            error('wrong input, right tangential directions (Rt) must be in complex conjugate pairs');
        end
    end
    s0_inp = double(s0_inp); % re-cast operation

    if ~isempty(s0_out)
        % sort expansion points & tangential directions
        s0_out = shiftVec(s0_out);
        s0_out = single(s0_out); % cast operation is important for roboustness of ismemeber command
        s0old = s0_out;
        if Opts.real 
            s0_out = cplxpair(s0_out); %make sure shifts can be paired 
        else
            s0_out = sort(s0_out);
        end
        if ~isempty(Lt)
            if size(Lt,2) ~= length(s0_out),   error('Inconsistent size of Lt');end
            if size(Lt,1) ~= size(sys.C,1),    error('Inconsistent size of Lt');end
            [~,cplxSorting] = ismember(s0_out,s0old); 
            Lt = Lt(:,cplxSorting);
            % check if tangential directions are complex conjugated
            if mod(size(find(imag(s0_inp)),2),2) ~= 0 && mod(size(find(imag(Rt(1,:))),2),2) ~= 0 &&...
               sum(sum(imag(Rt),2)) ~= 0
                error('wrong input, left tangential directions (Lt) must be in complex conjugate pairs');
            end
        end
    end
    s0_out = double(s0_out); % re-cast operation
    if length(s0_inp) > size(sys.A,1) || length(s0_out) > size(sys.A,1)
        error('sssMOR:arnoldi:reducedOrderExceedsOriginal','The desired reduced order exceeds the original order');
    end

    if ~isempty(s0_inp) && ~isempty(s0_out)
        % check if number of input/output expansion points matches
        if length(s0_inp) ~= length(s0_inp)
            error('Inconsistent length of expansion point vectors.');
        end
    end
end
   
%% RKSM Method
% built first subspace (two columns in case of cplx. conj. shifts) with arnoldi
switch input
    case 1
        [basis1] = arnoldi(sys.E,sys.A,sys.B,s0_inp(1,1:2),Opts); % basis1 is V 
        Opts.equation = 'control';
    case 2
        [basis1] = arnoldi(sys.E,sys.A,sys.B,s0_inp(1,1:2),Rt(:,1:2),Opts); % basis1 is V
        Opts.equation = 'control';
    case 3 
        [basis1] = arnoldi(sys.E',sys.A',sys.C',s0_out(1,1:2),Opts); % basis1 is W
        Opts.equation = 'observe';
    case 4
        [basis1] = arnoldi(sys.E,sys.A,sys.B,s0_inp(1,1:2),Opts); % basis1 is V 
        [basis2] = arnoldi(sys.E',sys.A',sys.C',s0_out(1,1:2),Opts); % basis1 is W
    case 5
        [basis1] = arnoldi(sys.E',sys.A',sys.C',s0_out(1,1:2),Lt(:,1:2),Opts); % basis1 is W
        Opts.equation = 'observe';
    case 6
        [basis1] = arnoldi(sys.E,sys.A,sys.B,s0_inp(1,1:2),Rt(:,1:2),Opts); % basis1 is V 
        [basis2] = arnoldi(sys.E',sys.A',sys.C',s0_out(1,1:2),Lt(:,1:2),Opts); % basis1 is W
end

% preprocessing: initialize some variables for programme, set function handles
newdir1 = zeros(size(sys.A,1),size(basis1,2)/2); % declare newdir-variable
newdir2 = newdir1;
nablaNorm_last = 0;
if isempty(s0_inp)
    nShifts = size(s0_out,2);
else
    nShifts = size(s0_inp,2);
end

% check purpose of crksm-function
if strcmp(Opts.purpose,'lyapunov')
    Opts.didlyap = 0;
    usage = @crksmLyap; 
    if strcmp(Opts.stopCrit,'residualLyap')
        Opts.tol = Opts.restolLyap;  
    else
        Opts.tol = Opts.rctol;  
    end
    if strcmp(Opts.equation,'control')
       pointerLyap = @lyapS; 
    else
       pointerLyap = @lyapR;   
    end
    clearFields = {'restolMOR'};
    Opts = rmfield(Opts,clearFields);
else
    usage = @crksmSysr;
    pointerLyap = []; 
    Z = [];
    if Opts.hermite && ~strcmp(Opts.shifts,'fixedCyclic') && strcmp(Opts.getShiftsStrategy,'adaptive')
        Opts.getShiftsStrategy = 'eigs';
        disp('Shift update strategy changed from adaptive to eigs because adaptive does not support the hermite, two-sided case');
    end
    clearFields = {'restolLyap'};
    Opts = rmfield(Opts,clearFields);
end

clear input Def varargin s0old A B C D E 

% first reduction step
if  ~exist('basis2','var') % one-sided projection (only basis1 exists)
    Ar = basis1'*sys.A*basis1;   Br = basis1'*sys.B;   Er = basis1'*sys.E*basis1;   Cr = sys.C*basis1;
    basis2 = [];
else % two-sided projection (basis1 and basis2 exist)
    Ar = basis2'*sys.A*basis1;   Br = basis2'*sys.B;   Er = basis2'*sys.E*basis1;   Cr = sys.C*basis1;
end
% build ssRed object
sysr = ssRed(Ar,Br,Cr,sys.D,Er);

% call usage handle function for the first solving step
[sysr,Z,data,Opts] = usage(sys,sysr,basis1,1,s0_inp,s0_out,pointerLyap,data,Opts);

% start iteration
if ~exist('data.out2','var')
    for ii = (size(basis1,2)/size(newdir1,2))+1:1:Opts.maxiter
        if size(basis1,2) == (ii-1)*size(newdir1,2)

            % get new shifts and tangetial directions
            if (~isempty(s0_inp) && ii > size(s0_inp,2)) || (~isempty(s0_out) && ii > size(s0_out,2)) 
                % determine residual or cholesky norm development
                if ii <= 2*nShifts && isempty(Rt) && isempty(Lt)
                    nablaNorm = norm(gradient(data.Norm(1+size(data.Norm,1)-nShifts:end))); 
                    if nablaNorm < 1e-2
                        nablaNorm = norm(log10(gradient(data.Norm(1+size(data.Norm,1)-nShifts:end))))/10; 
                    end
                else
                    nablaNorm = 0;
                end
                % recycle shifts or get new shifts
                if strcmp(Opts.shifts,'fixedCyclic') || nablaNorm > nablaNorm_last 
                    s0_inp = repmat(s0_inp,1,2);        Rt = repmat(Rt,1,2); 
                    s0_out = repmat(s0_out,1,2);        Lt = repmat(Lt,1,2); 
                else
                    [s0_inp,s0_out,Rt,Lt] = getShifts(sys,sysr,nShifts,s0_inp,Rt,s0_out,Lt,basis1,basis2,Opts);
                end
                nablaNorm_last = nablaNorm;
            end

            % get new direction, enlarge basis
            if isempty(basis2) || isempty(basis1)
                newdir1 = pointer(sys,basis1,basis2,s0_inp,s0_out,Rt,Lt,ii,size(newdir1,2),Opts);
                if Opts.real == 1
                    basis1 = [basis1 real(newdir1)];       % basis1 is either V or W 
                else
                    basis1 = [basis1 newdir1];
                end
            else
                if Opts.hermite
                    [newdir1,newdir2] = pointerV(sys,basis1,basis2,s0_inp,s0_out,Rt,Lt,ii,size(newdir1,2),Opts);
                else
                    newdir1 = pointerV(sys,basis1,basis2,s0_inp,s0_out,Rt,Lt,ii,size(newdir1,2),Opts);
                    newdir2 = pointerW(sys,basis2,basis1,s0_inp,s0_out,Rt,Lt,ii,size(newdir1,2),Opts);
                end
                if Opts.real == 1
                    basis1 = [basis1 real(newdir1)];    % basis1 is V  
                    basis2 = [basis2 real(newdir2)];    % basis2 is W
                else
                    basis1 = [basis1 newdir1];    % basis1 is V  
                    basis2 = [basis2 newdir2];    % basis2 is W
                end  
            end                   
       end % end of complex wait-sequence (if isreal(vnew))

       % orthogonalize new basis; note: for loop is neccessary for block
       if isempty(basis2) || isempty(basis1)
           hermite_gram_sch = 0;
           for jj=size(basis1,2)-(size(newdir1,2)-1):1:size(basis1,2)
                basis1 = gramSchmidt(jj,basis1,hermite_gram_sch,Opts);
           end
            basis1 = gramSchmidt(1,basis1,hermite_gram_sch,Opts);
       else
           if Opts.hermite
               hermite_gram_sch = 1; 
           else
               hermite_gram_sch = 0;
           end
           for jj=size(basis1,2)-(size(newdir1,2)-1):1:size(basis1,2)
               [basis1,~,basis2] = gramSchmidt(jj,basis1,basis2,hermite_gram_sch,Opts);
           end
       end
       % reduction step
       jnew = size(basis1,2)-(size(newdir1,2)-1);    jnew_last = size(basis1,2);  
       if ~isempty(basis2)
           Ar(1:jnew-1,jnew:jnew_last)     = basis2(:,1:jnew-1)'*sys.A*basis1(:,jnew:jnew_last);              % 1. step: new columns to the existing matrix Ar      
           Ar(jnew:jnew_last,1:jnew_last)  = basis2(:,jnew:jnew_last)'*sys.A*basis1(:,1:jnew_last);           % 2. step: rows are filled up so that a reduced square matrix results
           Er(1:jnew-1,jnew:jnew_last)     = basis2(:,1:jnew-1)'*sys.E*basis1(:,jnew:jnew_last);
           Er(jnew:jnew_last,1:jnew_last)  = basis2(:,jnew:jnew_last)'*sys.E*basis1(:,1:jnew_last);
           Br(jnew:jnew_last,:)            = basis2(:,jnew:jnew_last)'*sys.B;
           Cr(:,jnew:jnew_last)            = sys.C*basis1(:,jnew:jnew_last);
       else
           Ar(1:jnew-1,jnew:jnew_last)     = basis1(:,1:jnew-1)'*sys.A*basis1(:,jnew:jnew_last);        
           Ar(jnew:jnew_last,1:jnew_last)  = basis1(:,jnew:jnew_last)'*sys.A*basis1(:,1:jnew_last);   
           Er(1:jnew-1,jnew:jnew_last)     = basis1(:,1:jnew-1)'*sys.E*basis1(:,jnew:jnew_last);
           Er(jnew:jnew_last,1:jnew_last)  = basis1(:,jnew:jnew_last)'*sys.E*basis1(:,1:jnew_last);
           Br(jnew:jnew_last,:)            = basis1(:,jnew:jnew_last)'*sys.B;
           Cr(:,jnew:jnew_last)            = sys.C*basis1(:,jnew:jnew_last);
       end
       sysr = ssRed(Ar,Br,Cr,sys.D,Er);     % ssRed-object

       % call usage handle function for Lyapunov/sysr
       [sysr,Z,data,Opts] = usage(sys,sysr,basis1,ii,s0_inp,s0_out,pointerLyap,data,Opts);

       % quit for loop, show information
       if ~isempty(data.Shifts_Input) || ~isempty(data.Shifts_Output), break;   end

       % enlarge subspace with the imaginary part of the direction and eventually calculate new real directions
       if Opts.real == 1 && (~isreal(newdir1) || ~isreal(newdir2))
           if isempty(basis2) || isempty(basis1)
               basis1 = [basis1 imag(newdir1)];
           else
               if ~isreal(newdir1) && ~isreal(newdir2)
                   basis1 = [basis1 imag(newdir1)];
                   basis2 = [basis2 imag(newdir2)];
               elseif ~isreal(newdir1) && isreal(newdir2)
                   basis1 = [basis1 imag(newdir1)];

                   % calculate new real direction, if last vnew is real and last wnew is complex 
                   newdir2 = pointerW(sys,basis2,basis1,s0_inp,s0_out,Rt,Lt,ii,size(newdir1,2),Opts);
                   basis2 = [basis2 real(newdir2)];

               elseif isreal(newdir1) && ~isreal(newdir2)
                   basis2 = [basis2 imag(newdir2)];

                   % calculate new real direction, if last vnew is real and last wnew is complex
                   newdir1 = pointerV(sys,basis1,basis2,s0_inp,s0_out,Rt,Lt,ii,size(newdir1,2),Opts);
                   basis1 = [basis1 real(newdir1)];
               end
           end
           % empty newdir of basis1 and/or basis2 
           newdir1 = zeros(size(sys.A,1),size(newdir1,2));
           newdir2 = newdir1;
       end
    end  % end of for loop
    
    % create output
    if nnz(sysr.C) && ~isempty(s0_out) && isempty(basis2)
        W = basis1;  V = W;
    elseif isempty(basis2) && isempty(s0_out)
        V = basis1;  W = V;
    else
        V = basis1;  W = basis2;
    end
    % build final sysr-object
    Opts.originalOrder = sys.n;
    sysr = ssRed(sysr.A,sysr.B,sysr.C,sysr.D,sysr.E,'crksm',Opts,sys);
    
    if ii == Opts.maxiter
        warning('maximum number of iterations is reached without converging!' )
        if ~isempty(s0_inp), data.Shifts_Input  = s0_inp;    end     % this line is important for leaving the for loop
        if ~isempty(s0_out), data.Shifts_Output = s0_out;    end 
    end
end % end of: "if ~exist('data.out2','var')"
end


%% ***************************** AUXILIARY ********************************
function [V, TRv, W, TLw] = gramSchmidt(jCol,V,varargin)
%   Gram-Schmidt orthonormalization
%   Input:  jCol:  Column to be treated
%           V, W:  Krylov-Subspaces
%   Output: V, W:  orthonormal basis of Krylov-Subspaces
%           TRv, TLw: Transformation matrices

% input
if isa(varargin{end},'struct')
    Opts = varargin{end};
    varargin = varargin(1:end-1);
end

if length(varargin) == 1
    hermite_gram_sch = varargin{1};
    W = eye(size(V));
elseif length(varargin) == 2
    W = varargin{1};
    hermite_gram_sch = varargin{2};
end
    
% IP (inner product) is defined here
IP = @(x,y) (x.'*y);

TRv=eye(size(V,2));
TLw=eye(size(V,2));
if jCol>1
    switch Opts.orth
        case 'dgks'
            % iterates standard gram-schmidt
            orthError=1;
            count=0;
            while(orthError>Opts.dgksTol)
                h=IP(V(:,1:jCol-1),V(:,jCol));
                V(:,jCol)=V(:,jCol)-V(:,1:jCol-1)*h;
                TRv(:,jCol)=TRv(:,jCol)-TRv(:,1:jCol-1)*h;
                if hermite_gram_sch
                    h=IP(W(:,1:jCol-1),W(:,jCol));
                    W(:,jCol)=W(:,jCol)-W(:,1:jCol-1)*h;
                    TLw(:,jCol)=TLw(:,jCol)-TLw(:,1:jCol-1)*h;
                end
                orthError=norm(IP([V(:,1:jCol-1),V(:,jCol)/sqrt(IP(V(:,jCol),V(:,jCol)))],...
                    [V(:,1:jCol-1),V(:,jCol)/sqrt(IP(V(:,jCol),V(:,jCol)))])-speye(jCol),'fro');
                if count>50 % if dgksTol is too small, Matlab can get caught in the while-loop
                    error('Orthogonalization of the Krylov basis failed due to the given accuracy.');
                end
                count=count+1;
            end
        case 'mgs'
            for iCol=1:jCol-1
              h=IP(V(:,jCol),V(:,iCol));
              V(:,jCol)=V(:,jCol)-V(:,iCol)*h;
              TRv(:,jCol)=TRv(:,jCol)-h*TRv(:,iCol);
              if hermite_gram_sch
                h=IP(W(:,jCol),W(:,iCol));
                W(:,jCol)=W(:,jCol)-W(:,iCol)*h;
                TLw(:,jCol)=TLw(:,jCol)-h*TLw(:,iCol);
              end 
            end
       case '2mgs'
            for k=0:1
                for iCol=1:jCol-1
                  h=IP(V(:,jCol),V(:,iCol));
                  V(:,jCol)=V(:,jCol)-V(:,iCol)*h;
                  TRv(:,jCol)=TRv(:,jCol)-h*TRv(:,iCol);
                  if hermite_gram_sch
                    h=IP(W(:,jCol),W(:,iCol));
                    W(:,jCol)=W(:,jCol)-W(:,iCol)*h;
                    TLw(:,jCol)=TLw(:,jCol)-h*TLw(:,iCol);
                  end 
                end
            end
        otherwise
            error('Opts.orth is invalid.');
    end  
end

% normalize new basis vector
h = sqrt(IP(V(:,jCol),V(:,jCol)));
V(:,jCol)=V(:,jCol)/h;
TRv(:,jCol) = TRv(:,jCol)/h;
if hermite_gram_sch
    h = sqrt(IP(W(:,jCol),W(:,jCol)));
    W(:,jCol)=W(:,jCol)/h;
    TLw(:,jCol) = TLw(:,jCol)/h;
end
end

function [ Xnorm ] = NormFrobEfficient( Y,varargin )
% NormFrobEfficient: an efficient way of calculating the Frobenius Norm of
% a symmetric matrix X

% Detailed explanation
% Input: factors G, D of the low rank factorization of X=G*D*G'
% Output: Frobenius Norm x (scalar)
%
% Remark:
% If the Frobenius Norm of e.g. X=B*B' is needed then set G=B and D=I
% because X=G*D*G' = B*I*B' = B*B'.
% For detailed explanation on low rank factorization see paper "Efficient
% low rank solution of generalized Lyapunov equations"

% read in input
if nargin == 2
    X = varargin{1};
    if istril(X) || istriu(X)
        X_L = X;
        X_D = speye(size(X_L,1));
    end
else
    X_L = varargin{2};
    X_D = varargin{3};
end

% if no matrix decomposition is available
if ~exist('X_L','var')
    try    
        X_L = chol(X);
        X_D = speye(size(X));
    catch
        try 
            [X_L,X_D] = ldl(X);
        catch
            % if all other matrix decompositions failed use only the
            % symmetric part of X
            [X_L,X_D] = ldl((X+X')/2);
        end
    end
end
% QR-decomposition of X_L
[~,X_R] = qr(full(X_L),0);

% determine Frobenius norm: || X ||_F
if Y == 1
    Xnorm = norm(X_R*X_D*X_R');
else
    Xnorm = norm(X_R*X_D*X_R','fro');
end
end

% function for building new (block) directions for V-basis or for V and W-basis in hermite case
function [vnew,wnew] = blockV(sys,V,W,s0_inp,~,~,~,iter,colIndex,Opts)
    rhsB = V(:,size(V,2)-(colIndex-1):size(V,2));
    % preallocate memory, set unity matrix for tangential directions, set options
    Rt = eye(size(sys.B,2)); 
    vnew = zeros(size(sys.A,1),size(sys.B,2));
    if Opts.hermite
        wnew = zeros(size(sys.A,1),size(sys.C,1));
        Lt = eye(size(sys.C,1));
        rhsC = W(:,size(W,2)-(colIndex-1):size(W,2));
    end
    if s0_inp(1,iter) == s0_inp(1,iter-1)
        Opts.reuseLU = 1;
    else
       Opts.reuseLU = 0; 
    end
    % make SISO-system and calculate b-block column wise
    for ii = 1:1:size(sys.B,2)
        rhsB_ii = sys.E*rhsB; 
        if Opts.hermite
            rhsC_ii = sys.E'*rhsC;
            [v_ii,w_ii] = solveLse(sys.A,rhsB_ii,rhsC_ii',sys.E,s0_inp(1,iter),Rt(:,ii),Lt(:,ii),Opts);
            w_ii = w_ii(:,1);
            wnew(:,ii) = w_ii;
        else
            [v_ii] = solveLse(sys.A,rhsB_ii,sys.E,s0_inp(1,iter),Rt(:,ii),Opts);
            wnew = [];
        end
        v_ii = v_ii(:,1);  
        vnew(:,ii) = v_ii; 
        Opts.reuseLU = 1;
    end 
end

function [wnew] = blockW(sys,W,~,s0_inp,s0_out,~,~,iter,colIndex,~)
    rhsC = W(:,size(W,2)-(colIndex-1):size(W,2));
    % preallocate memory, set unity matrix for tangential directions, set options
    Lt = eye(size(sys.C,1)); 
    wnew = zeros(size(sys.A,1),size(sys.C,1));
    if (~isempty(s0_inp) && s0_inp(1,iter) == s0_out(1,iter)) || s0_out(1,iter) == s0_out(1,iter-1)
        Opts.reuseLU = 1;
    else
       Opts.reuseLU = 0; 
    end
    % make SISO-system and calculate b-block column wise
    for ii = 1:1:size(sys.C,1)
        rhsC_ii = sys.E'*rhsC;
        [w_ii] = solveLse(sys.A',rhsC_ii,sys.E',s0_out(1,iter),Lt(:,ii),Opts);
        if size(w_ii) > 1,  w_ii = w_ii(:,1);   end
        wnew(:,ii) = w_ii; 
        Opts.reuseLU = 1;
    end 
end

function [vnew,wnew] = tangentialV(sys,~,~,s0_inp,~,Rt,Lt,iter,~,Opts)
    if s0_inp(1,iter) == s0_inp(1,iter-1)
        Opts.reuseLU = 1;
        if Opts.hermite
           [vnew,wnew] = solveLse(sys.A,sys.B,sys.C,sys.E,s0_inp(1,iter),Rt(:,iter),Lt(:,iter),Opts); 
           wnew = wnew(:,1); 
        else
            vnew = solveLse(sys.A,sys.B,sys.E,s0_inp(1,iter),Rt(:,iter),Opts);
            wnew = []; 
        end
    else
        Opts.reuseLU = 0;
        if Opts.hermite
           [vnew,wnew] = solveLse(sys.A,sys.B,sys.C,sys.E,s0_inp(1,iter),Rt(:,iter),Lt(:,iter),Opts);
           wnew = wnew(:,1); 
        else
            vnew = solveLse(sys.A,sys.B,sys.E,s0_inp(1,iter),Rt(:,iter),Opts);
            wnew = []; 
        end
    end
    vnew = vnew(:,1);
end
    
function [wnew] = tangentialW(sys,~,~,s0_inp,s0_out,~,Lt,iter,~,Opts)
    if (~isempty(s0_inp) && s0_inp(1,iter) == s0_out(1,iter)) || s0_out(1,iter) == s0_out(1,iter-1)
        Opts.reuseLU = 1;
        wnew = solveLse(sys.A',sys.C',sys.E',s0_out(1,iter),Lt(:,iter),Opts);
    else
        Opts.reuseLU = 0;
        wnew = solveLse(sys.A',sys.C',sys.E',s0_out(1,iter),Lt(:,iter),Opts);
    end
    if size(wnew) > 1,  wnew = wnew(:,1);   end
end

function [sysr,Z,data,Opts] = crksmLyap(sys,sysr,basis1,iter,s0_inp,s0_out,pointerLyap,data,Opts)
    % solve controllability or observability Lyapunov equation
    [Z,Rnorm,Opts] = pointerLyap(sys,sysr,basis1,iter,Opts);
    if iter < 3,  data.Norm(1,1) = Rnorm; data.Norm(2,1) = Rnorm;  else, data.Norm(iter,1) = Rnorm;  end  
        
    % stop program
    if Rnorm < Opts.tol
       if Opts.lowrank == 1            % compute low rank factor of solution
           Z = basis1*(Z'); 
       end
       if ~isempty(s0_inp), data.Shifts_Input  = s0_inp;    end     % this line is important for leaving the for loop
       if ~isempty(s0_out), data.Shifts_Output = s0_out;    end 
       % show information of programme
       fprintf('CRKSM, usage Lyapunov, step:\t %d \t Convergence\n',iter);
    elseif size(Z,2) == size(sys.A,2)
       disp('The subspace has reached the dimension of the original system without converging!');
   end
end

function [sysr,Z,data,Opts] = crksmSysr(~,sysr,~,iter,s0_inp,s0_out,~,data,Opts)
    persistent stopCrit sysr_last 
    Z = [];
    if iter == 1, sysr_last = sss([],[],[]);  end
    % stopping criteria
    if all(real(eig(sysr))<0) && all(real(eig(sysr_last))<0)
       stopCrit = norm(sysr-sysr_last)/norm(sysr);
    else
       stopCrit = inf; %initialize in case the reduced model is unstable 
    end   
    if iter < 3,  data.Norm(1,1) = stopCrit; data.Norm(2,1) = stopCrit;  else, data.Norm(iter,1) = stopCrit;  end  
    
    if stopCrit < Opts.restolMOR
       if ~isempty(s0_inp), data.Shifts_Input  = s0_inp;    end     % this line is important for leaving the for loop
       if ~isempty(s0_out), data.Shifts_Output = s0_out;    end 
       fprintf('CRKSM, usage MOR, step: %d \t Convergence\n' ,iter);
    end
    sysr_last = sysr;
end

function [S,Rnorm,Opts] = lyapS(sys,sysr,basis1,iter,Opts)
persistent S_last
    if iter == 1, S_last = [];  end
    
    try
       S = lyapchol(sysr.A,sysr.B,sysr.E);
    catch
        % only show this warning once
        if Opts.didlyap == 0
            warning('Reduced system is unstable (iteration: %d), command "lyapchol" failed to solve for S',iter);
            fprintf('Programme continues solving the reduced Lyapunov equation with command "lyap" but an error may occur due to NaN or Inf entries in S \n');
            fprintf('For better stability behaviour try to perform crksm with one-sided projection only with V or W basis\n');
            fprintf('Try the call: [sysr,V,W,S,data] = CRKSM(sys, s0_inp) / [sysr,V,W,S,data] = CRKSM(sys, s0_inp, Rt)\n');
            fprintf('or [sysr,V,W,R,data] = CRKSM(sys, [], s0_out) / [sysr,V,W,R,data] = CRKSM(sys, [], s0_out, [], Lt)\n');
            Opts.didlyap = 1;
        end
        S = lyap(sysr.A,sysr.B*sysr.B',[],sysr.E);
        
    end
    % choose computation of norm/stopping criteria
    if strcmp(Opts.stopCrit,'residualLyap')
       % test determination (computation of residual after Panzer/Wolf)
       Er_inv_Br = solveLse(sysr.E,sysr.B);
       Opts.reuseLU = 0;
       Opts.Er_inv_Ar = solveLse(sysr.E,sysr.A,Opts);   

       % compute factors from residual
       Opts.B_ = sys.B-(sys.E*basis1)*Er_inv_Br;
       Cr_hat = solveLse(Opts.B_'*Opts.B_,Opts.B_'*((sys.A*basis1)-(sys.E*basis1)*Opts.Er_inv_Ar));
       F = sys.E*basis1*(Er_inv_Br+(S'*S)*Cr_hat');

       % compute residual norm (Euclidean Norm)
       if Opts.crksmNorm == 2
           res0  = norm(sysr.B' * sysr.B,2);
           Rnorm = max(abs(eig(full([Opts.B_'*Opts.B_+Opts.B_'*F, Opts.B_'*Opts.B_; F'*Opts.B_+F'*F, F'*Opts.B_])))) / res0; 
       else
           % Frobenius Norm
           res0  = norm(sysr.B' * sysr.B,'fro');
           Rnorm = sqrt(sum(eig(full([Opts.B_'*Opts.B_+Opts.B_'*F, Opts.B_'*Opts.B_; F'*Opts.B_+F'*F, F'*Opts.B_])))^2) / res0;
       end
    else
       if Opts.crksmNorm == 2
           X_lastnorm = NormFrobEfficient(1,S_last);
           X_norm = NormFrobEfficient(1,S);
       else
           X_lastnorm = NormFrobEfficient(0,S_last);
           X_norm = NormFrobEfficient(0,S);
       end
       Rnorm = abs(X_norm-X_lastnorm);
       S_last = S;
    end % end stopCrit
end

function [R,Rnorm,Opts] = lyapR(sys,sysr,basis1,iter,Opts)
persistent R_last 
    if iter == 1, R_last = [];  end
    
    try
        R = lyapchol(sysr.A',sysr.C',sysr.E');
    catch
        % only show this warning once
        if Opts.didlyap == 0
            warning('Reduced system is unstable (iteration: %d), command "lyapchol" failed to solve for R',iter);
            fprintf('Programme continues solving the reduced Lyapunov equation with command "lyap" but an error may occur due to NaN or Inf entries in R \n');
            fprintf('For better stability behaviour try to perform crksm with one-sided projection only with V or W basis\n');
            fprintf('Try the call: [sysr,V,W,S,data] = CRKSM(sys, s0_inp) / [sysr,V,W,S,data] = CRKSM(sys, s0_inp, Rt)\n');
            fprintf('or [sysr,V,W,R,data] = CRKSM(sys, [], s0_out) / [sysr,V,W,R,data] = CRKSM(sys, [], s0_out, [], Lt)\n');
            Opts.didlyap = 1;
        end
         R = lyap(sysr.A',sysr.C'*sysr.C,[],sysr.E');
    end
    
    % choose computation of norm/stopping criteria
    if strcmp(Opts.stopCrit,'residualLyap')
       % test determination (computation of residual after Panzer/Wolf), compute Er^-1*Br and Er^-1*Ar and other factors
       Er_invT_CrT = solveLse(sysr.E',sysr.C');
       Opts.reuseLU = 1;
       Opts.Er_invT_ArT = solveLse(sysr.E',sysr.A',Opts);   

       % compute factors from residual
       Opts.C_ = sys.C'-(sys.E'*basis1)*Er_invT_CrT;
       Cr_hat = solveLse(Opts.C_'*Opts.C_,Opts.C_'*((sys.A'*basis1)-(sys.E'*basis1)*Opts.Er_invT_ArT));
       F = sys.E'*basis1*(Er_invT_CrT+(R'*R)*Cr_hat');

       % compute residual norm (Euclidean Norm)
       if Opts.crksmNorm == 2
           res0  = norm(sysr.C * sysr.C',2);
           Rnorm = max(abs(eig(full([Opts.C_'*Opts.C_+Opts.C_'*F, Opts.C_'*Opts.C_; F'*Opts.C_+F'*F, F'*Opts.C_])))) / res0; 
       else
           % Frobenius Norm
           res0  = norm(sysr.B' * sysr.B,'fro');
           Rnorm = sqrt(sum(eig(full([Opts.C_'*Opts.C_+Opts.C_'*F, Opts.C_'*Opts.C_; F'*Opts.C_+F'*F, F'*Opts.C_])))^2) / res0;
       end
    else

       if Opts.crksmNorm == 2
           X_lastnorm = NormFrobEfficient(1,R_last);
           X_norm = NormFrobEfficient(1,R);
       else
           X_lastnorm = NormFrobEfficient(0,R_last);
           X_norm = NormFrobEfficient(0,R);
       end
       Rnorm = abs(X_norm-X_lastnorm);
       R_last = R;
    end % end stopCrit
end


