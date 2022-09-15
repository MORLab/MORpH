function [sysr,sysrVec] = cure(sys,Opts)
% CURE - CUmulative REduction framework
%
% Syntax:
%       sysr = CURE(sys)
%       sysr = CURE(sys,Opts)
%       [sysr, sysrVec] = CURE(sys,Opts)
%
% Description:
%       This function implements the CUmulative REduction framework
%       (CURE) introduced by Panzer, Wolf and Lohmann (see [1,2]).
%
%       Using the duality between Sylvester equation and Krylov subspaces, the 
%       error is factorized at each step of CURE and only the high-dimensional
%       factor is reduced in a subsequent step of the iteration.
%
%       Currently, this function supports following reduction strategies at
%       each step of CURE:
%       spark (def.), irka, rk+pork (pseudo-optimal reduction)
%
%       //Note: Currently CUREd SPARK works only for SISO systems.
%
% Input Arguments:
%       *Required Input Arguments:*
%       -sys: An sss-object containing the LTI system
%       *Optional Input Arguments:*
%       -Opts: A structure containing following fields
%           -.cure.redfun:  reduction algorithm
%                           [{'spark'} / 'irka' / 'rk+pork']
%           -.cure.nk:      reduced order at each iteration 
%                           [{'2'} / positive integer]
%           -.cure.fact:    factorization mode 
%                           [{'V'} / 'W']
%           -.cure.initMode:shift initialization mode 
%                           [{'zero'} / 'sm' / 'lm' / 'slm' / array]
%           -.cure.initN:   number of initial shifts
%                           [{stopval} (for stop='nmax'), {5*nk} (else),
%                           integer]
%           -.cure.stop:    stopping criterion
%                           [{'normROM'} / 'nmax' / 'h2Error']
%           -.cure.stopval: value according to which the stopping criterion is evaluated
%                           [{1e-6} / positive float]
%           -.cure.verbose: display text during cure 
%                           [{false} / true]
%           -.cure.test:    execute analysis code 
%                           [{false} / true]
%           -.cure.gif:     produce a .gif file of the CURE iteration
%                           [{false} / true]
%           -.cure.maxIter: maximum number of CURE iterations
%                           [{'20'} / positive integer]
%           -.warn:         show warnings
%                           [{false} / true]
%           -.w:            frequencies for analysis plots
%                           [{''} / '{wmin,wmax}' / vector of frequencies]
%           -.zeroThres:    value that can be used to replace 0 
%                           [{'1e-4'} / postivie float]
%
% Output Arguments:     
%       -sysr: Reduced system
%       -sysrVec: A cell array of all reduced models at each step
% Examples:
%       By default, cure reduces a given model sys to a reduced order of
%       sqrt(sys.n) by steps of nk = 2 using mespark (model function based
%       spark)
%> sys  = sss('building');
%> sysr = cure(sys); 
%> disp(sysr)
%
%       The behavior of the function can be highly customized using the
%       option structure Opts:
%
%> Opts.cure    = struct('nk',4, 'redfun', 'irka', 'verbose', 1, 'stopval',12);
%> sysr         = cure(sys,Opts); 
%> disp(sysr)
% 
% See Also: 
%       spark, rk, irka, porkV, porkW, getSylvester
%
% References:
%       * *[1] Panzer (2014)*, Model Order Reduction by Krylov Subspace Methods
%              with Global Error Bounds and Automatic Choice of Parameters
%       * *[2] Wolf (2014)*, H2 Pseudo-Optimal Moder Order Reduction
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
% Authors:      Heiko Panzer, Alessandro Castagnotto, Maria Cruz Varona,
%               Philipp Seiwald
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  23 Mar 2017
% Copyright (c) 2016,2017 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%% Parse input and load default parameters
    % default values
    Def.warn        = 0;%show warnings?
    Def.w           = []; %frequencies for bode plot
    Def.zeroThres   = 1e-4; % define the threshold to replace "0" by
        Def.cure.verbose    = 0; %show progress text?
        Def.cure.redfun     = 'spark'; %reduction algorithm
        Def.cure.nk         = 2; % reduced order at each step
        Def.cure.stop       = 'normROM'; %type of stopping criterion
        Def.cure.stopval    = 1e-6;
        Def.cure.initMode   = 'slm'; %shift initialization type
        Def.cure.initN      = 5*Def.cure.nk; %will be overwritten if stop=='nmax'
        Def.cure.fact       = 'V'; %error factorization
        Def.cure.test       = 0; %execute analysis code?
        Def.cure.gif        = 0; %produce a .gif of the CURE iteration
        Def.cure.maxIter    = 20; %maximum number of iterations
        Def.cure.checkEVB   = true; %check if [EV,B_] or dual have full rank
        Def.cure.sEVBTol    = 1e-16; %rank tolerance for [EV,B_] matrix (or dual)
        
    % create the options structure
    if ~exist('Opts','var') || isempty(Opts)
        Opts = Def;
    else
        Opts = parseOpts(Opts,Def);
    end              
    
    % make sure reduced order does not exceed original
    if Opts.cure.maxIter > sys.n/Opts.cure.nk
        Opts.cure.maxIter = floor(sys.n/Opts.cure.nk);
    end
    % initialize all shifts if Opts.cure.stop=='nmax'
    if strcmp(Opts.cure.stop,'nmax')
        Opts.cure.initN = Opts.cure.stopval;
    elseif Opts.cure.initN>sys.n
        %make sure the number of initial shifts is not greater than
        %original order
        warning('sssMOR:cure:numberOfInitialShifts','The number of initial shifts in CURE exceeds the full order and will be set to sys.n')
        Opts.cure.initN = sys.n;
    end
    
    % store name and the reductionParameters, if sys is of type ssRed
    name = sys.Name;
    reductionParameters = [];
    if isa(sys,'ssRed')
        reductionParameters = sys.reductionParameters;
    end
    
%%  Plot for testing
if Opts.cure.test
    fhOriginalSystem = figure('Name','CURE - Reduction of the original model');
    fhSystemBeingReduced = fhOriginalSystem; %the two coincide for the moment
    [m,w] = freqresp(sys,Opts.w);
    bode(frd(m,w),'-b'), hold on
    axH = findall(gcf,'type','axes');
    magHandle = axH(3); magLim = get(magHandle,'YLim');
    phHandle = axH(2); phLim  = get(phHandle,'YLim');
    drawnow
    if Opts.cure.gif, writeGif('create'); end
    drawnow
end
%%   Initialize some variables
[~,m] = size(sys.b);  p = size(sys.c,1);
Er_tot = []; Ar_tot = []; Br_tot = []; Cr_tot = []; B_ = sys.b; C_ = sys.c;
BrL_tot = zeros(0,p); CrL_tot = zeros(p,0); 
BrR_tot = zeros(0,m); CrR_tot = zeros(m,0);

sysr = sss(Ar_tot,Br_tot,Cr_tot,zeros(p,m),Er_tot);
Dr_tot = full(sys.d);

sysrVec = {}; % keep track of reduced models over iterations
%%   Start cumulative reduction
if Opts.cure.verbose, fprintf('\nBeginning CURE iteration...\n'); end

iCure = 0; %iteration counter
stop  = false;
while ~stop && iCure < Opts.cure.maxIter
    iCure = iCure + 1;
    %   Redefine the G_ system at each iteration
    sys = sss(sys.a,B_,C_,[],sys.e);
    
    %   Initializations
    [s0,Opts] = initializeShifts(sys,Opts,iCure);    
    if Opts.cure.verbose
        sStr = sprintf('%3.2e %+3.2ei \t %3.2e %+3.2ei',[real(s0(1)),imag(s0(1)),real(s0(2)),imag(s0(2))]);
        fprintf('\t\tstart shifts\t%s\n',sStr);
    end
    
	% 1) Reduction
    switch Opts.cure.fact
        case 'V'
            % V-based decomposition, if A*V - E*V*S - B*Rv = 0
            switch Opts.cure.redfun
                case 'spark'               
                    [sysrTemp,V,Sv,Rv,~] = spark(sys,s0,Opts);
                    [Ar,Br,Cr,~,Er] = dssdata(sysrTemp);
                case 'irka'
                    [sysrTemp,V,W,~,~,~,~,Sv,Rv] = irka(sys,s0);
                    [Ar,Br,Cr,~,Er] = dssdata(sysrTemp);
                    
                case 'rk+pork'
                    [sysrTemp, V, ~, ~, Sv, Rv] = rk(sys,s0);
                    [Ar,Br,Cr,Er] = porkV(V,Sv,Rv,C_);
                otherwise 
                    error('The selected reduction scheme (Opts.cure.redfun) is not availabe in cure');
            end
            n = size(V,2);
            Se = eig(Sv); % shifts for this iteration
            if Opts.cure.verbose
                sStr = sprintf('%3.2e %+3.2ei \t %3.2e %+3.2ei',[real(Se(1)),imag(Se(1)),real(Se(2)),imag(Se(2))]);
                fprintf('\t\tfinal shifts\t%s\n',sStr);
            end
        case 'W'
        % W-based decomposition, if A.'*W - E.'*W*Sw.' - C.'*Lw = 0
            switch Opts.cure.redfun
                case 'spark'               
                    Opts.spark.pork = 'W';
                    [sysrTemp,W,Sw,Lw,~] = spark(sys.',s0,Opts);
                    
                    [Ar,Br,Cr,~,Er] = dssdata(sysrTemp);
                case 'irka'
                    [sysrTemp,V,W,~,~,~,~,~,~,~,Sw,Lw] = irka(sys,s0);
                    
                    [Ar,Br,Cr,~,Er] = dssdata(sysrTemp);
                    
                case 'rk+pork'
                    [sysrTemp, ~, W, ~, ~, ~, ~, Sw, Lw] = rk(sys,[],s0);
                    
                    [Ar,Br,Cr,Er] = porkW(W,Sw,Lw,B_); 
                    
            end
            n = size(W,2);
            Se = eig(Sw); % shifts for this iteration
            if Opts.cure.verbose 
                sStr = sprintf('%3.2e %+3.2ei \t %3.2e %+3.2ei',[real(Se(1)),imag(Se(1)),real(Se(2)),imag(Se(2))]);
                fprintf('\t\tfinal shifts\t%s\n',sStr);
            end
    end
    
    %%  Cumulate the matrices and define sysr
	%Er = W.'*E*V;  Ar = W.'*A*V;  Br = W.'*B_;  Cr = C_*V;
    Er_tot = blkdiag(Er_tot, Er);
    Ar_tot = [Ar_tot, BrL_tot*Cr; Br*CrR_tot, Ar]; %#ok<*AGROW>
    Br_tot = [Br_tot; Br]; Cr_tot = [Cr_tot, Cr];
    
    if Opts.cure.fact=='V'
        B_ = B_ - sys.e*(V*(Er\Br));    % B_bot
        BrL_tot = [BrL_tot; zeros(n,p)];    BrR_tot = [BrR_tot; Br];
        CrL_tot = [CrL_tot, zeros(p,n)];    CrR_tot = [CrR_tot, Rv];
    elseif Opts.cure.fact=='W'
        C_ = C_ - Cr/Er*W.'*sys.e;		% C_bot
        BrL_tot = [BrL_tot; Lw.'];   BrR_tot = [BrR_tot; zeros(n,m)];
        CrL_tot = [CrL_tot, Cr];    CrR_tot = [CrR_tot, zeros(m,n)];
    end
    
    %%  Storing additional parameters
    %Stroring additional information about the reduction in the object 
    %containing the reduced model:
    %   1. Define a new field for the Opts struct and write the information
    %      that should be stored to this field
    %   2. Adapt the method "parseParamsStruct" of the class "ssRed" in such a
    %      way that the new defined field passes the check
      
    usedOpts = sysrTemp.reductionParameters(end).params;
    usedOpts.cure = Opts.cure;
    usedOpts.currentReducedOrder = sysr.n+Opts.cure.nk;
    usedOpts.originalOrder = sys.n;
    usedOpts.shifts = s0;

    
    if isa(sysr,'ssRed')
        sysr = ssRed(Ar_tot, Br_tot, Cr_tot, zeros(p,m), Er_tot, ...
                     strcat('cure_',Opts.cure.redfun),usedOpts,sysr);
    else            %first Iteration
        if ~isempty(reductionParameters)
            sysr = ssRed(Ar_tot, Br_tot, Cr_tot, zeros(p,m), Er_tot, ...
                         strcat('cure_',Opts.cure.redfun),usedOpts,reductionParameters);
        else
            sysr = ssRed(Ar_tot, Br_tot, Cr_tot, zeros(p,m), Er_tot, ...
                         strcat('cure_',Opts.cure.redfun),usedOpts);
        end
    end
    sysr.Name = name;

	%create a cell array of reduced models
    sysrVec{end+1} = sysr;
    
    % display
    if Opts.cure.test
        % frequency response
        sysr_bode = sysr; 
        figure(fhSystemBeingReduced);
        bode(sysr_bode,w,'--r');
        set(magHandle,'YLim', magLim); set(phHandle,'YLim', phLim);
        if Opts.cure.gif, writeGif('append'); end     
        drawnow
    end
    
    %% Evaluate stopping criterion
    [stop, stopCrit] = stoppingCriterion(sys,sysr,sysrVec,Opts);
    if length(Opts.cure.initMode)==Opts.cure.nk
        % initialization was only for one iteration
        % use the current (optimal) shifts for the next one
        Opts.cure.initMode = Se.';
    end
    
    if Opts.cure.verbose
        fprintf('CURE step %03u - Convergence (%s):\t%s \n',iCure,...
            Opts.cure.stop,sprintf('% 3.1e', stopCrit));
    end

end
%%   Was maxIter achieved?
if iCure >= Opts.cure.maxIter
    warning('Iterations count reached maxIter. You might want to increase maxIter or the convergence criterion')
end
%%   Add the feedthrough term before returning the reduced system
sysr.D = Dr_tot;
%%  Finishing execution
if Opts.cure.verbose
    fprintf('Stopping criterion satisfied. Exiting cure...\n\n');
else
    fprintf('CURE step %03u - Convergence (%s):\t%s \n',iCure,...
            Opts.cure.stop,sprintf('% 3.1e', stopCrit));
end
if Opts.cure.test
        sysr_bode = sysr;
        figure(fhOriginalSystem);
        bode(sysr_bode,w,'-g');
        drawnow
        title(sprintf('n_{red} = %i',size(sysr.a,1)));
        
        if Opts.cure.gif, writeGif('append'), end
end
%% --------------------------AUXILIARY FUNCTIONS---------------------------
function [stop,stopCrit] = stoppingCriterion(sys,sysr,sysrVec,opts)
%   computes the stopping criterion for CURE iteration
switch opts.cure.stop
    case 'h2Error'
        if sys.isBig
            warning('System size might be to large for stopping criterion');
        end
        if isempty(sys.h2Norm), sys.h2Norm = norm(sys,2); end
        if sysr.n>0 %avoid computing when initializing
            stopCrit = norm(sys-sysr,2)/sys.h2Norm;
            stop =  stopCrit <= opts.cure.stopval;
        else
            stopCrit    = NaN;
            stop = 0;
        end
    case 'normROM'
        if length(sysrVec)>2 %run at least three steps
            nNew = norm(sysr); nOld = norm(sysrVec{end-1});
            stopCrit = abs(nNew-nOld)/nNew;
            stop =  stopCrit < opts.cure.stopval;
        else
            stopCrit    = NaN;
            stop        = false;
        end
    case 'nmax'
        stopCrit = sysr.n;
        stop = (stopCrit >=opts.cure.stopval);
    otherwise
        error('The stopping criterion chosen does not exist or is not yet implemented');
end
function [s0,Opts] = initializeShifts(sys,Opts,iCure)
 % Compatibility with 0, 'zero'
 if Opts.cure.initMode == 0, Opts.cure.initMode = 'zero'; end

 % Check if a vector of values has been passed
 if isa(Opts.cure.initMode,'double')
     %==============================
     %  initial shifts were defined
     %==============================
     if length(Opts.cure.initMode) < Opts.cure.nk % not enough
         error('sssMOR:cure:wrongInitialization',...
             'The inital vector of frequencies is incompatible with the desired nk');
     elseif length(Opts.cure.initMode) == Opts.cure.nk % only 4 one iteration
         s0 = Opts.cure.initMode;
         s0 = reshape(s0,[1,length(s0)]); %make sure it's a row vector
     elseif length(Opts.cure.initMode)> Opts.cure.nk % more than one iteration
         firstIndex = (iCure-1)*Opts.cure.nk+1; %select value dep. on iCure
            if firstIndex + Opts.cure.nk -1 > length(Opts.cure.initMode)
                % repeat shifts
                Opts.cure.initMode = [Opts.cure.initMode, Opts.cure.initMode];
            end
            % new set of shifts
            s0 = Opts.cure.initMode(firstIndex:firstIndex+Opts.cure.nk-1);
     end
 else 
     %===================================
     % define initial shifts
     %===================================
     %  choose the number of shifts to compute
     ns0 = Opts.cure.initN;   
     %  compute the shifts
     switch Opts.cure.initMode
         case {'zero','zeros'} %zero initialization
             Opts.cure.initMode = Opts.zeroThres*ones(1,ns0);
         case 'sm' %smallest magnitude eigenvalues
             Opts.cure.initMode = -eigs(sys.a,sys.e,ns0,0, ...
                    struct('tol',1e-6,'v0',sum(sys.b,2))).';
         case 'lm' %largest magnitude eigenvalues
             Opts.cure.initMode = -eigs(sys.a,sys.e,ns0,'lm', ...
                    struct('tol',1e-6,'v0',sum(sys.b,2))).';
         case 'slm' % smallest and largest eigs
             %  decide how many 'lm' and 'sm' eigs to compute
             if ns0 <=4
                 nSm = 2; nLm = 2;
                 Opts.cure.initMode = [];
             else %ns0 > 4
                if isEven(ns0)
                    nSm = ns0/2; nLm = nSm;
                    Opts.cure.initMode = [];
                else
                    nSm = (ns0-1)/2; nLm = nSm;
                    Opts.cure.initMode = 0; %initialize the first shift at 0
                end
                if ~isEven(nSm), nSm = nSm +1; nLm = nLm -1;end
             end
             Opts.cure.initMode = [Opts.cure.initMode,...
                        -eigs(sys.a,sys.e,nSm,'sm', ...
                        struct('tol',1e-6,'v0',sum(sys.b,2))).',...
                        -eigs(sys.a,sys.e,nLm,'lm', ...
                        struct('tol',1e-6,'v0',sum(sys.b,2))).'];
         otherwise
             error('sssMOR:cure:undefinedInitialization',...
                 'The desired initialization for CURE is not defined');
     end
 s0 = Opts.cure.initMode(1:Opts.cure.nk);
end 
%   make sure the initial values for the shifts are complex conjugated
if mod(nnz(imag(s0)),2)~=0 %if there are complex valued shifts...
    % find the s0 which has no compl.conj. partner
    % (assumes that only one shifts has no partner)
    s0(abs(imag(s0)-imag(sum(s0)))<1e-16) = 0;
end
%     cplxpair(s0); %only checking
%   replace 0 with thresh, where threshold is a small number
%   (sometimes the optimizer complaints about cost function @0)
s0(s0==0)=Opts.zeroThres;  
reshape(s0,[1,length(s0)]); %make sure it's a row vector

function writeGif(gifMode)
    filename = 'CURE.gif';
    dt = 1.5;
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    switch gifMode
        case 'create'
            imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',dt);
        case 'append'
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',dt);
        otherwise
            error('Invalid gifMode')
    end   
function isEven = isEven(a)
    isEven = round(a/2)== a/2;
function y = TUM_Blau()
    y = [0 101 189]/255;
function y = TUM_Gruen()
    y = [162 173 0]/255;
function y = TUM_Orange()
    y = [227 114 34]/255;