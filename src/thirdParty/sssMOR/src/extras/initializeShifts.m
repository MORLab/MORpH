function [s0_inp,Rt,s0_out,Lt] = initializeShifts(sys,nShifts,nSets,Opts)
% INITIALIZESHIFTS - Initialize shifts for global and consecutive MOR
% 
% Syntax:
%		[s0_inp,Rt,s0_out,Lt] = INITIALIZESHIFTS(sys,Opts)
% 
% Description:
%       
%       Generates a 1 x nShift vector of shift frequencies s0 according to the 
%       chosen initial shift strategy, which can be defined in Opts.initShiftsStrategy.
%       s0 can be used as inputs to the reduction functions of sssMOR.
%       If necessary, a shift vector for a Krylov output space can be
%       built containing the same shifts as the s0_inp vector.
%
%       Besides, for strategies 'eigs' and 'ROM', right and left tangential
%       directions Rt and Lt are available. For all the other strategies (i.e.
%       'ADI', 'constant', 'linspaced', etc.) tangential directions with ones
%       are computed:
%       Rt = ones(sys.m,nShifts*nSets); Lt = ones(sys.p,nShifts*nSets);
%       
%       If nSets is specified and larger than one, then s0 is a 1 x nSets
%       cell array of shift vectors of size nShifts.
%       
% Input Arguments:
%		*Required Input Arguments:*
%       -sys:           An sss-object containing the LTI system
%       -nShifts:       Number of shifts
%       -nSets:         Number of sets of shifts to be initialized
%
%		*Optional Input Arguments:*
%		-Opts:              A structure containing following fields
%			-.initShiftsStrategy:  	strategy for shift generation;
%                                   [ADI / constant / ROM / {eigs} / 
%                                   linspaced / logspaced / random / lognrnd]
%                                   mixed strategy for real and imag part possible 
%			-.shiftType:            type of shifts;
%                                   [{'conj'} / 'real' / imag]
%			-.wmin:                 lower bound of generated shifts;
%                                   [{|eigs(sys,'sm')|}, positive double]
%			-.wmax:                 upper bound of generated shifts;
%                                   [{|eigs(sys,'lm')|}, positive double]
%           -.shifts.*              refer to opts.shifts in MESS_PARA for 
%                                   more info  
%			-.eigsTyp:              choice of eigenvalues for eigs and ROM;
%                                   [{'sm'} / 'lm' / 'li' / 'si'/ 'lr' / 'sr' / 'la' / 'sa']
%			-.constValue:           value for constant shift strategy;
%                                   [{0} / double]
%			-.offset:               offset for plain imag or real shifts 
%                                   [{0} / double]
%
% Output Arguments:
%       -s0_inp:            Vector (of sets) of shift frequencies for input Krylov subspace 
%       -Rt:                Matrix of right tangential directions
%       -s0_out:            Vector (of sets) of shift frequencies for output Krylov subspace
%                           (the same as for the input space s0_inp)
%       -Lt:                Matrix of left tangential directions
%
% Examples:
%       By default, initializeShifts generates shifts with the eigs
%       strategy. In this example, 10 shifts are generated by taken the
%       mirror eigenvalues with smallest magnitude.
%
%> sys = sss('building');
%> s0 = initializeShifts(sys,10); plot(complex(s0'),'x');
%
%       The behavior of the function can be customized using the
%       option structure Opts, e.g. mixed strategies with logspaced real 
%       and linspaced imag part of the shifts
%> Opts = struct('initShiftsStrategy',{{'logspaced','linspaced'}});
%> s0 = initializeShifts(sys,10,[],Opts); plot(complex(s0'),'x');
%
% See Also: 
%		cure, cirka, spark, crksm, sss/lyapchol
%
% References:
%		[1] Michael Ott (2016)*, Strategien zur Initialisierung der Entwicklungspunkte f�r H2-optimale Modellordnungsreduktion
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
% Authors:      Michael Ott, Siyang Hu, Alessandro Castagnotto
%               Maria Cruz Varona, Paul Heidenreich
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  28 Nov 2017
% Copyright (c) 2017 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------


%% Execution parameters
Def.initShiftsStrategy     = 'eigs';                % initialisation strategy
Def.constValueInp          = 0;                     % constant shift for input space
Def.constValueOut          = 0;                     % constant shift for output space
Def.wmin                   = abs(eigs(sys,1,'sm')); % lower bound  
Def.wmax                   = abs(eigs(sys,1));      % upper bound
Def.eigsType               = 'sm';                  % eigs parameter
Def.shiftTyp               = 'conj';                % plain imaginary shifts
Def.offset                 = 0;                     % global offset for shifts

% create the options structure
if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end

%% Check inputs
if ~exist('nSets','var') || isempty(nSets)
    nSets = 1;
end

% check for valid combination of strategies
if ~iscell(Opts.initShiftsStrategy)
    Opts.initShiftsStrategy = {Opts.initShiftsStrategy};
else
    if (any(strcmp(Opts.initShiftsStrategy{1},{'ADI','eigs','ROM','constant'}))...
            ||any(strcmp(Opts.initShiftsStrategy{end},{'ADI','eigs','ROM','constant'})))...
            && (length(Opts.initShiftsStrategy) > 1 || ~strcmp(Opts.shiftTyp,'conj'))
        error('invalid choice of initial shift strategy');
    end
end

% dealing with odd number Opts.nShifts
oddOrder = 0;

if mod(nShifts,2)
   nShifts = nShifts + 1;
   oddOrder = 1;
end

% initialize tangential directions with ones
Rt = ones(sys.m,nShifts*nSets);
Lt = ones(sys.p,nShifts*nSets);

%% Initialize shifts
switch Opts.initShiftsStrategy{1}
    case 'constant'
        s0_inp = Opts.constValueInp*ones(nSets,nShifts);  
        if nargout > 2
            s0_out = Opts.constValueOut*ones(nSets,nShifts);  
        end
    case 'eigs'
        [Rev,s0_inp] = (eigs(sys,nShifts*nSets,Opts.eigsType));
        s0_inp = -(diag(s0_inp));
        idxUnstable = real(s0_inp)<0;   % mirror shifts if unstable
        s0_inp(idxUnstable) = -s0_inp(idxUnstable);
        try
            cplxpair(s0_inp);
        catch
            s0_inp(end)=real(s0_inp(end));
        end
        
         % get right tangential directions
         if nargout > 1 && sys.isSiso == 0 
             Rt = full((Rev.'*sys.B)).';
             if nargout == 3
                s0_out = s0_inp;
             elseif nargout > 3
                % get s0_out and left tangential directions
                s0_out = s0_inp;
                if sys.issymmetric 
                    Lev = Rev;
                else
                    [Lev,~] = (eigs(sys',nShifts*nSets,Opts.eigsType)); 
                end 
                Lt = full(sys.C*Lev);
             end
         elseif nargout > 1 && sys.isSiso == 1
            s0_out = s0_inp;
         end
        
    case 'ADI'
        Def.shifts.method = 'heur';
        Def.shifts.b0              = ones(sys.n,1);         
        Def.shifts.info            = 0;                    
        Def.lse = 'gauss'; %lse (used only for adi)
        
        if ~exist('Opts','var') || isempty(Opts)
            Opts = Def;
        else
            Opts = parseOpts(Opts,Def);
        end
        
        if strcmp(Opts.shifts.method,'projection') || strcmp(Opts.shifts.method,'wachspress')
            Opts.shifts.method = 'heur';
            warning(['ADI shift method changed from projection/wachspress to heur, ',...
                     'because projection only computes one real or two complex conjugated shifts ',...
                     'and wachspress does not always return the desired number of shifts.'])
        end
        
        if ~sys.isDae
            % options for mess
            % eqn struct: system data
            eqn=struct('A_',sys.A,'E_',sys.E,'B',sys.B,'C',sys.C,'type','N','haveE',1);
            
            % opts struct: mess options            
            messOpts.shifts = Opts.shifts;
            
            % user functions: default
            if strcmp(Opts.lse,'gauss')
                oper = operatormanager('default');
            elseif strcmp(Opts.lse,'luChol')
                if sys.isSym
                    oper = operatormanager('chol');
                else
                    oper = operatormanager('lu');
                end
            end
        end
        
        % get adi shifts
        s0_inp = mess_para(eqn,messOpts,oper);
        
        % mirror shifts if unstable
        idxUnstable = real(s0_inp)<0;
        s0_inp(idxUnstable) = -s0_inp(idxUnstable);
        
        % truncation (of last real shift if length of s0 too long)
        if length(s0_inp) > nSets*nShifts
            reals0 = find((imag(s0_inp)==0));
            if ~isempty(reals0)
                s0_inp(reals0(end))=[];
            else
                error('Truncation does not work, s0 too long')
            end
        end
        
        % define output
        if nargout > 2
            s0_out = s0_inp;
        end
        
    case 'ROM'
        mineig=-eigs(sys,1,'sm'); % use mirrored eigenvalue with smallest magnitude
        
        if ~isreal(mineig)
            multip = ceil(nSets*nShifts/2);
            sysr = rk(sys,[mineig,conj(mineig);multip, multip],[mineig,conj(mineig);multip,multip]);
            s0_inp = single(-eig(sysr).');
            if length(s0_inp) > nSets*nShifts
                idx = find(imag(s0_inp)==0);
                if isempty(idx)
                    s0_inp=sort(s0_inp);
                else
                    [~,idx2] = sort(abs(s0_inp(idx))); %ascend
                    s0_inp(idx(idx2(end))) = [];
                end
            end
            
        else
            % always hermite case
            if size(sys.B,2) == size(sys.C,1)
                sysr = rk(sys,[mineig;nShifts],[mineig;nShifts]);
            else
                [~,V,~] = rk(sys,[mineig;nShifts]);
                [~,~,W] = rk(sys,[],[mineig;nShifts]);
                if size(V,2)<size(W,2)
                    V=[V,W(:,size(V,2)+1:size(W,2))];
                elseif size(V,2)>size(W,2)
                    W=[W,V(:,size(W,2)+1:size(V,2))];
                end
                sysr = projectiveMor(sys,V,W);
            end
        end
        
        [Rev,s0_inp] = (eigs(sysr,nShifts*nSets,Opts.eigsType));
        s0_inp = -(diag(s0_inp));
        idxUnstable = real(s0_inp)<0;   % mirror shifts if unstable
        s0_inp(idxUnstable) = -s0_inp(idxUnstable);
        try
            cplxpair(s0_inp);
        catch
            s0_inp(end)=real(s0_inp(end));
        end
        double(s0_inp);
        
         % get right tangential directions
         if nargout > 1 && sys.isSiso == 0 
             Rt = full((Rev'*sysr.B))';
             if nargout == 3
                 s0_out = s0_inp;
             elseif nargout > 3
                % get s0_out and left tangential directions
                s0_out = s0_inp;
                if sys.issymmetric 
                    Lev = Rev;
                else
                    [Lev,~] = (eigs(sysr',nShifts*nSets,Opts.eigsType)); 
                end 
                Lt = full(sysr.C*Lev);
             end
         elseif nargout > 1 && sys.isSiso == 1
            s0_out = s0_inp;
         end
    
    otherwise % grid and random based strategies
        
        % check for valid combination of strategies
        if length(Opts.initShiftsStrategy)==1 && strcmp(Opts.shiftTyp,'conj')
            Opts.initShiftsStrategy = [Opts.initShiftsStrategy; Opts.initShiftsStrategy];
        elseif ~strcmp(Opts.shiftTyp,'conj') && length(Opts.initShiftsStrategy)>1
            error('invalid choice of initial shift strategy');
        end
        
        % compute grid size and grid parameters
        if strcmp(Opts.shiftTyp,'conj')
            iSplit = ceil(sqrt(nShifts*nSets/2));
            NoS = nShifts*nSets/2;
        elseif strcmp(Opts.shiftTyp,'imag')
            iSplit = nShifts*nSets/2;
            NoS = iSplit;
        else
            iSplit = nShifts*nSets;
            if Opts.offset ~= 0
               iSplit = iSplit/2;
            end
            NoS = iSplit;
        end
        
        for ii = 1:length(Opts.initShiftsStrategy)
            switch Opts.initShiftsStrategy{ii}
                case 'linspaced'
                    s0_temp=linspace(Opts.wmin,Opts.wmax,iSplit);
                    s0_parts{ii} = s0_temp;
                    
                case 'logspaced'
                    s0_temp=logspace(log10(Opts.wmin),log10(Opts.wmax),iSplit);
                    s0_parts{ii} = s0_temp;
                    
                case 'random'
                    s0_temp=Opts.wmin+rand(1,NoS)*(Opts.wmax-Opts.wmin);
                    s0_parts{ii} = s0_temp;
                    
                case 'lognrnd'
                    s0_temp = lognrnd(log(max(Opts.wmin,1)),(log(Opts.wmax)-log(max(Opts.wmin,1)))/2.576,1,NoS);       
                    s0_parts{ii} = s0_temp;
                    
                otherwise
                    error('unknown initial shift strategy');
            end
        end
        
        % generate final s0 matrix
        if ~strcmp(Opts.shiftTyp,'real')
            if strcmp(Opts.shiftTyp,'conj')
                if  ~isempty(strfind(Opts.initShiftsStrategy{1},'spaced'))
                    s0_real = repelem(s0_parts{1},2*iSplit);
                    s0_real = s0_real(1:nShifts*nSets);
                else
                    s0_real = repelem(s0_parts{1},2);
                end
            else
                s0_real = Opts.offset;
            end
            if ~isempty(strfind(Opts.initShiftsStrategy{end},'spaced'))
                s0_imag = repmat(repelem(s0_parts{end},2).*repmat([1,-1],1,iSplit),1,iSplit)*1i;
                s0_imag = s0_imag(1:nShifts*nSets);
            else
                s0_imag = repelem(s0_parts{end},2).*repmat([1,-1],1,NoS)*1i;
            end
            s0_inp = s0_real + s0_imag;
        else
            if Opts.offset ~= 0
                s0_inp = repelem(s0_parts{:},2)+Opts.offset*repmat([1,-1],1,NoS)*1i;
            else
                s0_inp = s0_parts{:};
            end
        end
        
        % define output
        if nargout > 2
            s0_out = s0_inp;
        end
end

% plain real shifts at the end, when the number of shifts is odd
if oddOrder
    s0_inp(end) = [];  s0_inp(end) = real(s0_inp(end));
    if nargout > 1
        Rt(:,end) = [];  Rt(:,end) = real(Rt(:,end));
        if nargout == 3
            s0_out(:,end) = [];  s0_out(:,end) = real(s0_out(:,end));
        elseif nargout > 3
            s0_out(end) = [];  s0_out(end) = real(s0_out(end));
            Lt(:,end) = [];  Lt(:,end) = real(Lt(:,end));
        end
    end
    nShifts = nShifts - 1;
end

% Change s0_inp/s0_out & Rt/Lt from matrix to cell if several sets were computed
if nSets > 1
    s0_inp = reshape(s0_inp,1,nShifts*nSets);        s0_inp = mat2cell(s0_inp,1,nShifts*ones(1,nSets));
    if nargout > 1   
        Rt = reshape(Rt,sys.m,nShifts*nSets);        Rt = mat2cell(Rt,sys.m,nShifts*ones(1,nSets));
        if nargout == 3
            s0_out = reshape(s0_out,1,nShifts*nSets);    s0_out = mat2cell(s0_out,1,nShifts*ones(1,nSets));
        elseif nargout > 3
            s0_out = reshape(s0_out,1,nShifts*nSets);    s0_out = mat2cell(s0_out,1,nShifts*ones(1,nSets));
            Lt = reshape(Lt,sys.p,nShifts*nSets);        Lt = mat2cell(Lt,sys.p,nShifts*ones(1,nSets));
        end
    end
end

end
    