function [varargout] = solveLse(varargin)
% SOLVELSE - Solve linear system of equations
%
% Syntax:
%       SOLVELSE(A)
%       X                                = SOLVELSE(A,B)
%       [X,Y]                            = SOLVELSE(A,B,C)
%       [X,Y,Sx,Rx,Sy,Ly]                = SOLVELSE(A,B,E,s0)
%       [X,Sx,Rx]                        = SOLVELSE(A,B,E,s0,Rt)
%       [X,Y,Sx,Rx,Sy,Lx]                = SOLVELSE(A,B,C,E,s0,Rt,Lt)
%       [X,Y,Sx,Rx,Sy,Lx]                = SOLVELSE(sys)
%       [X,Y,Sx,Rx,Sy,Lx]                = SOLVELSE(sys,s0)
%       [X,Y,Sx,Rx,Sy,Lx]                = SOLVELSE(sys,s0,Rt,Lt)
%       [X,Y,Sx,Rx,Sy,Lx]                = SOLVELSE(A,B,C,E,s0,...,Opts)
%       [X,Sxj,Rxj]                      = SOLVELSE(jCol,X,A,B,E,s0,Rt)
%       [X,Y,Sxj,Rxj,Syj,Lyj]            = SOLVELSE(jCol,X,Y,A,B,C,E,s0,Rt,Lt)
%
% Description:
%       This function solves linear systems of equations X=(A-s0*E)\B
%       corresponding to shifts s0. The order of the shifts is crucial for
%       reusing already computed factorizations, so it is recommended to
%       sort the shifts in advance. If the output matrix C is passed in
%       addition, then SOLVELSE computes the solutions X=(A-s0*E)\B and
%       Y=C/(A-s0*E).'.
%
%       If the matrix E is empty or not specified, X=A\B is computed. If s0
%       is Inf, then the Markov parameter X=(A-s0*E)*B is computed.
%
%       If a system (ss, sss or ssRed) is passed to the function,
%       X=(A-s0*E)\B is computed if shifts s0 are specified, and X=(A-E)\B
%       otherwise.
%
%       If Opts.krylov is set, a Krylov subspace [1-3] is created from the
%       solutions. If C is specified, this means that input and output
%       Krylov subspaces corresponding to the same shifts are computed.
%       Optionally this function computes the Sylvester matrices Sv, Rv, Sy
%       and Ly. For more details, please refer to arnoldi.
%
%       In case of MIMO models, matrices of tangential directions Rt
%       (and Lt) have to be defined. They must have the same number of
%       columns as the shifts s0, so that for each tangential direction it
%       is clear to which shift it belongs.
%
%       If a shifted system has to be solved, sometimes the solution
%       directions need to be modified, e.g. orthogonalized, after every
%       computation. In this case, a vector of shift s0, a current solution
%       V and an index jCol can be passed to SOLVELSE, which will only
%       determine the solution corresponding to jCol and update this
%       column of V. SOLVELSE will take the previous shift and tangential
%       direction into account to reuse lu-factors.
%
%       If a system without shifts (e.g. A\B) has to be solved several
%       times, it is recommended to first call SOLVELSE(A). This will
%       initialize the function by storing the lu-factors of A. If
%       Opts.reuseLU=true is set afterwards, SOLVELSE will check if the
%       lu-factors belong to A and reuse them.
%
%
% Input Arguments:
%       *Required Input Arguments:*
%       -A/B:  System matrices or right/left side matrices
%       -sys:  an ss/sss/ssRed-object containing the LTI system
%       *Optional Input Arguments:*
%       -s0:                Vector of complex conjuate expansion points
%       -E/C:               System matrices
%       -Rt,Lt:             Matrix of right/left tangential directions
%       -X/Y:               Matrix containing lse solutions for jCol-1
%       -Opts:              a structure containing following options
%           -.lse:          use LU or hessenberg decomposition
%                           [{'sparse'} / 'full' / 'hess' / 'iterative' / 'gauss']
%           -.dgksTol:      tolerance for dgks orthogonalization
%                           [{1e-12} / positive float]
%           -.krylov:       standard or cascaded krylov basis
%                           [{0} / 'standardKrylov' / 'cascadedKrylov']
%           -.maxiterlse:   maximum number of iterations in iterSolve
%                           [{1e3} / positive integer]
%           -.tollse:       residual tolerance in iterSolve
%                           [{1e-6} / positive float]
%           -.solver:       preferred solver in iterSolve
%                           [{'cgs'} / 'bicgstab' / 'bicg']
%           -.verbose:      show warnings?
%                           [{1} / 0]
%           -.force:        force solve iteratively when not converging
%                           [{0} / 1]
%           -.refine:       iterative refinement
%                           [{false} / 'wilkinson' / 'cgs']
%           -.refTol:       iterative refinement tolerance
%                           [{1e-8} / positive float]
%           -.refMaxiter:   iterative refinement maximum number of iterations
%                           [{5} / positive integer]
%           -.reuseLU:      reuse previously computed lu-factors
%                           [{0} / 1]
%           -.reuseTol:     tolerance for reusing lu-factors norm(LU-A)<tol
%                           [{1e-10} / positive float]
%
% Output Arguments:
%       -X:        Lse solution corresp. to B/Orthonormal basis spanning the input Krylov subsp.
%       -Sx:       Matrix of input Sylvester Eq.
%       -Rx:       Right tangential directions of Sylvester Eq., (mxq) matrix
%       -Y:        Lse solution corresp. to C/Orthonormal basis spanning the output Krylov subsp.
%       -Sy:       Matrix of output Sylvester Eq.
%       -Ly:       Left tangential directions of Sylvester Eq., (pxq) matrix
%       -Sxj:      Column of input Sylvester Eq. corresponding to jCol
%       -Rxj:      Column of right tangential directions of Sylvester Eq. corresponding to jCol
%       -Syj:      Column of output Sylvester Eq. corresponding to jCol
%       -Ryj:      Column of right tangential directions of Sylvester Eq. corresponding to jCol
%
% Examples:
%       This code computes A^-1*B for two sparse random matrices A and B.
%
%> A=sprand(10,10,0.6); B=sprand(10,1,0.6);
%> x=solveLse(A,B); norm(x-A\B)
%
%       The computation is now repeated for several iterations. To reuse
%       the LU factors, first initialize solveLse and set the option 
%       'reuseLU' to true.
%
%> A=sprand(10,10,0.6); x=sprand(10,1,0.6); Opts.reuseLU=true;
%> solveLse(A,Opts); %initialize solveLse
%> for i=1:round(rand(1)*5)
%>     x=solveLse(A,x,Opts);
%> end
%
%       To solve a shifted system, pass a shift vector s0 to the function:
%
%> sys=sss('building'); s0=rand(8,1);
%> x=solveLse(sys,s0);
%
% See Also:
%       arnoldi, rk, irka, projectiveMor
%
% References:
%       * *[1] Duff et al. (1986)*, Direct methods for sparse matrices
%       * *[2] Saad (2003)*, Iterative methods for sparse linear systems
%       * *[3] Golub, Van Loan (1996)*, Matrix computations
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sss">sss</a>, a Sparse State-Space and System Analysis
% Toolbox developed at the Chair of Automatic Control in collaboration
% with the Professur fuer Thermofluiddynamik, Technische Universitaet Muenchen.
% For updates and further information please visit <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sss">sss</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Heiko Panzer, Alessandro Castagnotto, Maria Cruz Varona,
%               Lisa Jeschek
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  14 Aug 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

Def.lse         = 'sparse'; %use sparse or full LU or lse with Hessenberg decomposition {'sparse', 'full','hess','iterative', 'gauss'}
Def.dgksTol     = 1e-12;    %orthogonality tolerance: norm(V'*V-I,'fro')<tol
Def.krylov      = 0;        %standard or cascaded krylov basis (only for siso) {0,'cascade'}

Def.reuseLU     = false;    %reuse lu (false/true)
Def.reuseTol    = 1e-10;

Def.solver      = 'cgs';%first iterative solver to try
Def.maxiterlse  = 1000; %maximum number of iterations in iterative solver
Def.tollse      = 1e-6; %residual tolerance in iterSolve
Def.verbose     = 1;    %display warnings when iterative methods fail
Def.force       = 0;    %not converging in iterSolve leads to error (0) or warning (1)

Def.refine      = false; % iterative refinement (false, 'wilkinson', 'cgs')
Def.refTol      = 1e-8;
Def.refMaxiter  = 5;

%% parsing of inputs
if isa(varargin{end},'struct')
    Opts=varargin{end};
    varargin=varargin(1:end-1);
end

% create the options structure
if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end

% input of sys
if isa(varargin{1},'sss') || isa(varargin{1},'ss') || isa(varargin{1},'ssRed')
    sys=varargin{1};
    if isa(sys,'ssRed')
        Opts.lse='hess';
    end
    A=sys.A;
    B=sys.B;
    C=sys.C;
    E=sys.E;
    initLse=false;
    if nargout==1 || nargout==3
        hermite=false;
    else
        hermite=true;
    end
    if length(varargin)>1
        s0=varargin{2};
        if length(varargin)==3
            Rt=varargin{3};
            hermite=false;
        elseif length(varargin)==4
            Rt=varargin{3};
            Lt=varargin{4};
            hermite=true;
        end
    end
    
elseif length(varargin)>1
    % input of matrices
    A=varargin{1};
    B=varargin{2};
    initLse=false;
    
    switch length(varargin)
        case 2
            hermite=false;
            % check if Rt directions are necessary
            if size(B,2)>1
                Rt = speye(size(B,2));
                s0=zeros(1,size(B,2));
                withoutE=true;
            end
        case 3
            if ~isscalar(varargin{1}) && size(varargin{3},1)==size(A,1) && size(varargin{3},2)==size(A,2)
                error('Please specify s0.');
            elseif size(varargin{3},1)==size(B,2) && size(varargin{3},2)==size(B,1)
                C=varargin{3};
                hermite=true;
                if size(B,2)>1
                    Rt = speye(size(B,2));
                    Lt = Rt;
                    s0=zeros(1,size(B,2));
                    withoutE=true;
                end
            else
                error('Wrong input.');
            end
        case 4
            if size(varargin{3},1)==size(A,1) && size(varargin{3},2)==size(A,2) && (nargout==1 || nargout==3)
                E=varargin{3};
                s0=varargin{4};
                hermite=false;
            elseif size(varargin{4},1)==size(A,1) && size(varargin{4},2)==size(A,2) && (nargout==2 || nargout >3)
                C=varargin{3};
                E=varargin{4};
                s0=1;
                hermite=true;
            else
                error('Wrong input');
            end
        case 5
            if size(varargin{3},1)==size(A,1) && size(varargin{3},2)==size(A,2) && (nargout==1 || nargout==3)
                E=varargin{3};
                s0=varargin{4};
                Rt=varargin{5};
                hermite=false;
            elseif size(varargin{4},1)==size(A,1) && size(varargin{4},2)==size(A,2) && (nargout==2 || nargout >3)
                C=varargin{3};
                E=varargin{4};
                s0=varargin{5};
                hermite=true;
            else
                error('Wrong input');
            end
        case 6
            if isscalar(varargin{1})
                jCol=varargin{1};
                V=varargin{2};
                A=varargin{3};
                B=varargin{4};
                E=varargin{5};
                s0=varargin{6};
                hermite=false;
            else
                error('Wrong input');
            end
        case 7
            if isscalar(varargin{1}) && (nargout==1  || nargout==3)
                jCol=varargin{1};
                V=varargin{2};
                A=varargin{3};
                B=varargin{4};
                E=varargin{5};
                s0=varargin{6};
                Rt=varargin{7};
                hermite=false;
            elseif nargout==2 || nargout>3
                C=varargin{3};
                E=varargin{4};
                s0=varargin{5};
                Rt=varargin{6};
                Lt=varargin{7};
                hermite=true;
            else
                error('Wrong input.');
            end
        case 8
            if isscalar(varargin{1})
                jCol=varargin{1};
                V=varargin{2};
                W=varargin{3};
                A=varargin{4};
                B=varargin{5};
                C=varargin{6};
                E=varargin{7};
                s0=varargin{8};
                hermite=true;
            else
                error('Wrong input');
            end
        case 10
            if isscalar(varargin{1})
                jCol=varargin{1};
                V=varargin{2};
                W=varargin{3};
                A=varargin{4};
                B=varargin{5};
                C=varargin{6};
                E=varargin{7};
                s0=varargin{8};
                Rt=varargin{9};
                Lt=varargin{10};
                hermite=true;
            else
                error('Wrong input');
            end
        otherwise
            error('Wrong inputs');
    end
elseif length(varargin)==1
    A=varargin{1};
    initLse=true;
    hermite=false;
    withoutE=true;
else
    error('Wrong input.');
end


% check E-matrix, tangential directions and IP
if ~exist('E','var') || isempty(E)
    withoutE=true;
    if ~exist('s0','var') || isempty(s0)
        s0=0;
    end
else
    withoutE=false;
end

if ~exist('Rt','var') && ~initLse
    if size(B,2)==1 %siso
        Rt=ones(1,length(s0));
        if hermite && size(C,1)==1
            Lt=ones(1,length(s0));
        end
    else
        error('Please specify tangential directions.');
    end
end

if exist('jCol','var') && ~isempty(jCol) && strcmp(Opts.lse,'hess')
    error('jCol and Opts.lse=hess are not compatible');
end

% If the 'full' option is selected for LU, convert E,A once to full
if withoutE
    if strcmp(Opts.lse,'full')
        A = full(A);
    elseif strcmp(Opts.lse,'hess')
        [P,A] = hess(full(A)); B = P.'*B; if hermite, C = C*P; end
    elseif strcmp(Opts.lse,'sparse')
        A=sparse(A);
    end
else
    if strcmp(Opts.lse,'full')
        E = full(E); A = full(A);
    elseif strcmp(Opts.lse,'hess')
        [A,E,Q,Z] = hess(full(A),full(E)); B = Q*B; if hermite, C = C*Z; end
    elseif strcmp(Opts.lse,'sparse')
        E = sparse(E); A=sparse(A);
    end
end

if initLse
    if ~strcmp(Opts.lse,'sparse') && ~strcmp(Opts.lse,'full')
        error('Initialization of solveLse requires Opts.lse="sparse"/"full".');
    end
    Opts.krylov='init';
    nextDirection(1,0,zeros(size(A,1),1));
    
elseif exist('jCol','var') && ~isempty(jCol)
    
    if hermite
        [V, SRsylv, Rsylv, W, SLsylv, Lsylv]    = nextDirection(jCol, s0, V, W);
    else
        [V, SRsylv, Rsylv]                      = nextDirection(jCol, s0, V);
    end
    
    % output
    if hermite
        varargout{1}=V;
        varargout{2}=W;
        varargout{3}=SRsylv;
        varargout{4}=Rsylv;
        varargout{5}=SLsylv;
        varargout{6}=Lsylv;
    else
        varargout{1}=V;
        varargout{2}=SRsylv;
        varargout{3}=Rsylv;
    end
else
    % preallocate memory
    q=length(s0)+nnz(imag(s0));
    V=zeros(size(A,1),q);
    Rv=zeros(size(B,2),q);
    Sv=zeros(q);
    if hermite
        W = zeros(size(A,1),q);
        Lw = zeros(size(C,1),q);
        Sw=zeros(q);
    end
    
    for jCol=1:length(s0)
        if hermite
            [V, SRsylv, Rsylv, W, SLsylv, Lsylv]    = nextDirection(jCol, s0, V, W);
        else
            [V, SRsylv, Rsylv]                      = nextDirection(jCol, s0, V);
        end
        Sv(:,jCol) = SRsylv;
        Rv(:,jCol) = Rsylv*Rt(:,jCol);
        if hermite
            Sw(jCol,:) = SLsylv.';
            Lw(:,jCol) = Lsylv*Lt(:,jCol);
        end
    end
    
    if strcmp(Opts.lse,'hess')
        if withoutE
            V=P*V;
            if hermite
                W=P*W;
            end
        else
            V=Z*V;
            if hermite
                W=Q.'*W;
            end
        end
    end
    
    % output
    if hermite
        varargout{1}=V;
        varargout{2}=W;
        varargout{3}=Sv;
        varargout{4}=Rv;
        varargout{5}=Sw;
        varargout{6}=Lw;
    else
        varargout{1}=V;
        varargout{2}=Sv;
        varargout{3}=Rv;
    end
end

    function [V, SRsylv, Rsylv, W, SLsylv, Lsylv] = nextDirection(jCol, s0, V, W)
        %   Get the next direction by solving the lse
        %   Input:  jCol:  Column to be treated
        %           s0:    Vector containing the expansion points
        %           V, W:  solution of lse/Krylov subspaces
        %   Output: V, W:  Updated lse solutions/Krylov subspace
        %           SRsylv: update of column jCol of the Sylvester matrices
        %                  Sv (e.g. SRsylv(:,jCol)=SRsylv)
        %           Rsylv: update of column jCol of the Sylvester matrices
        %                  Rv (Rsylv either eye(size(B,2)) or
        %                  zeros(size(B,2)), e.g. Rsylv(:,jCol)=Rsylv*Rt(:,jCol)
        %           SLsylv: update of column jCol of the Sylvester matrices
        %                  Sw (e.g. SLsylv(:,jCol)=SLsylv)
        %           Lsylv: update of column jCol of the Sylvester matrices
        %                  Lw (Lsylv either eye(size(C,1)) or
        %                  zeros(size(C,1)), e.g. Lsylv(:,jCol)=Lsylv*Lt(:,jCol)
        
        SRsylv=zeros(size(V,2),1);
        if hermite
            SLsylv=zeros(size(W,2),1);
        end
        
        % build Krylov subspace or just solve lse with current s0, Rt and B
        switch Opts.krylov
            case 0
                tempV=B*Rt(:,jCol);
                newlu=1;
                newtan=1;
                SRsylv(jCol)=s0(jCol);
                Rsylv=eye(size(B,2));
                if hermite
                    SLsylv(jCol)=s0(jCol);
                    Lsylv=eye(size(C,1));
                    tempW = C.'*Lt(:,jCol);
                end
                % reuse old factors, but don't build Krylov subspace (tempV ~= E*V(:,jCol-1))
                if jCol>1
                    if s0(jCol)==s0(jCol-1)
                        newlu=0;
                        if Rt(:,jCol) == Rt(:,jCol-1)
                            newtan=0;
                        end
                    end
                end
                if hermite
                    [V(:,jCol), W(:,jCol)] = lse(newlu, newtan, jCol, s0, tempV, tempW);
                else
                    V(:,jCol) = lse(newlu, newtan, jCol, s0, tempV);
                end
            case 'standardKrylov'
                % new basis vector
                tempV=B*Rt(:,jCol); newlu=1; newtan=1;
                SRsylv(jCol)=s0(jCol);
                Rsylv=eye(size(B,2));
                if hermite
                    SLsylv(jCol)=s0(jCol);
                    Lsylv=eye(size(C,1));
                    tempW = C.'*Lt(:,jCol);
                end
                if jCol>1
                    if s0(jCol)==s0(jCol-1)
                        newlu=0;
                        if Rt(:,jCol) == Rt(:,jCol-1)
                            % Higher order moments, for the SISO and MIMO case
                            newtan = 0;
                            tempV = V(:,jCol-1); %overwrite
                            SRsylv(jCol-1)=1;
                            Rsylv=zeros(size(B,2));
                            if hermite
                                SLsylv(jCol-1)=1;
                                Lsylv=zeros(size(C,1));
                                tempW = W(:,jCol-1);
                            end
                        else
                            newtan = 1;
                        end
                    end
                end
                if hermite
                    [V(:,jCol), W(:,jCol)] = lse(newlu, newtan, jCol, s0, tempV, tempW);
                else
                    V(:,jCol) = lse(newlu, newtan, jCol, s0, tempV);
                end
            case 'cascadedKrylov'
                if size(B,2)==1
                    newlu=1; newtan=1;
                    SRsylv(jCol)=s0(jCol);
                    if hermite
                        SLsylv(jCol)=s0(jCol);
                    end
                    if jCol==1
                        tempV=B;
                        Rsylv=1;
                        if hermite
                            tempW=C.';
                            Lsylv=1;
                        end
                    else
                        if s0(jCol)==s0(jCol-1)
                            newlu=0;
                            tempV=V(:,jCol-1);
                            if hermite
                                tempW=W(:,jCol-1);
                            end
                        else
                            tempV=E*V(:,jCol-1);
                            if hermite
                                tempW=E.'*W(:,jCol-1);
                            end
                        end
                        Rsylv=0;
                        SRsylv(jCol-1)=1;
                        if hermite
                            Lsylv=0;
                            SLsylv(jCol-1)=1;
                        end
                    end
                    if hermite
                        [V(:,jCol), W(:,jCol)] = lse(newlu, newtan, jCol, s0, tempV, tempW);
                    else
                        V(:,jCol) = lse(newlu, newtan, jCol, s0, tempV);
                    end
                else
                    error('sss:solveLse:cascadeSiso','A cascaded Krylov basis is only available for SISO systems.');
                end
            case 'init'
                if ~initLse
                    error('Wrong Opts.krylov.');
                end
                lse(1, 1, 1, 0, V);
            otherwise
                error('Opts.krylov is invalid.');
        end
    end

    function [tempV, tempW] = lse(newlu, newtan, jCol, s0, tempV, tempW)
        %   Solve linear system of equations to obtain the new direction
        %   Moment matching:  newlu=0: tempV=(A-s0_old*E)^-1*(E*tempV)
        %                     newlu=1: tempV=(A-s0*E)^-1*tempV
        %   Markov parameter: newlu=0: tempV=E^-1*(A*tempV)
        %                     newlu=1: tempW=E^-1*tempV
        %   Input:  newlu: new lu decomposition required
        %           newtan: new tangential direction required
        %           jCol: Column to be treated
        %           s0: Vector containing the expansion points
        %           tempV, tempW: previous direction
        %   Output: tempV, tempW: new direction
        persistent R S L U a o;
        
        if isinf(s0(jCol)) %Realization problem (match Markov parameters)
            % Note: s0=[] for initLse, reuseLU, withoutE -> no isinf(s0) possible
            
            % Krylov subspace
            if (newlu==0 && (size(B,2)==1 || newtan==0) && strcmp(Opts.krylov,'standardKrylov')) || strcmp(Opts.krylov,'cascadedKrylov')
                tempV=A*tempV;
                if hermite
                    tempW=A*tempW;
                end
            end
            
            % check if lu-factors are not empty
            if newlu==0 && (isempty(R) && isempty(S) && isempty(L) && isempty(U) && isempty(a) && isempty(o))
                newlu=1;
            end
            
            % compute new lu
            if newlu==1
                try
                    % compute Cholesky factors of E
                    U=[];
                    [R,p,S] = chol(E);
                    if p~=0 % different factorization?
                        error('solveLse:cholp','chol: p~=0');
                    end
                catch err
                    if strcmp(err.identifier,'MATLAB:posdef') || strcmp(err.identifier,'solveLse:cholp')|| ~strcmp(Opts.lse,'sparse')
                        % E is not pos. def -> use LU instead
                        switch Opts.lse
                            case 'sparse'
                                [L,U,a,o,S]=lu(E,'vector');
                            case 'full'
                                [L,U]=lu(E);
                        end
                    else
                        rethrow(err);
                    end
                end
            end
            
            % solve lse
            if ~isempty(U) || strcmp(Opts.lse,'hess') || strcmp(Opts.lse,'gauss') % not cholesky
                switch Opts.lse
                    case 'sparse'
                        oldV=tempV;
                        tempV(o,:) = U\(L\(S(:,a)\tempV)); %LU x(o,:) = S(:,a)\b
                        if hermite
                            oldW=tempW;
                            tempW(o,:) = U\(L\(S(:,a)\tempW));
                        end
                        % iterative refinement
                        if Opts.refine
                            if hermite
                                [tempV,tempW] = iterativeRefinement(E, oldV, tempV, oldW, tempW);
                            else
                                tempV = iterativeRefinement(E,oldV, tempV);
                            end
                        end
                    case 'full'
                        tempV = U\(L\tempV);
                        if hermite
                            tempW = U\(L\tempW);
                        end
                    case 'hess'
                        tempV = E\tempV;
                        if hermite
                            tempW = E\tempW;
                        end
                    case 'gauss'
                        tempV = E\tempV;
                        if hermite
                            tempW = E\tempW;
                        end
                    otherwise
                        error('Lse method not implemented.');
                end
            else % cholesky
                tempV = S*(R\(R.'\(S.'*tempV)));
                if hermite
                    tempW = S*(R\(R.'\(S.'*tempW)));
                end
            end
            
        else %Rational Krylov
            if ~strcmp(Opts.lse,'iterative') %direct methods
                
                % Krylov subspace
                if newlu==0
                    if strcmp(Opts.krylov,'standardKrylov') || strcmp(Opts.krylov,'cascadedKrylov')
                        if size(B,2)==1 %SISO
                            tempV=E*tempV;
                            if hermite, tempW = E.'*tempW; end
                        elseif newtan==0
                            % Tangential matching of higher order moments
                            tempV=E*tempV;
                            if hermite, tempW = E.'*tempW; end
                        end
                    end
                end
                
                % check if lu-factors are not empty
                if newlu==0 && (isempty(R) && isempty(S) && isempty(L) && isempty(U) && isempty(a) && isempty(o))
                    newlu=1;
                end
                
                % check if old lu can be reused (only for full/sparse LU)
                if newlu==1 && Opts.reuseLU && withoutE && ~initLse
                    if strcmp(Opts.lse,'full') && ~isempty(L) && ~isempty(U) && norm(L*U - A,'fro')<Opts.reuseTol
                        newlu=0;
                    elseif (strcmp(Opts.lse,'sparse') && ~isempty(L) && ~isempty(U) ...
                            && ~isempty(a) && ~isempty(o) && ~isempty(S) && norm(L*U - S(:,a)\A(:,o),'fro')<Opts.reuseTol)
                        newlu=0;
                    end
                end
                
                if newlu==1
                    if withoutE
                        switch Opts.lse
                            case 'sparse'
                                % vector LU for sparse matrices
                                [L,U,a,o,S]=lu(A,'vector');
                            case 'full'
                                [L,U] = lu(A);
                        end
                    else
                        switch Opts.lse
                            case 'sparse'
                                % vector LU for sparse matrices
                                [L,U,a,o,S]=lu(A-s0(jCol)*E,'vector');
                            case 'full'
                                [L,U] = lu(A-s0(jCol)*E);
                        end
                    end
                end
                
                if ~initLse
                    % Solve the linear system of equations
                    switch Opts.lse
                        case 'sparse'
                            oldV=tempV;
                            tempV(o,:) = U\(L\(S(:,a)\tempV)); %LU x(o,:) = S(:,a)\b
                            if hermite
                                oldW=tempW;
                                tempW = (S(:,a)).'\(L.'\(U.'\(tempW(o,:))));
                            end %U'L'S(:,a) x = c'(o,:)
                            % iterative refinement
                            if Opts.refine
                                if withoutE
                                    if hermite
                                        [tempV,tempW] = iterativeRefinement(A, oldV, tempV, oldW, tempW);
                                    else
                                        tempV = iterativeRefinement(A,oldV, tempV);
                                    end
                                else
                                    if hermite
                                        [tempV,tempW] = iterativeRefinement(A-s0(jCol)*E, oldV, tempV, oldW, tempW);
                                    else
                                        tempV = iterativeRefinement(A-s0(jCol)*E,oldV, tempV);
                                    end
                                end
                            end
                        case 'full'
                            tempV = U\(L\tempV);
                            if hermite, tempW = (L.'\(U.'\(tempW))); end
                        case 'hess'
                            if withoutE
                                tempV = A\tempV;
                            else
                                tempV = (A-s0(jCol)*E)\tempV;
                            end
                            if hermite
                                if withoutE
                                    tempW = A.'\tempW;
                                else
                                    tempW = (A-s0(jCol)*E).'\tempW;
                                end
                            end
                        case 'gauss'
                            if withoutE
                                tempV = A\tempV;
                            else
                                tempV = (A-s0(jCol)*E)\tempV;
                            end
                            if hermite
                                if withoutE
                                    tempW = A.'\tempW;
                                else
                                    tempW = (A-s0(jCol)*E).'\tempW;
                                end
                            end
                        otherwise
                            error('Lse method not implemented.');
                    end
                end
                
            else %iterative methods
                if withoutE
                    [tempV,flag,method] = iterSolve(A,tempV,newlu,hermite);
                    if Opts.verbose, disp(method); end
                    if hermite
                        [tempW,flag,method] = iterSolve((A).',tempW,0,hermite);
                        if Opts.verbose, disp('entering hermite section'); end
                    end
                else
                    if newlu
                        [tempV,flag,method] = iterSolve(A-s0(jCol)*E,tempV,newlu,hermite);
                        if Opts.verbose, disp(method); end
                        if hermite
                            [tempW,flag,method] = iterSolve((A-s0(jCol)*E).',tempW,0,hermite);
                            if Opts.verbose, disp('entering hermite section'); end
                        end
                    elseif strcmp(Opts.krylov,'standardKrylov') || strcmp(Opts.krylov,'cascadedKrylov')
                        [tempV,flag,method] = iterSolve(A-s0(jCol)*E,E*tempV,newlu,hermite);
                        if Opts.verbose, disp(method); end
                        if hermite
                            [tempW,flag,method] = iterSolve((A-s0(jCol)*E).',E.'*tempW,0,hermite);
                            if Opts.verbose, disp('entering hermite section'); end
                        end
                    else % Opts.krylov==0 and newlu==0 -> tempV instead of E*tempV -> no krylov subspace
                        [tempV,flag,method] = iterSolve(A-s0(jCol)*E,tempV,newlu,hermite);
                        if Opts.verbose, disp(method); end
                        if hermite
                            [tempW,flag,method] = iterSolve((A-s0(jCol)*E).',tempW,0,hermite);
                            if Opts.verbose, disp('entering hermite section'); end
                        end
                    end
                end
            end
        end
        
        if jCol==length(s0) && ~Opts.reuseLU && ~initLse
            clear R S L U a o first sym pd flag method failed nolu solver
        end
        
        function [tempV, tempW, rNormVecV, rNormVecW] = iterativeRefinement(AsE, oldV, tempV, oldW, tempW)
            %   ITERATIVEREFINEMENT - Use iterative LSE solver to refine LSE solution
            %
            %   This function takes as input the approximate solution tempV to the
            %   linear system AsE*tempV = oldV and tries to improve its accuracy. oldV is also
            %   allowed to be a matrix, so that tempV is the solution of a set of lse
            %   with common AsE matrix.
            %
            %   This can be done either by direct methods (compare [1,2]) or
            %   by feeding the data to an iterative solver (so far, only cgs is
            %   implemented). In both cases, the algorihm stops when the desired
            %   accuracy, measured as the relative norm of the residual
            %   r = oldV-AsE*tempV is achieved.
            %
            %   This can be useful when subspecting that a solution X obtained
            %   through direct methods (e.g. lu or chol) may be inaccurate due to
            %   roundoff errors.
            %   Input:    AsE        A-s0(jCol)*E
            %             oldV       tempV before solving lse
            %             tempV      current lse solution
            %   Output:   tempV      improved lse solution
            if Opts.refine
                % Determine if the condition number of the problem is low enough for
                % refinement to make sense
                kU = log10(condest(U)); %it seems that cond(U) is close to cond(AsE)
                if log10(Opts.refTol) < log10(eps) + kU
                    warning(['The condition number of the problem is too high for ',...
                             'iterative refinement to make sense (at least with given tolerance). ',...
                             'Switching Opts.refine to 0.'])
                    Opts.refine = false;
                else   
                    switch Opts.refine
                        case 'wilkinson'
                    oldVNorm = norm(oldV,'fro');
                    k = 0; rNormVecV = zeros(1,Opts.refMaxiter);
                    
                    resV = oldV-AsE*tempV; %residual
                    rNormV = norm(resV,'fro')/oldVNorm;
                    
                    if hermite
                        oldWNorm = norm(oldW,'fro');
                        rNormVecW = zeros(1,Opts.refMaxiter);
                        
                        resW = oldW-AsE.'*tempW; %residual
                        rNormW = norm(resW,'fro')/oldWNorm;
                    end
                    
                    %   Refine!
                    while (rNormV > Opts.refTol || (hermite && rNormW > Opts.refTol)) && k <= Opts.refMaxiter
                        k = k+1;
                        
                        %             D = o*(U\(L\(a*(S\R)))); %lu(AsE)
                        Dv = zeros(size(resV));
                        Dv(o,:) = U\(L\(S(:,a)\resV));
                        tempV = tempV + Dv;
                        
                        resV = oldV-AsE*tempV;
                        rNormV = norm(resV,'fro')/oldVNorm;
                        rNormVecV(k) = rNormV;
                        
                        if hermite
                            Dw = (S(:,a)).'\(L.'\(U.'\(resW(o,:))));
                            tempW = tempW + Dw;
                            resW = oldW-AsE.'*tempW;
                            rNormW = norm(resW,'fro')/oldWNorm;
                            rNormVecW(k) = rNormW;
                        end
                    end
                    
                    if k <= Opts.refMaxiter
                        rNormVecV(k+1:end) = [];
                        if hermite
                            rNormVecW(k+1:end) = [];
                        end
                    end
                    
                        case 'cgs'
                    warning('off','MATLAB:cgs:tooSmallTolerance');
                    for iCol = 1:size(tempV,2)
                        %Solve LSE with desired accuracy
                        [tempv,exitFlagv,rNormVecV] = cgs(AsE,oldV(:,iCol),Opts.refTol,Opts.refMaxiter,L,U,tempV(:,iCol));
                        if hermite
                            [tempw,exitFlagw,rNormVecW] = cgs(AsE.',oldW(:,iCol),Opts.refTol,Opts.refMaxiter,U.',L.',tempW(:,iCol));
                        end
                        %replace column if refinement worked
                        if exitFlagv == 0,
                            tempV(:,iCol) = tempv;
                            if hermite && exitFlagw == 0
                                tempW(:,iCol) = tempw;
                            else
                                warning('sss:iterativeRefinement:cgsNoConvergence','iterative refinement did not work.');
                            end
                        else
                            warning('sss:iterativeRefinement:cgsNoConvergence','iterative refinement did not work.');
                        end
                    end
                    warning('on','MATLAB:cgs:tooSmallTolerance');
                        otherwise
                    error('Iterative refinement method not implemented.');
                    end
                end
            end
        end
    end

    function [ x, varargout ] = iterSolve( A, b, newlu, hermite)
        %ITERSOLVE function to solve the linear system A*x=b iteratively
        %   Detailed explanation goes here
        
        persistent first L U sym pd flag method failed nolu solver;
        
        if isempty(first);
            first = 1;
        end
        
        if first == 1           % set persistent variables on first call of function
            
            failed = {};
            
            if hermite==1 && newlu==0;             % LU factorization for transposed system: A~L*U => A'~U'*L'
                if ~isempty(L) && ~isempty(U)      % ilu, not ichol is already computed
                    temp = L;
                    L = U.';
                    U = temp.';
                    clear temp;
                end
                
                if sym==1 && pd==1
                    A = -A;
                    b = -b;
                end
                
            end
            
            
        end
        
        if newlu==1 && first==1            % new L and U required (different s0)
            
            nolu = 0;
            L = [];
            U = [];
            
            if isempty(solver)
                solver = Opts.solver;
            end
            
            % check for symmetry and definiteness
            
            if norm(A-A','fro')/norm(A,'fro') < 1e-10
                sym = 1;
            else
                sym = 0;
            end
            
            if ispd(-A)
                A = -A;
                b = -b;
                pd = 1;
            else
                pd = 0;
            end
        end
        
        first = 0;                   % first round over
        
        solver = selectSolver(solver, failed);
        
        %     analyse = ['maxiterlse: ',num2str(Opts.maxiterlse),', tol: ',num2str(Opts.tollse)];
        %     disp(analyse);
        
        if sym && pd      % pcg
            if nolu~=1
                if isempty(L) && isempty(U)
                    try
                        L = ichol(A);
                        [x,flag] = pcg(A,b,Def.tollse,Opts.maxiterlse,L,L.',b);
                        method = 'pcg with ichol';
                    catch
                        try
                            [L,U] = ilu(A);
                            [x,flag] = pcg(A,b,Opts.tollse,Opts.maxiterlse,L,U,b);
                            method = 'pcg with ilu';
                        catch
                            [x,flag] = pcg(A,b,Opts.tollse,Opts.maxiterlse,[],[],b);
                            method = 'cg (no preconditioner)';
                            nolu = 1;
                        end
                    end
                elseif isempty(U)   % ichol already available
                    [x,flag] = pcg(A,b,Def.tollse,Opts.maxiterlse,L,L.',b);
                    method = 'pcg with ichol';
                else                % ilu already available
                    [x,flag] = pcg(A,b,Opts.tollse,Opts.maxiterlse,L,U,b);
                    method = 'pcg with ilu';
                end
            else
                [x,flag] = pcg(A,b,Opts.tollse,Opts.maxiterlse,[],[],b);
                method = 'cg (no preconditioner)';
            end
        else
            switch solver
                case 'cgs'
                    if nolu ~= 1                          % LU factorization not yet failed
                        if isempty(L) && isempty(U)        % LU factorization not yet computed
                            try
                                [L,U] = ilu(A);
                                method = 'cgs with ilu';
                            catch
                                method = 'cgs (no preconditioner)';
                                nolu = 1;
                            end
                        else
                            method = 'cgs with ilu';
                        end
                    else
                        method = 'cgs (no preconditioner)';
                    end
                    
                    [x,flag] = cgs(A,b,Opts.tollse,Opts.maxiterlse,L,U,b);
                    
                    if flag ~= 0
                        failed{length(failed)+1} = 'cgs';
                        
                        if Opts.verbose == 1              % print failure
                            msg = ['cgs did not converge.'];
                            warning(msg);
                        end
                    end
                    
                case 'bicgstab'
                    if nolu ~= 1                          % LU factorization not yet failed
                        if isempty(L) && isempty(U)        % LU factorization not yet computed
                            try
                                [L,U] = ilu(A);
                                method = 'bicgstab with ilu';
                            catch
                                method = 'bicgstab (no preconditioner)';
                                nolu = 1;
                            end
                        else
                            method = 'bicgstab with ilu';
                        end
                    else
                        method = 'bicgstab (no preconditioner)';
                    end
                    
                    [x,flag] = bicgstab(A,b,Opts.tollse,ceil(Opts.maxiterlse/2),L,U,b);
                    
                    if flag ~= 0
                        failed{length(failed)+1} = 'bicgstab';
                        if Opts.verbose == 1              % print failure
                            msg = ['bicgstab did not converge.'];
                            warning(msg);
                        end
                    end
                    
                case 'bicg'
                    if nolu ~= 1                           % LU factorization not yet failed
                        if isempty(L) && isempty(U)        % LU factorization not yet computed
                            try
                                [L,U] = ilu(A);
                                method = 'bicg with ilu';
                            catch
                                method = 'bicg (no preconditioner)';
                                nolu = 1;
                            end
                        else
                            method = 'bicg with ilu';
                        end
                    else
                        method = 'bicg (no preconditioner)';
                    end
                    
                    [x,flag] = bicg(A,b,Opts.tollse,ceil(Opts.maxiterlse/2),L,U,b);
                    
                    if flag ~= 0
                        failed{length(failed)+1} = 'bicg';
                        if Opts.verbose == 1              % print failure
                            msg = ['bicg did not converge.'];
                            warning(msg);
                        end
                        
                    end
            end
            
            if flag ~= 0 && length(failed) < 3
                iterSolve(A,b,newlu,hermite);
            end
            
        end
        
        first = [];
        
        if(nargout == 2)
            varargout{1} = flag;
        end
        
        if(nargout == 3)
            varargout{1} = flag;
            varargout{2} = method;
        end
        
        if flag ~= 0
            msg = ['Could not achieve desired tolerance (',num2str(Opts.tollse),') within ',num2str(Opts.maxiterlse),' iterations. Results may be inaccurate!'];
            if Opts.force == 1;
                warning(msg);
            else
                error(msg);
            end
            solver = Opts.solver;
        end
        
    end

    function [solver] = selectSolver( solver, failed )
        solvers = {'cgs','bicgstab','bicg'};
        switch length(failed)
            case 0
            case 1
                index = strcmp(failed, solvers);
                solver = solvers{find(index==0,1)};
                
            case 2
                index1 = strcmp(failed{1}, solvers);
                index2 = strcmp(failed{2}, solvers);
                solver = solvers{index1==index2};
        end
        
    end

end