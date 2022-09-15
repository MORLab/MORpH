function [redSys, V, W, nLU] = arnoldiPH(sys, s0, varargin)
% ARNOLDIPH - reduces the system sys by the rational Arnoldi algorithm with
%           interpolation points specified in s0.
%
% Syntax:
%   redSys = ARNOLDIPH(sys, s0)
%   redSys = ARNOLDIPH(sys, s0, Opts)
%   redSys = ARNOLDIPH(sys, s0, b)
%   redSys = ARNOLDIPH(sys, s0, b, Opts)
%   [redSys,V,W,nLU] = ARNOLDIPH(sys, ...)
%
% Description:
%       redSys = arnoldiPH(sys, s0) returns the reduced system redSys. It
%       computes the ROM by interpolation of the transfer function with
%       Krylov subspaces. The Arnoldi method has been adapted to PH systems
%       (see References).
%       The reduced system order will be as high as the number of
%       interpolation points in s0 for real shifts. Imaginary shifts can be
%       passed as complex conjugate values - if they are unique, the
%       algorithm will complement the complex conjugate shifts.
%
%       redSys = arnoldiPH(sys, s0) can be used for MIMO systems as well. It
%       then computes the block-Krylov subspace for projective model order
%       reduction.
%
%       redSys = arnoldiPH(sys, s0, b) is used for MIMO systems. The matrix b
%       contains the tangent directions for the Krylov subspace as columns.
%       There must be as many tangent directions as there are interpolation
%       points in s0.
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:		phs object
%       - s0:       vector of interpolation points
%       *Optional Input Arguments:*
%       - b:		matrix containing the tangent directions as columns
%       - Opts:     structure with execution parameters
%           - .structurePreservation:   determines how structure is
%                       preserved by the algorithm:
%                       Available options (see also structurePreservation.m):
%                       > 'scaling': use scaled energy coordinates and W = V; [1]
%                       > 'specialInverse': W = Q*V*inv(V'*Q*V); [4]
%                       > 'Cholesky': V'*Q*V = L'*L, V <- V*inv(L), W <- QV; [3]
%                       > 'Cholesky+': V'*Q'*E*V = L'*L, V <- V*inv(L), W <- QV
%                       > 'QV':       V <- V, W <- QV [5]
%                       Default: 'Cholesky'
%           - .constAlg:    use this option to specify how the projection
%                       subspace matrix V is determined 
%                       [{'phs'}, 'sss']
%                       > 'phs': use implementation from MORpH toolbox
%                       > 'sss': use implementation from sssMOR toolbox
%           - .orth:    orthogonalization method that is used for V:
%                       [{'2mgs'} / false / 'mgs']
%                       > false: no orthogonalization
%                       > 'mgs': single Gram-Schmidt
%                       > '2mgs': double application of Gram-Schmidt
%                       See [5] for justification.
%           - .reorthog:  Reorthogonalization of V:
%                         see structurePreservation.m
%           - .verbose: > true: Output messages and warnings during
%                       computation will be shown
%                       > false: Output messages and warnings during
%                       computation will not be shown
%                       [{true} / false]
%           - .arnoldi.*:   other options that will be passed on to the
%                       used construction algorithm (sssMOR->arnoldi);
%                       Please refer to documentation of the respective
%                       algorithm.
%           - .phs.*    Options that will be passed to the phs class
%
%
% Output Arguments:
%       - redSys: 	phsRed object of the reduced system
%       - V:        matrix that was used to obtain the reduced order model
%       - W:        matrix that was used to obtain the reduced order model
%                   (depends on structure preservation option)
%       - nLU:      number of LU decompositions
%
%
% See Also:
%       structurePreservation, phs, arnoldi, demo_arnoldiPH
%
% References:
%       [1] Polyuga (2010), "Model Reduction of Port-Hamiltonian Systems",
%           Dissertation, 2010.
%       [2] S. Gugercin, R. Polyuga, C. Beattie, and A. van der
%           Schaft, "Interpolation-based H2 model reduction for port-
%           Hamiltonian systems," in Proceedings of the 48h IEEE Conference 
%           on Decision and Control (CDC), 2009, pp. 5362–5369.
%       [3] S. Gugercin, R. V. Polyuga, C. Beattie, and A. van der Schaft,
%           "Structure-preserving tangential interpolation for model reduction 
%           of port-Hamiltonian systems," Automatica, 48(9), pp. 1963–1974, 2012.
%       [4] T. Wolf, B. Lohmann, R. Eid, and P. Kotyczka. Passivity and 
%           structure preserving order reduction of linear port-Hamiltonian 
%           systems using Krylov subspaces. Eur. J. Control, 16(4):401–406, 2010.
%       [5] Giraud et al., "The loss of orthogonality in the Gram-Schmidt
%           orthogonalization process", Computers & Mathematics with
%           applications, 50(7), pp. 1069-1075, 2005.
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

%% Parse inputs
narginchk(2,4)
[s0, b, Opts] = parseInputs(sys, s0, varargin{:});

%% Initialize system matrices for certain structure preservation options
if strcmp(Opts.structurePreservation, 'scaling')
    sys = scaling(sys);
end
if strcmp(Opts.structurePreservation, 'Cholesky')
    if sys.isImplicit
        % check symmetry of Q
        if ~(issymmetric(sys.Q)) && norm(sys.Q-sys.Q')/norm(sys.Q) > 1e-12
            warning('MORpH:arnoldiPH:Cholesky_QNotSymmetric',...
                ['Q is not symmetric. The structure preservation strategy '...
                '''Cholesky'' might fail. This function will proceed anyway...\n'...
                'Consider chosing ''Cholesky+'' instead of ''Cholesky''']);
        end
    end
end

%% Determine Krylov subspace matrix V
if isequal(Opts.constAlg, 'sss')
    % sss arnoldi

    % sorting shifts for sssMOR->arnoldi
    [s0_sort,b_sort] = sortShiftsForSSSMOR(s0,b);

    if ~sys.isMIMO || isempty(b)
        [V,~,~,~,~,~,nLU] = arnoldi(sys.E,(sys.J-sys.R)*sys.Q,sys.G-sys.P,s0_sort,Opts);
    else
        [V,~,~,~,~,~,nLU] = arnoldi(sys.E,(sys.J-sys.R)*sys.Q,sys.G-sys.P,s0_sort,b_sort,Opts);
    end


elseif isequal(Opts.constAlg, 'phs')
    % phs version

    if strcmp(Opts.orth,'dgks')
        error("MORpH:arnoldiPH:bad_option_combination",...
            "Construction algorithm 'phs' not compatible with orthogonalization option 'dgks'.");
    end

    if ~sys.isMIMO
        [V,nLU] = computeV((sys.J-sys.R)*sys.Q, sys.G-sys.P, sys.E, s0, Opts, 'SISO');
    else
        [V,nLU] = computeV((sys.J-sys.R)*sys.Q, sys.G-sys.P, sys.E, s0, Opts, 'MIMO', b);
    end

else
    error('MORpH:arnoldiPH:unknownAlgorithm','The specified algorithm in Opts.constAlg is not supported')
end

%% Apply structure preservation and create reduced order model
[J_red, R_red, Q_red, G_red, E_red, P_red, S_red, N_red, V, W] = structurePreservation(sys, V, Opts);

% create phsRed object
redSys = phsRed(J_red, R_red, Q_red, G_red, E_red, P_red, S_red, N_red, Opts.phs);
redSys.method = @arnoldiPH;
redSys.parameters = Opts;
redSys.parameters.shifts = s0;
if sys.isMIMO
    redSys.parameters.tangentDirections = b;
end

end

%% ------------------------ AUXILIARY FUNCTIONS -------------------------

function [V,nLU] = computeV(A,B,E,s0,Opts,flag_SISO_MIMO, b)
% Computes the Krylov-subspace matrix V for matrices A,B,E and shifts s0.
% Additional options are passed by Opts. flag_SISO_MIMO can be 'SISO' or
% 'MIMO'. b (tangent directions matrix) may be provided optionally.

narginchk(6,7);
if nargin == 6
    b = [];
end

% Sort shifts and tangent directions
if strcmp(flag_SISO_MIMO,'SISO')
    % SISO system
    [s0, ~, n, redOrder] = processShifts(s0);
elseif isempty(b)
    % MIMO system without tangent directions
    [s0, ~, n, redOrder] = processShifts(s0, b);
    n_input = size(B,2);
else
    % MIMO system with tangent directions
    [s0, b, n, redOrder] = processShifts(s0, b);
end

% Set nLU (number of LU decompositions "needed" --> LU decomposition
% is actually only performed if shift appears more than once)
nLU = length(s0);

% Preallocate memory for V
V = zeros(size(A,1),redOrder);

% Transform A, B, E to sparse matrices
if strcmp(Opts.lse,'sparse')
    A = sparse(A);
    B = sparse(B);
    E = sparse(E);
end

% Calculate Krylov subspace drections for each shift

% Initialize indices
indexS = 1; % index for sigmas (used for relating sigma with its corresponding n)
indexV = 0; % index for column in V (keeps track of the column which was updated last)

% Iterate through all shifts and determine the related subspace
% directions
for sigma = s0
    % for all shifts sigma in s0

    % Initialize v_old (tangent direction)
    if strcmp(flag_SISO_MIMO,'SISO')
        % SISO system
        v_old = B;
    elseif isempty(b)
        % MIMO system without tangent directions
        v_old = B;
    else
        % MIMO system with tangent directions
        v_old = B*b(:,indexS);
    end

    % Perform LU decomposition if shift will be used more than once
    if n(indexS) > 1 && sigma ~= Inf
        % compute LU decomposition
        switch Opts.lse
            case 'full'
                [L,U] = lu(A-sigma*E);
                p = []; q = []; D = [];
            case 'sparse'
                % inv(D(:,p))*A(:,q) = L*U
                % A*x = b <=> x(q) = inv(U)*inv(L)*inv(D(:,p))*b
                % see documentation of lu for more information
                [L, U, p, q, D] = lu(A - sigma*E, 'vector');
            otherwise
                error('Opts.lse is not well defined. (Can be ''full'' or ''sparse''.)');
        end
    else
        L = []; U = []; p = []; q = []; D = [];
    end

    for j = 1:n(indexS)
        % for each occurence of sigma

        % Real sigma
        if isreal(sigma) && isreal(v_old)
            % Calculate new subspace direction with sigma
            v = nextSubspaceDirection(A,E,sigma,v_old,j,L,U,p,q,D);
            % Find orthogonal vector for extension of V
            switch Opts.orth
                case false

                case 'mgs'
                    v = gramSchmidt(v,V(:,1:indexV));
                case '2mgs'
                    v = gramSchmidt(v,V(:,1:indexV)); v = gramSchmidt(v,V(:,1:indexV));
            end

            % Extend V
            if strcmp(flag_SISO_MIMO,'SISO')
                % SISO system
                V(:, indexV+1) = v;
                indexV = indexV+1;
            elseif isempty(b)
                % MIMO system without tangent directions
                V(:, indexV+1:indexV+n_input) = v;
                indexV = indexV+n_input;
            else
                % MIMO system with tangent directions
                V(:, indexV+1) = v;
                indexV = indexV+1;
            end
            v_old = v;

            % Complex sigma
        else
            % Instead of complex directions use complex conjugate vector (also in V) to obtain 2 real vectors!
            % Calculate new subspace direction with sigma
            v = nextSubspaceDirection(A,E,sigma,v_old,j,L,U,p,q,D);   % complex!
            v_old = v;  % v_old needs to be assigned here since direction v changes by splitting and orthogonalizing the complex vector

            % Real part
            switch Opts.orth
                case false
                    v_real = real(v);
                case 'mgs'
                    v_real = gramSchmidt(real(v),V(:,1:indexV));
                case '2mgs'
                    v_real = gramSchmidt(real(v),V(:,1:indexV)); v_real = gramSchmidt(v_real,V(:,1:indexV));
            end
            % Extend V
            if strcmp(flag_SISO_MIMO,'SISO')
                % SISO system
                V(:, indexV+1) = v_real;
                indexV = indexV+1;
            elseif isempty(b)
                % MIMO system without tangent directions
                V(:, indexV+1:indexV+n_input) = v_real;
                indexV = indexV+n_input;
            else
                % MIMO system with tangent directions
                V(:, indexV+1) = v_real;
                indexV = indexV+1;
            end

            % Imaginary part
            switch Opts.orth
                case false
                    v_imag = imag(v);
                case 'mgs'
                    v_imag = gramSchmidt(imag(v),V(:,1:indexV));
                case '2mgs'
                    v_imag = gramSchmidt(imag(v),V(:,1:indexV)); v_imag = gramSchmidt(v_imag,V(:,1:indexV));
            end
            % Extend V
            if strcmp(flag_SISO_MIMO,'SISO')
                % SISO system
                V(:, indexV+1) = v_imag;
                indexV = indexV+1;
            elseif isempty(b)
                % MIMO system without tangent directions
                V(:, indexV+1:indexV+n_input) = v_imag;
                indexV = indexV+n_input;
            else
                % MIMO system with tangent directions
                V(:, indexV+1) = v_imag;
                indexV = indexV+1;
            end

        end %if isreal(sigma)

    end % for j = 1:n(indexS)

    % Update indexS
    indexS = indexS + 1;

end % --> next shift
end % of function computeV

function v = nextSubspaceDirection(A,E,sigma,v,j,L,U,p,q,D)
% Calculates the jth direction of the Krylov subspace of (A-sigma*E) and v
% where v is the (j-1)th direction of the subspace.
% L, U are LU decompositions of (A-sigma*E) (can be empty if not used);
% p, q, D are additional vectors/matrices for sparse LU decomposition
%   (can be empty if not used).
% The Krylov subspace of interest is: K(inv(A-sigma*E)*E,(A-sigma*E)*B)

narginchk(10,10)
if isempty(L)
    %% Standard matlab solver (\)
    if sigma ~= inf
        % v = ((A - sigma*E)^(-1) * E)^(j-1) * (A-sigma*E)^(-1) * b
        if j == 1
            v = (A - sigma*E) \ v;
        else
            v = (A - sigma*E) \ E*v;
        end
    else    % match Markov parameters
        if j == 1
            return;
        else
            v = A*v;
        end
    end
elseif isempty(p)
    %% Use full LU decomposition for faster computation
    if sigma ~= inf
        % v = ((A - sigma*E)^(-1) * E)^(j-1) * (A-sigma*E)^(-1) * b
        if j == 1
            v = U\(L\v);
        else
            v = U\(L\(E*v));
        end
    else    % match Markov parameters
        if j == 1
            return;
        else
            v = A*v;
        end
    end
else
    %% Use sparse LU decomposition for even faster computation
    % p and q are permutation vectors for sparse LU decomposition
    % inv(D(:,p))*(A-sigma*E)(:,q) = L*U
    % (A-sigma*E)*v = b <=> v(q) = inv(U)*inv(L)*inv(D(:,p))*b
    % See documentation of lu for more information
    if sigma ~= inf
        % v = ((A - sigma*E)^(-1) * E)^(j-1) * (A-sigma*E)^(-1) * b
        if j == 1
            v(q,:) = U\(L\(D(:,p)\v));
        else
            v(q,:) = U\(L\(D(:,p)\(E*v)));
        end
    else    % match Markov parameters
        if j == 1
            return;
        else
            v = A*v;
        end
    end
end
end % of function nextSubspaceDirection

function v = gramSchmidt(v,V)
% Applies Gram-Schmidt orthogonalization to find v* with range([v,V]) = range([v*,V])
% which is orthogonal to all columns of V

if size(v,2) == 1   % single direction

    if ~isempty(V)  % check if V already has nonzero elements
        for i = 1:size(V,2)
            v = v - (v'*V(:,i)/norm(V(:,i)))*V(:,i);
        end
        v = v/norm(v);
    else
        v = v/norm(v);
    end

else            % multiple directions

    if ~isempty(V)  % check if V already has nonzero elements
        for j = 1:size(v,2)
            for i = 1:length(V(1,:))
                v(:,j) = v(:,j) - (v(:,j)'*V(:,i)/norm(V(:,i)))*V(:,i);
                v(:,j) = v(:,j)/norm(v(:,j));
            end
            V = [V, v(:,j)];
        end
    else
        v(:,1) = v(:,1)/norm(v(:,1));
        V = v(:,1);
        for j = 2:size(v,2)
            for i = 1:length(V(1,:))
                v(:,j) = v(:,j) - (v(:,j)'*V(:,i)/norm(V(:,i)))*V(:,i);
                v(:,j) = v(:,j)/norm(v(:,j));
            end
            V = [V, v(:,j)];
        end
    end

end

end % of function gramSchmidt

function [s0sorted, bSorted, multiplicities, redOrder] = processShifts(s0, b)
% Processes the row-vector of shifts and the column-matrix of tangent directions.
% Duplicates of the shift/tangent vector pairings will be removed.
% Complex conjugate pairs are represented by one instance.
% The multiplicities are stored in the row-vector multiplicities.
% The order which the reduced system will have is returned by redOrder.

if nargin < 2
    b = [];
end
% Matrix with shift & tangent direction pairings as ROWS
combined = [rowVector(s0); b]';
% (rows is required to use ismember 'rows' option)

% Matrix with shift and tangent direction pairings as ROWS
sorted = zeros(size(combined));
multiplicities = zeros(1, size(combined,1));

iRow = 0; % next free index in sorted
while ~isempty(combined)
    iRow = iRow + 1;
    row = combined(1,:);
    % Add new row to sorted
    sorted(iRow,:) = row;
    % Get indeces of duplicate and complex conjugate elements
    duplicateIndices = ismember(combined, row, 'rows');
    duplicateComplexIndices = ismember(combined, conj(row), 'rows');
    % Set multiplicities
    multiplicities(iRow) = max(sum(duplicateIndices),sum(duplicateComplexIndices));
    % Remove row and duplicates
    combined(max(duplicateIndices, duplicateComplexIndices),:) = [];
end

% Prepare output
s0sorted = sorted(1:iRow,1)';
bSorted = sorted(1:iRow,2:end)';
multiplicities = multiplicities(1:iRow);

redOrder = 0;
for i=1:length(multiplicities)
    if isreal(s0sorted(i)) && (isempty(bSorted) || isreal(bSorted(:,i)))
        redOrder = redOrder + multiplicities(i);
    else
        redOrder = redOrder + 2*multiplicities(i);
    end
end
end % of function sortShifts

function [s0_sort,b_sort] = sortShiftsForSSSMOR(s0,b)
% Sorting shifts for sssMOR->arnoldi: returns a sorted vector s0_sort in
% which complex shifts always appear in complex conjugate pairs.
% This function first uses 'accumlate' to get sorted s0 with multiplicities
% n. This information is then used to repeat the shifts (and b, if
% necessary) with 'repmat' according to multiplicities and complex
% properties.

if isempty(b)
    [s0,~,n] = processShifts(s0,b);
else
    [s0,n,b] = processShifts(s0,b);
end

indexComplex = find(imag(s0)~=0); % returns index of complex values
n_shifts = 0;

if isempty(b)
    % without tangent directions
    b_sort = [];
    for i = 1:length(n)
        if isempty(indexComplex)
            n_shifts = n_shifts+n(i);
            continue
        end
        if (max(i==indexComplex)==0)
            n_shifts = n_shifts+n(i);
        else
            n_shifts = n_shifts+2*n(i);
        end
    end
    s0_sort = zeros(1,n_shifts);
    idx = 0;
    for i = 1:length(s0)
        if isempty(indexComplex)
            s0_sort(idx+1:idx+n(i)) = repmat(s0(i),1,n(i));
            idx = idx + n(i);
            continue
        end
        if max(i==indexComplex)==0
            s0_sort(idx+1:idx+n(i)) = repmat(s0(i),1,n(i));
            idx = idx + n(i);
        else
            s0_sort(idx+1:idx+2*n(i)) = repmat([s0(i),conj(s0(i))],1,n(i));
            idx = idx + 2*n(i);
        end
    end
else
    % with tangent directions
    for i = 1:length(n)
        if isempty(indexComplex)
            n_shifts = n_shifts+n(i);
            continue
        end
        if (max(i==indexComplex)==0)
            n_shifts = n_shifts+n(i);
        else
            n_shifts = n_shifts+2*n(i);
        end
    end
    s0_sort = zeros(1,n_shifts);
    b_sort = zeros(size(b,1),n_shifts);
    idx = 0;
    for i = 1:length(s0)
        if isempty(indexComplex)
            s0_sort(idx+1:idx+n(i)) = repmat(s0(i),1,n(i));
            b_sort(:,idx+1:idx+n(i)) = repmat(b(:,i),1,n(i));
            idx = idx + n(i);
            continue
        end
        if max(i==indexComplex)==0
            s0_sort(idx+1:idx+n(i)) = repmat(s0(i),1,n(i));
            b_sort(:,idx+1:idx+n(i)) = repmat(b(:,i),1,n(i));
            idx = idx + n(i);
        else
            s0_sort(idx+1:idx+2*n(i)) = repmat([s0(i),conj(s0(i))],1,n(i));
            b_sort(:,idx+1:idx+2*n(i)) = repmat([b(:,i),conj(b(:,i))],1,n(i));
            idx = idx + 2*n(i);
        end
    end
end
end % of function sortShiftsForSSSMOR

function [s0, b, Opts] = parseInputs(sys, s0, varargin)

% Check phs input type
if ~isa(sys,'phs')
    error('MORpH:arnoldiPH:wrongInput', 'Original model is not an object of the phs-class');
end

% Catch DAE systems
if sys.isDAE
    error('MORpH:arnoldiPH:wrongInput', ...
        'ArnoldiPH currently only supports pHODE systems.');
end

% Opts
if ~isempty(varargin) && isstruct(varargin{end})
    Opts = varargin{end};
    varargin(end) = [];
else
    Opts = struct();
end

OptsAdmissible.constAlg = {'phs','sss'};
OptsAdmissible.orth = {'2mgs',false,'mgs','dgks'};
OptsAdmissible.structurePreservation = {'Cholesky','scaling','specialInverse','Cholesky+','QV',false};
OptsAdmissible.lse = {'sparse','full'};
OptsAdmissible.verbose = {true,false};
OptsAdmissible.phs = struct();
OptsAdmissible.phs.inputValidation = {true,false};
Opts = phsMOR_parseOpts(Opts,OptsAdmissible);

% s0 and b
if min(size(s0)) ~= 1
    error("MORpH:arnoldiPH:wrongShiftVector","s0 must be a row or column vector");
end
s0 = rowVector(s0);
if length(varargin) >= 1
    b = varargin{1};
    assert(size(b,1) == size(sys.G, 2),...
        'MORpH:arnoldiPH:dimGwrong',...
        'The tangent directions must have the same dimensions as the system input.');
    assert(size(s0,2) == size(b,2),...
        'MORpH:arnoldiPH:dimGwrong',...
        'You must provide the same number of tangent vectors and shifts.');
else
    b = [];
end

end % of function parseInputs