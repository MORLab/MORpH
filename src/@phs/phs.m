classdef phs
    % PHS (class) - class and methods for handling port-Hamiltonian systems
    %
    % Syntax:
    %   sys = PHS(J, R, Q, G)
    %   sys = PHS(J, R, Q, G, E)
    %   sys = PHS(J, R, Q, G, Opts)
    %   sys = PHS(J, R, Q, G, E, Opts)
    %   sys = phs(J, R, Q, G, E, P, S, N)
    %   sys = phs(J, R, Q, G, E, P, S, N, Opts)
    %
    % Description:
    %       sys = phs(J, R, Q, G) returns a phs-object representing the
    %       port-Hamiltonian system sys defined by
    %           dx/dt = (J - R)*Q*x(t) + G*u(t)
    %               y = G'*Q*x
    %            H(x) = 0.5*x'*Q*x      (Hamiltonian)
    %               f = -dx/dt     (flow)
    %               e = grad_x(H)  (effort)        [1,2]
    %
    %       sys = phs(J, R, Q, G, E) returns a phs-object representing the
    %       port-Hamiltonian system sys defined by
    %           E*dx/dt = (J - R)*Q*x(t) + G*u(t)
    %                 y = G'*Q*x
    %              H(x) = 0.5*x'*E'*Q*x   (Hamiltonian)
    %                 f = -dx/dt     (flow)
    %                 e = grad_x(H)  (effort)      [3,4]
    %
    %       sys = phs(J, R, Q, G, E, P, S, N) returns a phs-object representing
    %       the port-Hamiltonian system sys defined by
    %           E*dx/dt = (J - R)*Q*x(t) + (G - P)*u(t)
    %                 y = (G + P)'*Q*x + (S + N)*u
    %              H(x) = 0.5*x'*E'*Q*x   (Hamiltonian)
    %                 f = -dx/dt     (flow)
    %                 e = grad_x(H)  (effort)      [3,4]
    %
    %       For detailed information on matrix properties, see [3, Definition 5].
    %
    %
    % Input Arguments:
    %       *Required Input Arguments:*
    %       - J:        real skew-symmetric matrix describing the interconnection
    %                   structure of the system
    %       - R:        real symmetric positive semidefinite matrix describing the
    %                   dissipative behaviour of the system
    %       - Q:        real symmetric positive semidefinite matrix describing the
    %                   energy storage behaviour of the system. (For
    %                   descriptor systems, restrictions differ slightly. -> see
    %                   'E')
    %       - G:        real input/output matrix
    %
    %       *Optional Input Arguments:*
    %       - E:        real (non-)singular descriptor matrix (Q'*E must be
    %                   symmetric and positive semi-definite)
    %       - P:        additional real input/output matrix
    %       - S:        real symmetric positive semidefinite feedthrough matrix
    %       - N:        real skew-symmetric feedthrough matrix
    %       - Opts:     structure with options:
    %           - .inputValidation:     switch input validation on (true) or off
    %                                   (false); This may speed up instantiation
    %                                   but faulty input is not detected any
    %                                   more.
    %                                   [{true} / false]
    %           - .inputTolerance:      Tolerance for input validation. If,
    %                                   for example, Q has eigenvalues close to
    %                                   zero in the complex left half plane, it will
    %                                   still be considered positive semidefinite
    %                                   if their (negative) real part is within
    %                                   Opts.inputTolerance. 
    %                                   [{1e-10} / positive double]
    %           - .verbose:             Turn off warnings by setting to false.
    %                                   [{true} / false]
    %
    % Output Arguments:
    %       - sys:      phs-object
    %
    % Examples:
    %       The following code creates a simple phs object:
    %
    %       J = [0 -1; 1 0]; R = [0 0; 0 1]; Q = eye(2); G = [1; 0];
    %       sys = phs(J, R, Q, G);
    %
    % See Also:
    %       ss, sss, demo_phs_class
    %
    % References:
    %       [1] V. Duindam, A. Macchelli, S. Stramigioli, and H. Bruyninckx. Modeling and control of
    %            complex physical systems: the port-Hamiltonian approach. Springer-Verlag, Berlin, Heidelberg, 2009.
    %       [2] A. van der Schaft and D. Jeltsema. Port-Hamiltonian systems theory: An introductory overview.
    %           Foundations and Trends in Systems and Control, 1(2-3):173–378, 2014.
    %       [3] C. Beattie, V. Mehrmann, H. Xu, and H. Zwart. Port-Hamiltonian descriptor systems. Math.
    %           Control Signals Systems, 30(17):1–27, 2018.
    %       [4] C. Mehl, V. Mehrmann, and M. Wojtylak. Linear algebra properties of dissipative Hamiltonian
    %           descriptor systems. SIAM J. Matrix Anal. Appl., 39(3):1489–1519, 2018.
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

    %% Properties
    properties (Access = public)
        J
        R
        Q
        G
        E
        P
        S
        N
        Opts
    end
    properties (Access = private)
        flag_descriptor
        flag_MIMO
        flag_DAE
        flag_hasStaircase
    end
    properties (Dependent)
        dim (1,1)
        isMIMO
        isImplicit
        isDAE
        hasStaircase
    end
    %end of properties

    %% Methods
    methods

        function sys = phs(J, R, Q, G, varargin)
            % Creates phs object

            % Empty constructor (required by MATLAB)
            if nargin == 0
                sys.Opts = struct();
                sys.Opts.inputValidation = true;
                sys.Opts.inputTolerance = 1e-10;
                sys.Opts.verbose = true;
                sys.J = []; sys.R = []; sys.Q = []; sys.G = [];
                sys.E = []; sys.P = []; sys.S = []; sys.N = [];
                return
            end

            % Parse inputs
            [J, R, Q, G, E, P, S, N, Opts] = phs.parseInputs(J, R, Q, G, varargin{:});

            % Save system properties
            % (Make matrices sparse to save computation time and memory)
            sys.Opts = Opts;
            sys.J = sparse(J);
            sys.R = sparse(R);
            sys.Q = sparse(Q);
            sys.G = sparse(G);
            sys.E = sparse(E);
            sys.P = sparse(P);
            sys.S = sparse(S);
            sys.N = sparse(N);

            if Opts.inputValidation
                phs.inputValidation(sys);
                % inputValidation function is static (see methods (Static))
                % and specified in external file in @phs class folder
            else
                if Opts.verbose
                    warning('MORpH:phs:noInputValidation',...
                        'No input validation for correctness of port-Hamiltonian system!')
                end
            end

        end


        %% Getter and Setter

        function J = get.J(sys); J = sys.J; end
        function R = get.R(sys); R = sys.R; end
        function Q = get.Q(sys); Q = sys.Q; end
        function G = get.G(sys); G = sys.G; end
        function E = get.E(sys); E = sys.E; end
        function P = get.P(sys); P = sys.P; end
        function S = get.S(sys); S = sys.S; end
        function N = get.N(sys); N = sys.N; end
        function dim = get.dim(sys); dim = length(sys.E); end

        function sys = set.Opts(sys,Opts)
            OptsAdmissible.inputValidation = {true, false};
            OptsAdmissible.verbose = {true, false};
            OptsAdmissible.inputTolerance = 1e-10;
            sys.Opts = phsMOR_parseOpts(Opts,OptsAdmissible);
        end

        function sys = set.J(sys, J)
            if ~isempty(sys.J) && sys.Opts.verbose
                warning('phs:phs:ChangedProperty',...
                    'Property has been changed successfully.\nConsider running phs.inputValidation(sys) to validate new properties');
            end
            sys.J = J;
            sys = sys.updateFlags('J');     % Set flags depending on J
        end

        function sys = set.R(sys, R)
            if ~isempty(sys.R) && sys.Opts.verbose
                warning('phs:phs:ChangedProperty',...
                    'Property has been changed successfully.\nConsider running phs.inputValidation(sys) to validate new properties');
            end
            sys.R = R;
            sys = sys.updateFlags('R');     % Set flags depending on R
        end
        function sys = set.Q(sys, Q)
            if ~isempty(sys.Q) && sys.Opts.verbose
                warning('phs:phs:ChangedProperty',...
                    'Property has been changed successfully.\nConsider running phs.inputValidation(sys) to validate new properties');
            end
            sys.Q = Q;
            sys = sys.updateFlags('Q');     % Set flags depending on Q
        end

        function sys = set.G(sys, G)
            if ~isempty(sys.G) && sys.Opts.verbose
                warning('phs:phs:ChangedProperty',...
                    'Property has been changed successfully.\nConsider running phs.inputValidation(sys) to validate new properties');
            end
            sys.G = G;
            sys = sys.updateFlags('G');     % Set flags depending on G
        end

        function sys = set.E(sys, E)
            if ~isempty(sys.E) && sys.Opts.verbose
                warning('phs:phs:ChangedProperty',...
                    'Property has been changed successfully.\nConsider running phs.inputValidation(sys) to validate new properties');
            end
            sys.E = E;
            sys = sys.updateFlags('E');     % Set flags depending on E
        end

        function sys = set.P(sys, P)
            if ~isempty(sys.P) && sys.Opts.verbose
                warning('phs:phs:ChangedProperty',...
                    'Property has been changed successfully.\nConsider running phs.inputValidation(sys) to validate new properties');
            end
            sys.P = P;
            sys = sys.updateFlags('P');     % Set flags depending on P
        end

        function sys = set.S(sys, S)
            if ~isempty(sys.S) && sys.Opts.verbose
                warning('phs:phs:ChangedProperty',...
                    'Property has been changed successfully.\nConsider running phs.inputValidation(sys) to validate new properties');
            end
            sys.S = S;
            sys = sys.updateFlags('S');     % Set flags depending on S
        end

        function sys = set.N(sys, N)
            if ~isempty(sys.N) && sys.Opts.verbose
                warning('phs:phs:ChangedProperty',...
                    'Property has been changed successfully.\nConsider running phs.inputValidation(sys) to validate new properties');
            end
            sys.N = N;
            sys = sys.updateFlags('N');     % Set flags depending on N
        end

        function [J, R, Q, G, E, P, S, N] = getMatrices(sys)
            % GETMATRICES - returns the matrices of the system
            %   Syntax:
            %       [J, R, Q, G, E, P, S, N] = GETMATRICES(sys)
            J = sys.J;
            R = sys.R;
            Q = sys.Q;
            G = sys.G;
            E = sys.E;
            P = sys.P;
            S = sys.S;
            N = sys.N;
        end

        function sys = set.dim(sys, ~)
            error('phs:phs:cannotSetProperty', 'This property (dim) is set automatically.');
        end

        function isImplicit = get.isImplicit(sys)
            isImplicit =  sys.flag_descriptor;
        end

        function sys = set.isImplicit(sys, ~)
            error('phs:phs:cannotSetProperty', 'This property (isImplicit) is set automatically.');
        end

        function isMIMO = get.isMIMO(sys)
            isMIMO = sys.flag_MIMO;
        end

        function sys = set.isMIMO(sys, ~)
            error('phs:phs:cannotSetProperty', 'This property (isMIMO) is set automatically.');
        end

        function isDAE = get.isDAE(sys)
            isDAE = sys.flag_DAE;
        end

        function sys = set.isDAE(sys, ~)
            error('phs:phs:cannotSetProperty', 'This property (isDAE) is set automatically.');
        end

        function hasStaircase = get.hasStaircase(sys)
            hasStaircase = sys.flag_hasStaircase;
        end

        function sys = set.hasStaircase(sys, ~)
            error('phs:phs:cannotSetProperty', 'This property (hasStaircase) is set automatically.');
        end

        %% Converter

        function sys = phs2ss(sys)
            % PHS2SS - transforms a phs system to regular ss system
            %   Syntax:
            %       sys_ss = phs2ss(sys)
            %   Input arguments:    sys - phs-object
            %   Output arguments:   sys_ss - ss-object
            warning('off','Control:ltiobject:SSSparse2Full');
            if ~sys.flag_descriptor
                sys = ss((sys.J - sys.R)*sys.Q, sys.G-sys.P, (sys.G+sys.P)'*sys.Q, sys.S+sys.N);
            else
                % if sys.Opts.verbose, disp('System is descriptor system'), end
                sys = dss((sys.J - sys.R)*sys.Q, sys.G-sys.P, (sys.G+sys.P)'*sys.Q, sys.S+sys.N, sys.E);
            end
            warning('on','Control:ltiobject:SSSparse2Full');
        end

        function sys = ss(sys)
            % SS - transforms a phs system to regular ss system
            %   Syntax:
            %       sys_ss = ss(sys)
            %   Input arguments:    sys - phs-object
            %   Output arguments:   sys_ss - ss-object
            sys = phs2ss(sys);
        end

        function sys = phs2sss(sys)
            % PHS2SSS - transforms a phs system to a sss system
            %   Syntax:
            %       sys_sss = phs2sss(sys)
            %   Input arguments:    sys - phs-object
            %   Output arguments:   sys_sss - sss-object
            sys = sss((sys.J - sys.R)*sys.Q, sys.G-sys.P, (sys.G+sys.P)'*sys.Q, sys.S+sys.N, sys.E);
        end

        function sys = sss(sys)
            % SSS - transforms a phs system to a sss system
            %   Syntax:
            %       sys_sss = sss(sys)
            %   Input arguments:    sys - phs-object
            %   Output arguments:   sys_sss - sss-object
            sys = phs2sss(sys);
        end

        function sys = makeFull(sys)
            % MAKEFULL - returns a system that has only full system matrices
            %   Syntax:
            %       sys_full = makeFull(sys)
            %   Input arguments:    sys - phs-object
            %   Output arguments:   sys_full - phs-object
            sys.J = full(sys.J); sys.R = full(sys.R); sys.Q = full(sys.Q); sys.G = full(sys.G);
            sys.E = full(sys.E); sys.P = full(sys.P); sys.S = full(sys.S); sys.N = full(sys.N);
        end

        function sys = makeSparse(sys)
            % MAKESPARSE - returns a system that has only sparse system matrices
            %   Syntax:
            %       sys_full = makeSparse(sys)
            %   Input arguments:    sys - phs-object
            %   Output arguments:   sys_full - phs-object
            sys.J = sparse(sys.J); sys.R = sparse(sys.R); sys.Q = sparse(sys.Q); sys.G = sparse(sys.G);
            sys.E = sparse(sys.E); sys.P = sparse(sys.P); sys.S = sparse(sys.S); sys.N = sparse(sys.N);
        end

        %% Operator overloading

        sys = minus(sys1,sys2);     % sys1-sys2; in external file
        sys = mtimes(sys1,sys2);    % sys1*sys2; in external file
        sys = plus(sys1,sys2);      % sys1+sys2; in external file

    end

    %% Private methods
    methods (Access = private)
        function sys = updateFlags(sys, changedPropertyKey)
            % Update private flags of the system according to the property that
            % was changed. changedPropertyKey refers to the changed property
            % (so far 'E' and 'G' are implemented).
            switch changedPropertyKey
                case 'E'
                    % Check if system is descriptor or DAE system
                    if isequal(sys.E, eye(size(sys.E)))
                        sys.flag_descriptor = false;
                        sys.flag_DAE = false;
                    else
                        sys.flag_descriptor = true;
                        if svds(sys.E, 1, 'smallest') < 1e-12
                            sys.flag_DAE = true;
                        else
                            sys.flag_DAE = false;
                        end
                    end
                    sys.flag_hasStaircase = staircaseCheck(sys);
                case 'G'
                    % Check if system is MIMO system
                    if size(sys.G,2) > 1
                        sys.flag_MIMO = true;
                    else
                        sys.flag_MIMO = false;
                    end
                    sys.flag_hasStaircase = staircaseCheck(sys);
                case {'J','R','Q','P'}
                    sys.flag_hasStaircase = staircaseCheck(sys);
                otherwise
            end
        end
    end

    %% Static methods
    methods (Static)

        function [J, R, Q, G, E, P, S, N, Opts] = parseInputs(J, R, Q, G, varargin)
            narginchk(4, 9)
            s = size(G);

            % Opts
            Opts = struct();
            if ~isempty(varargin) && isstruct(varargin{end})
                Opts = varargin{end};
                varargin = varargin(1:end-1);
            end

            OptsAdmissible.inputValidation = {true, false};
            OptsAdmissible.verbose = {true, false};
            OptsAdmissible.inputTolerance = 1e-10;
            Opts = phsMOR_parseOpts(Opts,OptsAdmissible);

            % Optional system matrices
            switch length(varargin)
                case 0
                    E = speye(length(J));
                    P = zeros(s);
                    S = zeros(s(2),s(2));
                    N = zeros(s(2),s(2));
                case 1
                    E = varargin{:};
                    P = zeros(s);
                    S = zeros(s(2),s(2));
                    N = zeros(s(2),s(2));
                case 2
                    [E,P] = varargin{:};
                    S = zeros(s(2),s(2));
                    N = zeros(s(2),s(2));
                case 3
                    [E, P, S] = varargin{:};
                    N = zeros(s(2),s(2));
                case 4
                    [E, P, S, N] = varargin{:};
                otherwise
                    error('MORpH:phs:parseInput', 'Input combination is not specified!')
            end
        end

    end

    %% Signatures of functions in separate files (see folder @phs)
    methods
        [varargout] = bode(varargin);                       % Bode-function
        bodemag(varargin);                                  % Bodemag-function
        [varargout] = eig(sys, varargin);                   % Eigenvalues
        [varargout] = eigs(sys, varargin);                  % Eigenvalues
        [sys,sys_old,changes] = enforcePHStructure(sys, proceedAnyway)  % Make sure the system is PH
        fresp = evalfr(sys, fr);                            % Transfer function evaluation
        fbSys = feedback(sys_plant, sys_contr)              % Feedback interconnection for pH systems
        sys_frd = frd(sys, varargin);                       % Returns frd (Frequency Response Data) object
        [H, W] = freqresp(sys, varargin);                   % Frequency response
        sysTr = gyrator(sys1, sys2, M);                     % Transformer interconnection for pH systems
        [varargout] = impulse(varargin);                    % Impulse-function (time-domain simulation)
        [varargout] = lsim(varargin);                       % Simulation of a system with custom input
        [R, S] = lyapchol(sys);                             % Solve systems Lyapunov equation
        sys = makeExplicit(sys,varargin);                   % Transform implicit to explicit system (E=I)
        [n, varargout] = norm(varargin);                    % H_2- and L_infinity-norm
        [varargout] = pzmap(varargin);                      % Poles and zeros
        [r, p, d] = residue(sys);                           % Poles, residue, and feedthrough
        dims = staircaseDims(sys);                          % Computes state decomposition for staircase system
        hasStaircase = staircaseCheck(sys,varargin);        % Checks if system is in staircase form
        sys = scaling(sys, varargin)                        % Transform system to scaled energy coordinates (Q=I)
        [tseries, Y, X] = simulate(sys, u, tseries, x0, Opts);  % Simulate the system response
        [varargout] = step(varargin);                       % Step function (time-domain simulation)
        figureHandle = spy(sys,figureHandle);               % Plot spy for all system matrices
        sysTr = transformer(sys1, sys2, M);                 % Transformer interconnection for pH systems
        G = transferFunction(sys);                          % Transfer function of sys
    end
    methods (Static)
        result = isPositiveDefinite(A,tol,semi,varargin);   % Returns true if  A is positive definite
        isPH = inputValidation(sys);                        % Check for PH structure of the system
    end

end