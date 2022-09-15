classdef test_phs < matlab.unittest.TestCase
    % TEST test_phs - Tests the functionalities of the phs class
    %
    % Press <F5> or enter "runtests("test_phs")" to run this testscript
    %
    % Description:
    %   This test verifies the functionalities of the phs class.
    %   These scenarios are considered:
    %       - Input: Test if all input combinations are processed correctly
    %       - Output: Test if all (optional) outputs are provided
    %       - Options: Test options
    %       - Input/Options: Test if options are effective to affect
    %               input processing
    %       - Class Functions: Test functionality of the phs class functions
    %               which are meant for use by user
    %       - Class Functions/Subunits: Test functionality of the phs class
    %               functions which are not meant for direct use by user
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

    properties
        % standard matrices which should be set by the phs class if no
        % parameter is provided
        standardJ = [0 -1; 1 0];
        standardR = [1 0; 0 0];
        standardQ = eye(2);
        standardG = [1; 0];
        standardE = eye(2);
        standardP = [0;0];
        standardS = 0;
        standardN = 0;
    end
    properties (TestParameter)
        % default test input matrices
        J = {[0 -1; 1 0]};
        R = {[1 0; 0 0]};
        Q = {eye(2)};
        G = {[1; 0]};
        E = {eye(2)};
        P = {[0;0]};
        S = {0};
        N = {0};
    end

    %% Setup / teardown of test environment
    methods (TestClassSetup) % once before all tests
        function lessOutput(testCase)
            % disable warnings to show cleaner output
            warning('off','all');
            testCase.addTeardown(@warning,'on','all');
        end
        function includeSystemSetupFunctions(testCase)
            % make sure that system folder with system setups is included
            if ~(exist('setup_MassSpringDamperSystem.m','file') == 2)
                error("You need to include ""...\MORpH\demos\test_systems"" in your MATLAB path to run the tests!")
            end
        end
    end

    methods (TestClassTeardown) % once after all tests

    end

    methods (TestMethodSetup) % before each test

    end

    methods (TestMethodTeardown) % after each test

    end

    %% TESTS

    methods (Test, TestTags = {'Input'})
        % Check input combinations

        function JRQGonly(testCase,J,R,Q,G,E,P,S,N)
            % basic input

            testSys = phs(J,R,Q,G);
            testCase.verifyProperties(testSys);
            testCase.verifyEqualSystemMatrices(testSys,J,R,Q,G,...
                testCase.standardE,testCase.standardP,testCase.standardS,testCase.standardN);
            testCase.verifyFalse(testSys.isMIMO);
            testCase.verifyFalse(testSys.isImplicit);
            testCase.verifyFalse(testSys.isDAE);
        end

        function MIMO(testCase,J,R,Q,G,E,P,S,N)
            % MIMO system input with matrices G, P, S, and N

            G = [1 0; 0 3];
            P = [0 0;0 0];
            S = [0 0;0 0];
            N = [0 0;0 0];

            testSys = phs(J,R,Q,G);
            testCase.verifyProperties(testSys);
            testCase.verifyEqualSystemMatrices(testSys,J,R,Q,G,...
                testCase.standardE,P,S,N);

            testSys = phs(J,R,Q,G,E,P,S,N);
            testCase.verifyProperties(testSys);
            testCase.verifyEqualSystemMatrices(testSys,J,R,Q,G,E,P,S,N)
        end

        function withE(testCase, J, R, Q, G, E, P, S, N)
            % SISO system with E specified

            testSys = phs(J,R,Q,G,E);
            testCase.verifyProperties(testSys);
            testCase.verifyEqualSystemMatrices(testSys,J,R,Q,G,E,P,S,N);
        end

        function JNotSkew(testCase, J, R, Q, G, E, P, S, N)
            % wrong input: J not skew symmetric --> phs should error

            J = [0 1; 1 0];
            Opts = struct();
            verifyError(testCase, @()phs(J,R,Q,G,E,P,S,N,Opts),'MORpH:phs:wrongInput');
        end

        function RNotSymm(testCase,J,R,Q,G,E,P,S,N)
            % wrong input: R not symmetric --> phs should error

            R = [1 0; 1 0];
            Opts = struct();
            verifyError(testCase, @()phs(J,R,Q,G,E,P,S,N,Opts),'MORpH:phs:wrongInput');
        end

        function RNotPositiveDefinite(testCase, J, R, Q, G, E, P, S, N)
            % wrong input: R not positive definite --> should error

            R = [-1 0; 0 0];
            Opts = struct();
            verifyError(testCase, @()phs(J,R,Q,G,E,P,S,N,Opts),'MORpH:phs:wrongInput');
        end

        function QNotPositiveSemiDefinite(testCase, J, R, Q, G, E, P, S, N)
            % wrong input: Q not positive semi-definite --> phs should error

            Q = -eye(2);
            Opts = struct();
            verifyError(testCase, @()phs(J,R,Q,G,E,P,S,N,Opts),'MORpH:phs:wrongInput');
        end

        function SNotSymm(testCase, J, R, Q, G, E, P, S, N)
            % wrong input: S not symmetric --> phs should error

            S = [1 0, 1 0];
            Opts = struct();
            verifyError(testCase, @()phs(J,R,Q,G,E,P,S,N,Opts),'MORpH:phs:wrongInput');
        end

        function wrongDimensionG(testCase, J, R, Q, G, E, P, S, N)
            % wrong input: G has different dimension --> phs should error

            G = [1; 0; 3];
            Opts = struct();
            verifyError(testCase, @()phs(J,R,Q,G,E,P,S,N,Opts),'MORpH:phs:wrongInput');
        end

        function descriptor(testCase, J, R, Q, G, E, P, S, N)
            % descriptor system input with E ~= identity matrix

            E = eye(2)*2;
            testSys = phs(J,R,Q,G,E,P,S,N);
            testCase.verifyProperties(testSys);
            testCase.verifyTrue(testSys.isImplicit);
            testCase.verifyEqualSystemMatrices(testSys,J,R,Q,G,E,P,S,N)
        end

        function DAE(testCase,J,R,Q,G,E,P,S,N)
            % input of a differential algebraic system with E not full-ranked

            E = [1 0; 0 0];
            testSys = phs(J,R,Q,G,E,P,S,N);
            testCase.verifyProperties(testSys);
            testCase.verifyTrue(testSys.isImplicit);
            testCase.verifyTrue(testSys.isDAE);
            testCase.verifyEqualSystemMatrices(testSys,J,R,Q,G,E,P,S,N);
        end

    end

    methods (Test, TestTags = {'Output'})
        % Check (optional) output

    end

    methods (Test, TestTags = {'Options'})
        % Check options

    end

    methods (Test, TestTags = {'Input','Options'})
        % Check input-option combinations

        function withOpts(testCase,J,R,Q,G,E)
            % additinal input: Opts-struct

            Opts.inputValidation = true;
            Opts.inputTolerance = 1e-5;
            testSys = phs(J,R,Q,G,E,Opts);
            testCase.verifyProperties(testSys);
            testCase.verifyEqual(testSys.Opts.inputTolerance, 1e-5);
            testCase.verifyEqualSystemMatrices(testSys,J,R,Q,G,...
                testCase.standardE,testCase.standardP,testCase.standardS,testCase.standardN);
        end

        function withMaximumInput(testCase,J,R,Q,G,E,P,S,N)
            % all possible input parameters provided

            Opts = struct();
            testSys = phs(J,R,Q,G,E,P,S,N,Opts);
            testCase.verifyProperties(testSys);
            testCase.verifyEqual(testSys.Opts.inputTolerance, 1e-10);
            testCase.verifyEqualSystemMatrices(testSys,J,R,Q,G,E,P,S,N);
        end

        function JNotSkewButWithoutValidation(testCase,J,R,Q,G,E,P,S,N)
            % wrong input: J not skew symmetric BUT validation deactivated -->
            % no errors expected!

            J = [0 1; 1 0];
            Opts = struct();
            Opts.inputValidation = false;

            try
                testSys = phs(J,R,Q,G,E,P,S,N,Opts);
                errored = false;
            catch
                errored = true;
            end
            if errored
                error('J has been validated');
            end
            testCase.verifyProperties(testSys);
            testCase.verifyEqualSystemMatrices(testSys,J,R,Q,G,E,P,S,N);
        end

        function highInputTolerance(testCase,J,R,Q,G,E,P,S,N)
            % Q is not positive semi-definite but tolerance has been increased

            Q = [1 0; 1e-3, -1]; % not positive semi-definite
            Opts = struct();
            Opts.inputTolerance = 1;

            try
                testSys = phs(J,R,Q,G,E,P,S,N,Opts);
                errored = false;
            catch
                errored = true;
            end
            if errored
                error('Validation errored despite of high tolerance');
            end
            testCase.verifyProperties(testSys);
            testCase.verifyEqualSystemMatrices(testSys,J,R,Q,G,E,P,S,N);
        end

        function sparsity(testCase, J, R, Q, G, E, P, S, N)

            sys = phs(J,R,Q,G,E,P,S,N);
            testCase.verifyTrue(issparse(sys.J));
            testCase.verifyTrue(issparse(sys.R));
            testCase.verifyTrue(issparse(sys.Q));
            testCase.verifyTrue(issparse(sys.G));
            testCase.verifyTrue(issparse(sys.E));
            testCase.verifyTrue(issparse(sys.P));
            testCase.verifyTrue(issparse(sys.S));
            testCase.verifyTrue(issparse(sys.N));

            sys.J = J;
            sys.R = R;
            sys.Q = Q;
            sys.G = G;
            sys.E = E;
            sys.P = P;
            sys.S = S;
            sys.N = N;
            testCase.verifyFalse(issparse(sys.J));
            testCase.verifyFalse(issparse(sys.R));
            testCase.verifyFalse(issparse(sys.Q));
            testCase.verifyFalse(issparse(sys.G));
            testCase.verifyFalse(issparse(sys.E));
            testCase.verifyFalse(issparse(sys.P));
            testCase.verifyFalse(issparse(sys.S));
            testCase.verifyFalse(issparse(sys.N));

            sys = sys.makeSparse();
            testCase.verifyTrue(issparse(sys.J));
            testCase.verifyTrue(issparse(sys.R));
            testCase.verifyTrue(issparse(sys.Q));
            testCase.verifyTrue(issparse(sys.G));
            testCase.verifyTrue(issparse(sys.E));
            testCase.verifyTrue(issparse(sys.P));
            testCase.verifyTrue(issparse(sys.S));
            testCase.verifyTrue(issparse(sys.N));

            sys = sys.makeFull();
            testCase.verifyFalse(issparse(sys.J));
            testCase.verifyFalse(issparse(sys.R));
            testCase.verifyFalse(issparse(sys.Q));
            testCase.verifyFalse(issparse(sys.G));
            testCase.verifyFalse(issparse(sys.E));
            testCase.verifyFalse(issparse(sys.P));
            testCase.verifyFalse(issparse(sys.S));
            testCase.verifyFalse(issparse(sys.N));
        end
    end

    methods (Test, TestTags = {'Class Functions'})
        % Check functions of phs-class

        function save_and_load(testCase, J, R, Q, G, E, P, S, N)
            % Save and load phs systems to/from file (this should not result in
            % empty system properties!)
            testSys = phs(J, R, Q, G, E, P, S, N);
            save('test.mat', 'testSys');
            clearvars testSys
            testCase.verifyEqual(exist('testSys','var'), 0)
            load('test.mat', 'testSys');
            testCase.verifyEqual(exist('testSys','var'), 1)
            testCase.verifyEqualSystemMatrices(testSys, J, R, Q, G, E, P, S, N);
            delete test.mat
        end

        function setter(testCase,J,R,Q,G,E,P,S,N)
            % set phs-properties after creation

            testSys = phs(); % Empty phs object
            testSys.J = J;
            testSys.R = R;
            testSys.Q = Q;
            testSys.G = G;
            testSys.E = E;
            testSys.P = P;
            testSys.S = S;
            testSys.N = N;
            testSys.inputValidation(testSys);
            testCase.verifyProperties(testSys);
            testCase.verifyEqual(testSys.Opts.inputTolerance, 1e-10);
            testSys = testSys.makeSparse();
            testCase.verifyEqualSystemMatrices(testSys,J,R,Q,G,E,P,S,N);
        end

        function setterWithWrongInput(testCase,J,R,Q,G,E,P,S,N)
            % set phs-properties after creation with WRONG inputs --> should
            % error

            J = [0 1; 1 0]; % not skew symmetric
            R = [-1 0; 0 0]; % not positive semi definite

            testSys = phs(); % empty phs object
            testSys.J = J;
            testSys.R = R;
            testSys.Q = Q;
            testSys.G = G;
            testSys.E = E;
            testSys.P = P;
            testSys.S = S;
            testSys.N = N;
            testCase.verifyError(@()testSys.inputValidation(testSys),'MORpH:phs:wrongInput');

            testSys = phs();
            testSys.Opts.verbose = false;
            testSys.Opts.inputValidation = false;
            testSys.J = J;
            testSys.R = R;
            testSys.Q = Q;
            testSys.G = G;
            testSys.E = E;
            testSys.P = P;
            testSys.S = S;
            testSys.N = N;
            % should error anyway (since it is explicitely called by user)
            testCase.verifyError(@()testSys.inputValidation(testSys),'MORpH:phs:wrongInput');
        end

        function testMakeExplicit(testCase,J,R,Q,G,E,P,S,N)
            % test phs-function "makeExplicit"

            E = 2*[1 0; 0 1]; % E is not identity matrix, E is non-singular
            testSys = phs(J,R,Q,G,E,P,S,N);
            testSys.Opts.verbose = false;
            testCase.verifyProperties(testSys);
            testCase.verifyTrue(testSys.isImplicit);
            testCase.verifyEqualSystemMatrices(testSys,J,R,Q,G,E,P,S,N);

            testSysExp = makeExplicit(testSys);
            % check that original system has not changed
            testCase.verifyProperties(testSys);
            testCase.verifyTrue(testSys.isImplicit, true);
            testCase.verifyEqualSystemMatrices(testSys,J,R,Q,G,E,P,S,N);
            % check that new system is now explicit
            testCase.verifyProperties(testSysExp);
            testCase.verifyFalse(testSysExp.isImplicit);
            testCase.verifyNotEqual(testSys.Q, testSysExp.Q);
            % check that transfer function has not changed
            sys_ss_original = ss(testSys);
            sys_ss_explicit = ss(testSysExp);
            testCase.verifyEqual(norm(tf(sys_ss_original)-tf(sys_ss_explicit)),0,'AbsTol',1e-13);
        end

        function convertToSS(testCase,J,R,Q,G,E,P,S,N)
            % test conversion from phs-object to ss-object

            Opts.verbose = false;
            testSys = phs(J,R,Q,G,E,P,S,N,Opts);
            sys_ss = ss(testSys);
            testCase.verifyClass(sys_ss, ?ss);
            testCase.verifyEqual(sys_ss.A, full((testSys.J-testSys.R)*testSys.Q));
            testCase.verifyEqual(sys_ss.B, full(testSys.G-testSys.P));
            testCase.verifyEqual(sys_ss.C, full((testSys.G+testSys.P)'*testSys.Q));
            testCase.verifyEqual(sys_ss.D, full(testSys.S+testSys.N));
        end

        function convertToSSS(testCase,J,R,Q,G,E,P,S,N)
            % test conversion from phs-object to sss-object

            testSys = phs(J,R,Q,G,E,P,S,N);
            sys_sss = sss(testSys);
            testCase.verifyClass(sys_sss, ?sss);
            testCase.verifyEqual(sys_sss.A, (testSys.J-testSys.R)*testSys.Q);
            testCase.verifyEqual(sys_sss.B, testSys.G-testSys.P);
            testCase.verifyEqual(sys_sss.C, (testSys.G+testSys.P)'*testSys.Q);
            testCase.verifyEqual(sys_sss.D, testSys.S+testSys.N);
        end

    end

    methods (Test, TestTags = {'Class Functions', 'Subunits'})
        % Test phs class functions which are not intended for explicit use by user

        function checkPositiveDefiniteness(testCase)
            % check additional phs-function "isPositiveDefinite"

            A = rand(10);
            A = tril(A)*tril(A)';   % this should be a positive definite matrix
            Opts.verbose = false;

            testCase.verifyTrue(phs.isPositiveDefinite(A,1e-13,0));
            testCase.verifyTrue(phs.isPositiveDefinite(A,0,0,Opts));

            A = [A, zeros(length(A),1);zeros(1,length(A)+1)]; % this should be a positive semi-definite matrix

            testCase.verifyFalse(phs.isPositiveDefinite(A,0,0));
            testCase.verifyFalse(phs.isPositiveDefinite(A,0,0,Opts));
            testCase.verifyTrue(phs.isPositiveDefinite(A,1e-13,1));
            testCase.verifyTrue(phs.isPositiveDefinite(A,1e-13,1,Opts));
            testCase.verifyTrue(phs.isPositiveDefinite(A,0,1,Opts));

            A(end,end) = -1e-12; % this should be an indefinite matrix
            testCase.verifyFalse(phs.isPositiveDefinite(A,0,0));
            testCase.verifyFalse(phs.isPositiveDefinite(A,0,0,Opts));
            testCase.verifyFalse(phs.isPositiveDefinite(A,1e-13,1));
            testCase.verifyFalse(phs.isPositiveDefinite(A,1e-13,1,Opts));
            testCase.verifyFalse(phs.isPositiveDefinite(A,0,1,Opts));

            testCase.verifyTrue(phs.isPositiveDefinite(A,1e-11,1));
        end

        function checkEnforcePHStructure(testCase)
            low = 1e-14;
            high = 1;
            Opts.inputTolerance = 1e-10;

            warning('on','all')

            % perturb J
            sys = setup_MassSpringDamperSystem(20,2,1,1,'MIMO');
            [J,R,Q,G,E,P,S,N] = sys.getMatrices();
            sys.Opts = Opts;
            sys.J(end,1) = sys.J(end,1) + low;
            phs.inputValidation(sys);   % should not error (within tolerance)
            phs.inputValidation(sys.enforcePHStructure);
            testCase.verifyWarningFree(@() sys.enforcePHStructure);
            sys.J = J;
            sys.J(end,1) = sys.J(end,1) + high;
            testCase.verifyError(@() phs.inputValidation(sys),'MORpH:phs:wrongInput');
            testCase.verifyWarning(@() sys.enforcePHStructure, "phs:enforcePHstructure:changedAboveTolerance");

            % perturb R
            sys = setup_MassSpringDamperSystem(20,2,1,1,'MIMO');
            [J,R,Q,G,E,P,S,N] = sys.getMatrices();
            sys.Opts = Opts;
            sys.R(end,1) = sys.R(end,1) + low;
            phs.inputValidation(sys);   % should not error (within tolerance)
            phs.inputValidation(sys.enforcePHStructure);
            testCase.verifyWarningFree(@() sys.enforcePHStructure);
            sys.R = R;
            sys.R(end,1) = sys.R(end,1) + high;
            testCase.verifyError(@() phs.inputValidation(sys),'MORpH:phs:wrongInput');
            testCase.verifyWarning(@() sys.enforcePHStructure, "phs:enforcePHstructure:changedAboveTolerance");

            % perturb Q
            sys = setup_MassSpringDamperSystem(20,2,1,1,'MIMO');
            [J,R,Q,G,E,P,S,N] = sys.getMatrices();
            sys.Opts = Opts;
            sys.Q(end,1) = sys.Q(end,1) + low;
            phs.inputValidation(sys);   % should not error (within tolerance)
            phs.inputValidation(sys.enforcePHStructure);
            testCase.verifyWarningFree(@() sys.enforcePHStructure);
            sys.Q = Q;
            sys.Q(end,1) = sys.Q(end,1) + high;
            testCase.verifyError(@() phs.inputValidation(sys),'MORpH:phs:wrongInput');
            testCase.verifyWarning(@() sys.enforcePHStructure, "phs:enforcePHstructure:changedAboveTolerance");

            % perturb G
            sys = setup_MassSpringDamperSystem(20,2,1,1,'MIMO');
            [J,R,Q,G,E,P,S,N] = sys.getMatrices();
            sys.Opts = Opts;
            sys.G(4) = sys.G(4) + high;
            phs.inputValidation(sys);   % should not error
            phs.inputValidation(sys.enforcePHStructure);
            [~,changes] = sys.enforcePHStructure();
            testCase.verifyEqual(changes.abs.G,0);
            testCase.verifyEqual(changes.abs.J,0);
            testCase.verifyEqual(changes.abs.R,0);
            testCase.verifyEqual(changes.abs.Q,0);
            testCase.verifyEqual(changes.abs.E,0);
            testCase.verifyEqual(changes.abs.P,0);
            testCase.verifyEqual(changes.abs.S,0);
            testCase.verifyEqual(changes.abs.N,0);
            testCase.verifyWarningFree(@() sys.enforcePHStructure);

            % perturb E
            sys = setup_MassSpringDamperSystem(20,2,1,1,'MIMO');
            [J,R,Q,G,E,P,S,N] = sys.getMatrices();
            sys.Opts = Opts;
            sys.E(end,1) = sys.E(end,1) + low;
            phs.inputValidation(sys);   % should not error (within tolerance)
            phs.inputValidation(sys.enforcePHStructure);
            testCase.verifyWarningFree(@() sys.enforcePHStructure);
            sys.E = E;
            sys.E(end,1) = sys.E(end,1) + high;
            testCase.verifyError(@() phs.inputValidation(sys),'MORpH:phs:wrongInput');
            testCase.verifyWarning(@() sys.enforcePHStructure, "phs:enforcePHstructure:changedAboveTolerance");

            % perturb P
            sys = setup_MassSpringDamperSystem(20,2,1,1,'MIMO');
            [J,R,Q,G,E,P,S,N] = sys.getMatrices();
            sys.Opts = Opts;
            sys.P(4) = sys.P(end,1) + low;
            phs.inputValidation(sys);   % should not error (within tolerance)
            phs.inputValidation(sys.enforcePHStructure);
            testCase.verifyWarningFree(@() sys.enforcePHStructure);
            sys.P = P;
            sys.P(4) = sys.P(4) + high;
            testCase.verifyError(@() phs.inputValidation(sys),'MORpH:phs:wrongInput');
            testCase.verifyWarning(@() sys.enforcePHStructure, "phs:enforcePHstructure:changedAboveTolerance");


            % perturb S
            sys = setup_MassSpringDamperSystem(20,2,1,1,'MIMO');
            [J,R,Q,G,E,P,S,N] = sys.getMatrices();
            sys.Opts = Opts;
            sys.S(end,1) = sys.S(end,1) + low;
            phs.inputValidation(sys.enforcePHStructure);
            testCase.verifyWarning(@() sys.enforcePHStructure, "phs:enforcePHstructure:changedAboveTolerance");
            sys.S = S;
            sys.S(end,1) = sys.S(end,1) + high;
            testCase.verifyError(@() phs.inputValidation(sys),'MORpH:phs:wrongInput');
            testCase.verifyWarning(@() sys.enforcePHStructure, "phs:enforcePHstructure:changedAboveTolerance");

            % perturb N
            sys = setup_MassSpringDamperSystem(20,2,1,1,'MIMO');
            [J,R,Q,G,E,P,S,N] = sys.getMatrices();
            sys.Opts = Opts;
            sys.N(end,1) = sys.N(end,1) + low;
            phs.inputValidation(sys.enforcePHStructure);
            testCase.verifyWarning(@() sys.enforcePHStructure, "phs:enforcePHstructure:changedAboveTolerance");
            sys.N = N;
            sys.N(end,1) = sys.N(end,1) + high;
            testCase.verifyError(@() phs.inputValidation(sys),'MORpH:phs:wrongInput');
            testCase.verifyWarning(@() sys.enforcePHStructure, "phs:enforcePHstructure:changedAboveTolerance");

            % Descriptor system (perturb Q and E
            sys = setup_MassSpringDamperSystem(20,2,1,1,'MIMO');
            sys.E = 2*sys.E;
            sys.Opts = Opts;
            [J,R,Q,G,E,P,S,N] = sys.getMatrices();
            phs.inputValidation(sys);

            sys.E(end,1) = sys.E(end,1) + low;
            sys.Q(end,1) = sys.Q(end,1) + low;
            phs.inputValidation(sys);   % should not error (within tolerance)
            phs.inputValidation(sys.enforcePHStructure);
            testCase.verifyWarningFree(@() sys.enforcePHStructure);
            sys.E = E;
            sys.Q = Q;
            sys.E(end,1) = sys.E(end,1) + high;
            sys.Q(end,1) = sys.Q(end,1) + high;
            testCase.verifyError(@() phs.inputValidation(sys),'MORpH:phs:wrongInput');
            testCase.verifyWarning(@() sys.enforcePHStructure, "phs:enforcePHstructure:changedAboveTolerance");


            warning('off','all')
        end

        function checkScaling(testCase, J, R, G)
            % Check scaling function (from utility) which ensures Q == I
            Q = ones(size(J)) + eye(size(J));
            sys = phs(J, R, Q, G);
            sys_scaled = scaling(sys);
            testCase.verifyEqual(sys_scaled.Q, speye(size(Q)));
            testCase.verifyLessThan(norm(sys-sys_scaled), 1e-13);
        end

    end

    %% Supporting functions
    methods
        % Test support functions for reuse
        function verifyProperties(testCase,testSys)
            % Check if properties have correct type
            testCase.verifyTrue(ismatrix(testSys.J));
            testCase.verifyTrue(ismatrix(testSys.R));
            testCase.verifyTrue(ismatrix(testSys.Q));
            testCase.verifyTrue(ismatrix(testSys.G));
            testCase.verifyTrue(ismatrix(testSys.E));
            testCase.verifyTrue(ismatrix(testSys.P));
            testCase.verifyTrue(ismatrix(testSys.S));
            testCase.verifyTrue(ismatrix(testSys.N));
            testCase.verifyTrue(isstruct(testSys.Opts));
        end
        function verifyEqualSystemMatrices(testCase, testSys, J, R, Q, G, E, P, S, N)
            % Check if properties have correct value
            testCase.verifyEqual(testSys.J,sparse(J));
            testCase.verifyEqual(testSys.R,sparse(R));
            testCase.verifyEqual(testSys.Q,sparse(Q));
            testCase.verifyEqual(testSys.G,sparse(G));
            testCase.verifyEqual(testSys.P,sparse(P));
            testCase.verifyEqual(testSys.E,sparse(E));
            testCase.verifyEqual(testSys.S,sparse(S));
            testCase.verifyEqual(testSys.N,sparse(N));
        end
    end
end