classdef test_irkaPH < matlab.unittest.TestCase
    % TEST test_irkaPH - Tests the functionalities of irkaPH
    %
    % Press <F5> or enter "runtests("test_irkaPH")" to run this testscript
    %
    % Description:
    %   This script tests the function 'irkaPH' of the MORpH-toolbox. The
    %   following test scenarios are considered:
    %       - Input: Check if function can run all intended input combinations
    %       - Output: Check if each additional output can be obtained
    %       - Options: Check if additional options work correctly
    %       - Functionality: Try more complex reduction to simulate real-world
    %                       application
    %
    %-----------------------------------------------------------------------
    % This file is part of
    %
    % <a href="https://github.com/MORLab/MORpH">MORpH</a> - a MATLAB toolbox to store, analyze,
    % interconnect and reduce large-scale port-Hamiltonian models
    %
    % Authors:      Julius Durmann
    % E-Mail:       <a href="mailto:morlab.rt@ed.tum.de">morlab.rt@ed.tum.de</a>
    % Website:      <a href="https://www.epc.ed.tum.de/en/rt/home">www.epc.ed.tum.de/rt</a>
    % Copyright :   © 2022 Chair of Automatic Control, TUM (see LICENSE.md)
    %-----------------------------------------------------------------------

    %% Properties
    properties (TestParameter)
        sys = {setup_MassSpringDamperSystem(40,2,1,1)};   % standard MSD system for most of the tests
        sysMIMO = {setup_MassSpringDamperSystem(40,2,1,1,'MIMO')}; % MIMO MSD system for MIMO testing
        shifts = {logspace(-5,5,5)}; % standard shifts
    end

    properties
        tolerance = 1e-10;
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


    %% Tests
    methods (Test, TestTags = {'Input'})
        % check input combinations

        function standardInput(testCase,sys)
            % Test standard input (system + reduced order)

            redSys = irkaPH(sys,5);
            testCase.verifyEqual(redSys.dim, 5);
        end

        function standardInputWithOpts(testCase,sys)
            % Test standard input with simple options

            opts.arnoldiPH.constAlg = 'phs';
            opts.arnoldiPH.structurePreservation = 'QV';
            opts.tol = 1e-3;
            opts.arnoldiPH.phs.verbose = true;
            redSys = irkaPH(sys,5,opts);
            testCase.verifyEqual(redSys.dim, 5);
            testCase.verifyEqual(redSys.parameters.arnoldiPH.constAlg, opts.arnoldiPH.constAlg);
            testCase.verifyEqual(redSys.parameters.arnoldiPH.structurePreservation, opts.arnoldiPH.structurePreservation);
            testCase.verifyEqual(redSys.parameters.tol, opts.tol);
            testCase.verifyEqual(redSys.Opts.verbose, opts.arnoldiPH.phs.verbose);
        end

        function initialShifts(testCase,sys,shifts)
            % Supply initial shifts

            redSys = irkaPH(sys, shifts);
            testCase.verifyEqual(redSys.dim, 5);
            testCase.verifyEqual(redSys.parameters.initShifts, 'custom');
            testCase.verifyEqual(redSys.parameters.startShifts, shifts);
        end

        function MIMO(testCase,sysMIMO)
            % Test MIMO system reduction
            opts.constAlg = 'phs';
            redSys = irkaPH(sysMIMO,20,opts);
            testCase.verifyEqual(redSys.dim, 20);
        end

        function MIMOWithInitialShifts(testCase,sysMIMO)
            % Test MIMO system reduction with initial shifts supplied

            shifts = logspace(-10,10,20);
            opts.constAlg = 'phs';
            redSys = irkaPH(sysMIMO,shifts,opts);
            testCase.verifyEqual(redSys.dim, 20);
            testCase.verifyEqual(redSys.parameters.initShifts, 'custom');
            testCase.verifyEqual(redSys.parameters.startShifts, shifts);
        end

        function MIMOWithInitialShiftsAndTangent(testCase,sysMIMO)
            % Test MIMO system reduction with initial shifts and tangent
            % directions supplied

            shifts = logspace(-10,10,5);
            tangent = rand(2,5);
            opts.constAlg = 'phs';
            redSys = irkaPH(sysMIMO,shifts,tangent,opts);
            testCase.verifyEqual(redSys.dim, 5);
            testCase.verifyEqual(redSys.parameters.initShifts, 'custom');
            testCase.verifyEqual(redSys.parameters.startShifts, shifts);
            testCase.verifyEqual(redSys.parameters.startTangent, tangent);
        end

        function MIMOWithComplexShiftsAndTangent(testCase,sysMIMO)
            % Test MIMO system reduction with complex initial shifts and tangent
            % directions supplied

            shifts = logspace(-10,10,5) + 1i*logspace(-10,10,5);
            tangent = rand(2,5);
            opts.constAlg = 'phs';
            redSys = irkaPH(sysMIMO,shifts,tangent,opts);
            testCase.verifyEqual(redSys.dim, 10);
            testCase.verifyEqual(redSys.parameters.initShifts, 'custom');
            testCase.verifyEqual(redSys.parameters.startShifts, shifts);
            testCase.verifyEqual(redSys.parameters.startTangent, tangent);
        end

    end

    methods (Test, TestTags = {'Output'})
        % check (optional) output

        function standardInputWithOptionalOutput(testCase,sys)
            % Verify that optional outputs are provided

            [redSys, V, s0, b, W, nLU] = irkaPH(sys,5);
            testCase.verifyEqual(redSys.dim, 5);
            testCase.verifySize(V,[sys.dim,5]);
            testCase.verifyEqual(rank(V), 5);
            testCase.verifySize(s0,[5,1]);
            testCase.verifyEqual(b, ones(1,5));
            testCase.verifySize(W, [sys.dim,5]);
            testCase.verifyGreaterThan(nLU, 0)
        end

        function returnedShiftsAreCorrect(testCase, sys)
            % Verify that the (optionally) returned shift vector corresponds to
            % the reduced system

            opts.arnoldiPH.constAlg = 'phs';
            opts.arnoldiPH.structurePreservation = 'Cholesky+';
            [redSys, V, s0] = irkaPH(sys,5,opts);

            [redSys_reconstructed,V_reconstructed] = arnoldiPH(sys,s0,opts.arnoldiPH);

            testCase.verifyEqual(redSys.J, redSys_reconstructed.J);
            testCase.verifyEqual(redSys.R, redSys_reconstructed.R);
            testCase.verifyEqual(redSys.Q, redSys_reconstructed.Q);
            testCase.verifyEqual(redSys.G, redSys_reconstructed.G);
            testCase.verifyEqual(V, V_reconstructed);
        end

        function returnedShiftsAndTangentAreCorrect(testCase, sys)
            % Verify that the (optionally) returned shift vector and tangent
            % vector matrix corresponds to the reduced system

            opts.arnoldiPH.constAlg = 'phs';
            opts.arnoldiPH.structurePreservation = 'Cholesky+';
            [redSys, V, s0, b] = irkaPH(sys,5,opts);

            [redSys_reconstructed,V_reconstructed] = arnoldiPH(sys,s0,b,opts.arnoldiPH);

            testCase.verifyEqual(redSys.J, redSys_reconstructed.J);
            testCase.verifyEqual(redSys.R, redSys_reconstructed.R);
            testCase.verifyEqual(redSys.Q, redSys_reconstructed.Q);
            testCase.verifyEqual(redSys.G, redSys_reconstructed.G);
            testCase.verifyEqual(V, V_reconstructed);
        end

        function redSysMethod(testCase, sys)
            redSys = irkaPH(sys,5);
            testCase.verifyEqual(redSys.method, @irkaPH);
        end

        function redSysParameters(testCase,sys)
            % Verify that options are correctly displayed in redSys.parameters
            opts.arnoldiPH.constAlg = 'phs';
            opts.arnoldiPH.structurePreservation = 'Cholesky+';
            redSys = irkaPH(sys,5,opts);

            testCase.verifyEqual(redSys.parameters.arnoldiPH.constAlg, opts.arnoldiPH.constAlg);
            testCase.verifyEqual(redSys.parameters.arnoldiPH.structurePreservation, ...
                opts.arnoldiPH.structurePreservation);
        end

        function redSysInfo(testCase, sys)
            % Verify that redSys object has additional information iter and
            % stopCritEvolution
            opts.arnoldiPH.constAlg = 'phs';
            opts.arnoldiPH.structurePreservation = 'Cholesky+';
            redSys = irkaPH(sys,5,opts);

            testCase.verifyTrue(isstruct(redSys.info));
            testCase.verifyTrue(isfield(redSys.info, "iter"));
            testCase.verifyTrue(isfield(redSys.info, "stopCritEvolution"));
            testCase.verifyEmpty(redSys.info.arnoldiPH);
        end
    end

    methods (Test, TestTags = {'Options'})
        % check if options are applied correctly

        function scaling(testCase,sys)
            % Test structure preserving option 'scaling' of arnoldiPH

            opts.arnoldiPH.constAlg = 'phs';
            opts.arnoldiPH.structurePreservation = 'scaling';
            [redSys,V] = irkaPH(sys,5,opts);
            testCase.verifyEqual(redSys.parameters.arnoldiPH.structurePreservation, 'scaling')
            testCase.verifyEqual(rank(V),  5);
            testCase.verifyLessThan(norm(eye(size(V'*V)) - V'*V), 1e-14);
            testCase.verifyEqual(redSys.dim, 5);
        end

        function specialInverse(testCase,sys)
            % Test structure preserving option 'specialInverse' of arnoldiPH

            opts.arnoldiPH.constAlg = 'phs';
            opts.arnoldiPH.structurePreservation = 'specialInverse';
            [redSys,V] = irkaPH(sys,5,opts);
            testCase.verifyEqual(redSys.parameters.arnoldiPH.structurePreservation, 'specialInverse')
            testCase.verifyEqual(rank(V),  5);
            testCase.verifyLessThan(norm(eye(size(V'*V)) - V'*V), 1e-14);
            testCase.verifyEqual(redSys.dim, 5);
        end

        function Cholesky(testCase,sys)
            % Test structure preserving option 'Cholesky' of arnoldiPH

            opts.arnoldiPH.constAlg = 'phs';
            opts.arnoldiPH.structurePreservation = 'Cholesky';
            [redSys,V] = irkaPH(sys,5,opts);
            testCase.verifyEqual(redSys.parameters.arnoldiPH.structurePreservation, 'Cholesky')
            testCase.verifyEqual(rank(V),  5);
            testCase.verifyEqual(redSys.dim, 5);
        end

        function Cholesky_plus(testCase,sys)
            % Test structure preserving option 'Cholesky+' of arnoldiPH

            opts.arnoldiPH.constAlg = 'phs';
            opts.arnoldiPH.structurePreservation = 'Cholesky+';
            [redSys,V] = irkaPH(sys,5,opts);
            testCase.verifyEqual(redSys.parameters.arnoldiPH.structurePreservation, 'Cholesky+')
            testCase.verifyEqual(rank(V),  5);
            testCase.verifyEqual(redSys.dim, 5);
        end

        function QV(testCase,sys)
            % Test structure preserving option 'QV' of arnoldiPH

            opts.arnoldiPH.constAlg = 'phs';
            opts.arnoldiPH.structurePreservation = 'QV';
            [redSys,V] = irkaPH(sys,5,opts);
            testCase.verifyEqual(redSys.parameters.arnoldiPH.structurePreservation, 'QV')
            testCase.verifyEqual(rank(V),  5);
            testCase.verifyLessThan(norm(eye(size(V'*V)) - V'*V), 1e-14);
            testCase.verifyEqual(redSys.dim, 5);
        end

        function stoppingCriteria(testCase,sys)
            % Test stopping criteria for IRKA iteration

            opts = struct();
            redSys = irkaPH(sys,5,opts);
            testCase.verifyEqual(redSys.dim, 5);
            testCase.verifyEqual(redSys.parameters.stopCrit,'s0');

            opts.stopCrit = 's0';
            redSys = irkaPH(sys,5,opts);
            testCase.verifyEqual(redSys.dim, 5);
            testCase.verifyEqual(redSys.parameters.stopCrit, 's0');

            opts.stopCrit = 'sysr';
            redSys = irkaPH(sys,5,opts);
            testCase.verifyEqual(redSys.dim, 5);
            testCase.verifyEqual(redSys.parameters.stopCrit, 'sysr');

            opts.stopCrit = 'combAll';
            redSys = irkaPH(sys,5,opts);
            testCase.verifyEqual(redSys.dim, 5);
            testCase.verifyEqual(redSys.parameters.stopCrit, 'combAll');

            opts.stopCrit = 'combAny';
            redSys = irkaPH(sys,5,opts);
            testCase.verifyEqual(redSys.dim, 5);
            testCase.verifyEqual(redSys.parameters.stopCrit, 'combAny');

            opts.stopCrit = 's0+tanDir'; % not meaningful for SISO but should work anyway
            opts.degTol = 10;
            redSys = irkaPH(sys,5,opts);
            testCase.verifyEqual(redSys.dim, 5);
            testCase.verifyEqual(redSys.parameters.stopCrit, 's0+tanDir');
            testCase.verifyEqual(redSys.parameters.degTol, opts.degTol);
        end

        function stoppingCriteriaMIMO(testCase,sysMIMO)
            % Test stopping criteria of IRKA iteration for MIMO system

            opts = struct();
            redSys = irkaPH(sysMIMO,5,opts);
            testCase.verifyEqual(redSys.dim, 5);
            testCase.verifyEqual(redSys.parameters.stopCrit, 's0');

            opts.stopCrit = 's0';
            redSys = irkaPH(sysMIMO,5,opts);
            testCase.verifyEqual(redSys.dim, 5);
            testCase.verifyEqual(redSys.parameters.stopCrit, 's0');

            opts.stopCrit = 'sysr';
            redSys = irkaPH(sysMIMO,5,opts);
            testCase.verifyEqual(redSys.dim, 5);
            testCase.verifyEqual(redSys.parameters.stopCrit, 'sysr');

            opts.stopCrit = 'combAll';
            redSys = irkaPH(sysMIMO,5,opts);
            testCase.verifyEqual(redSys.dim, 5);
            testCase.verifyEqual(redSys.parameters.stopCrit, 'combAll');

            opts.stopCrit = 'combAny';
            redSys = irkaPH(sysMIMO,5,opts);
            testCase.verifyEqual(redSys.dim, 5);
            testCase.verifyEqual(redSys.parameters.stopCrit, 'combAny');

            opts.stopCrit = 's0+tanDir';
            opts.degTol = 10;
            redSys = irkaPH(sysMIMO,5,opts);
            testCase.verifyEqual(redSys.dim, 5);
            testCase.verifyEqual(redSys.parameters.stopCrit, 's0+tanDir');
            testCase.verifyEqual(redSys.parameters.degTol, opts.degTol);

            testCase.verifyTrue(redSys.isMIMO);
        end

        function decreasedTolerance(testCase,sys)
            % Test decreasing tolerance and increasing iteration number by
            % comparing to a reference iteration result

            opts.maxIter = 100;
            [redSys_ref,V] = irkaPH(sys,5,opts); % for reference (tol = 1e-3)
            opts.tol = 1e-7;
            [redSys,V] = irkaPH(sys,5,opts);
            testCase.verifyGreaterThan(redSys.info.iter, redSys_ref.info.iter);
            testCase.verifyLessThanOrEqual(redSys.info.stopCritEvolution(end), opts.tol)
        end

        function shiftInitialization(testCase,sys)
            % Test different shift initialization strategies:
            %   - By initShifts function
            %   - With zeros
            %   - With mirrored eigenvalues
            %   - Custom: initial shifts provided
            %   - Bad input --> error

            sys = setup_MassSpringDamperSystem(100, 2, 1, 1);
            sysMIMO = setup_MassSpringDamperSystem(100, 2, 1, 1, 'MIMO');

            Opts = struct();
            redSys = irkaPH(sys, 20);
            redSysMIMO = irkaPH(sysMIMO, 20);
            testCase.verifyFalse(redSys.isMIMO)
            testCase.verifyTrue(redSysMIMO.isMIMO)
            testCase.verifyEqual(redSys.parameters.initShifts, 'eig_circle')

            strategies = {'eig_circle', 'zeros', 'linear', 'logarithmic', 'eig_large', 'eig_small', 'diag'};

            % Use initShifts
            for iStrategy = 1:length(strategies)
                Opts.initShifts = strategies{iStrategy};
                redSys = irkaPH(sys, 20, Opts);
                redSysMIMO = irkaPH(sysMIMO, 20, Opts);
                testCase.verifyFalse(redSys.isMIMO)
                testCase.verifyTrue(redSysMIMO.isMIMO)
                testCase.verifyEqual(redSys.parameters.initShifts, Opts.initShifts)
                testCase.verifyEqual(redSysMIMO.parameters.initShifts, Opts.initShifts)
            end

            % Provide initial shifts
            shifts = zeros(20,1);
            for iStrategy = 1:length(strategies)
                Opts.initShifts = strategies{iStrategy};
                redSys = irkaPH(sys, shifts);
                redSysMIMO = irkaPH(sysMIMO, shifts);
                testCase.verifyFalse(redSys.isMIMO)
                testCase.verifyTrue(redSysMIMO.isMIMO)
                testCase.verifyEqual(redSys.parameters.initShifts, 'custom')
                testCase.verifyEqual(redSysMIMO.parameters.initShifts, 'custom')
                testCase.verifyEqual(redSys.parameters.startShifts, shifts)
                testCase.verifyEqual(redSysMIMO.parameters.startShifts, shifts)
            end

            % Provide initial shifts and tangent directions
            shifts = rand(20,1);
            b = 2*ones(2, 20);
            for iStrategy = 1:length(strategies)
                redSysMIMO = irkaPH(sysMIMO, shifts, b);
                testCase.verifyTrue(redSysMIMO.isMIMO)
                testCase.verifyEqual(redSysMIMO.parameters.initShifts, 'custom')
                testCase.verifyEqual(redSysMIMO.parameters.startShifts, shifts)
                testCase.verifyEqual(redSysMIMO.parameters.startTangent, b)
            end

            Opts.initShifts = 'unknownStrategy';
            testCase.verifyError(@() irkaPH(sys, 20, Opts),...
                'MORpH:phsMOR_parseOpts:input_not_admissible');
        end

    end

    methods (Test, TestTags = {'Options','verbose'})
        % Test verbose mode
        function verboseOff(test)
            % Test verbose mode
            % Verify correct behavior by looking at the command window output!

            warning('on','all')
            sys = setup_MassSpringDamperSystem(40,2,1,1);

            opts.verbose = false;
            disp("Starting non-verbose mode");
            [redSys,V] = irkaPH(sys,5,opts);
            disp("Finished computation");
            warning('off','all')
        end
    end

    methods (Test, TestTags = {'Functionality'})
        % Verify correct outputs for more realistic scenarios

        function accuracy(testCase)
            % Check accuracy for order 200 model

            sys = setup_MassSpringDamperSystem(200,2,1,1);

            redSys = irkaPH(sys,20);
            testCase.verifyLessThan(norm(redSys-sys)/norm(sys), 1e-3);

            redSys = irkaPH(sys,50);
            testCase.verifyLessThan(norm(redSys-sys)/norm(sys), 1e-6);

            redSys = irkaPH(sys,100);
            testCase.verifyLessThan(norm(redSys-sys)/norm(sys), 1e-8);
        end
    end

end