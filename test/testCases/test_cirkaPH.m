classdef test_cirkaPH < matlab.unittest.TestCase
    % TEST test_cirkaPH - Tests the functionalities of cirkaPH
    %
    % Press <F5> or enter "runtests("test_cirkaPH")" to run this testscript
    %
    % Description:
    %   This script tests the function 'cirkaPH' of the MORpH-toolbox. The
    %   following test scenarios are considered:
    %       - Input: Check if function can run all intended input combinations
    %       - Output: Check if each additional output can be obtained
    %       - Options: Check if additional options work correctly
    %       - Functionality: Try more complex reduction to simulate real-world
    %                        application
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

    properties (TestParameter)
        sys = {setup_MassSpringDamperSystem(40,2,1,1)};   % standard MSD system for most of the tests
        sysMIMO = {setup_MassSpringDamperSystem(40,2,1,1,'MIMO')}; % MIMO MSD system for MIMO testing
        shifts = {logspace(-5,5,5)};
    end

    properties
        tolerance = 1e-7;
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
        % Test input combinations

        function standardInput(test,sys)
            % Very basic input: system and reduced order only

            redSys = cirkaPH(sys,5);
            test.verifyEqual(redSys.dim, 5);
            test.verifyEqual(redSys.method, @cirkaPH);
        end

        function standardInputWithOpts(testCase,sys)
            % Pass opts struct
            opts.tol = 1e-5;
            opts.irkaPH.tol = 1e-6;
            opts.irkaPH.arnoldiPH.structurePreservation = 'QV';
            redSys = cirkaPH(sys,5,opts);
            testCase.verifyEqual(redSys.dim, 5);
            testCase.verifyEqual(redSys.parameters.tol,opts.tol)
            testCase.verifyEqual(redSys.parameters.irkaPH.tol,opts.irkaPH.tol)
            testCase.verifyEqual(redSys.parameters.irkaPH.arnoldiPH.structurePreservation,...
                opts.irkaPH.arnoldiPH.structurePreservation)
        end

        function initialShifts(testCase,sys,shifts)
            % Supply initial shifts for first iteration instead of reduced
            % order

            redSys = cirkaPH(sys,shifts);
            testCase.verifyEqual(redSys.dim, length(shifts));
            testCase.verifyEqual(redSys.parameters.startShifts, shifts)
        end

        function MIMO(testCase,sysMIMO)
            % Test MIMO system reduction (with opts)

            opts.irkaPH.arnoldiPH.constAlg = 'phs';
            redSys = cirkaPH(sysMIMO,20,opts);
            testCase.verifyEqual(redSys.dim, 20);
            testCase.verifyEqual(redSys.parameters.irkaPH.arnoldiPH.constAlg,...
                opts.irkaPH.arnoldiPH.constAlg);
        end

        function MIMOWithInitialShifts(testCase,sysMIMO)
            % Supply initial shifts for first iteration (MIMO reduction)

            shifts = logspace(-10,10,20);
            opts.irkaph.constAlg = 'phs';
            redSys = cirkaPH(sysMIMO,shifts,opts);
            testCase.verifyEqual(redSys.dim, 20);
            %             testCase.verifyEqual(redSys.parameters.startShifts,shifts)
        end

        function MIMOWithInitialShiftsAndTangent(testCase,sysMIMO)
            % Supply initial shifts and tangent directions

            shifts = logspace(-10,10,20);
            tangent = ones(size(sysMIMO.G,2),length(shifts));
            opts.irkaph.constAlg = 'phs';
            redSys = cirkaPH(sysMIMO,shifts,tangent,opts);
            testCase.verifyEqual(redSys.dim, 20);
            testCase.verifyEqual(redSys.parameters.startShifts,shifts)
            testCase.verifyEqual(redSys.parameters.startTangent,tangent)
        end
    end

    methods (Test, TestTags = {'Output'})
        % Check if optional output parameters are provided correctly

        function standardInputWithOptionalOutput(testCase,sys)
            % Provide standard input and request V and s0
            % V should be a sys.dim x 5 matrix
            % s0 should be a 1 x 5 row vector

            [redSys, V, W, s0, Rt, sysm, s0mTot, nLU] = cirkaPH(sys,5);
            testCase.verifyEqual(redSys.dim, 5);
            testCase.verifyEqual(rank(V), 5);
            testCase.verifySize(V,[sys.dim,5]);
            testCase.verifySize(s0,[5,1]);
            redSys.parameters.irkaPH.arnoldiPH
            redSys.parameters.irkaPH.arnoldiPH.structurePreservation
            redSys.parameters.modelFctPH.structurePreservation
            J = W'*sys.J*W;
            Q = V'*sys.Q*V;
            testCase.verifyLessThan(norm(redSys.J-J), 1e-12);
            testCase.verifyLessThan(norm(redSys.Q-Q), 1e-12);
        end

        function outputObjectParameters(testCase, sys)
            % Check if the output object redSys provides additional information
            % in struct parameters (see phsRed class)

            redSys = cirkaPH(sys,20);

            % property "parameters" of the returned object should have the
            % following variables:
            testCase.verifyTrue(isfield(redSys.parameters,'startShifts'));
            testCase.verifyTrue(isfield(redSys.parameters,'startTangent'));
            testCase.verifyTrue(isfield(redSys.parameters,'s0m'));
            testCase.verifyTrue(isfield(redSys.parameters,'Rtm'));
            testCase.verifyTrue(isfield(redSys.parameters,'tol'));
            testCase.verifyTrue(isfield(redSys.parameters,'maxIter'));
            testCase.verifyTrue(isfield(redSys.parameters,'degTol'));
            testCase.verifyTrue(isfield(redSys.parameters,'stopCrit'));
            testCase.verifyTrue(isfield(redSys.parameters,'clearInit'));
            testCase.verifyTrue(isfield(redSys.parameters,'irkaPH'));
            testCase.verifyTrue(isfield(redSys.parameters,'modelFctPH'));
            testCase.verifyTrue(isstruct(redSys.parameters.irkaPH));
            testCase.verifyTrue(isstruct(redSys.parameters.modelFctPH));
            testCase.verifyTrue(isfield(redSys.parameters.modelFctPH,'updateModel'));

            testCase.verifyTrue(isfield(redSys.info,'s0'));
            testCase.verifyTrue(isfield(redSys.info,'Rt'));
            testCase.verifyTrue(isfield(redSys.info,'s0m'));
            testCase.verifyTrue(isfield(redSys.info,'Rtm'));
            testCase.verifyTrue(isfield(redSys.info,'kIrka'));
            testCase.verifyTrue(isfield(redSys.info,'finalStopCrit'));
            testCase.verifyTrue(isfield(redSys.info,'relH2err'));
            testCase.verifyTrue(isfield(redSys.info,'nHighDimLU'));
            testCase.verifyTrue(isfield(redSys.info,'modelFctOrder'));
            testCase.verifyTrue(isfield(redSys.info,'originalOrder'));
            testCase.verifyTrue(isfield(redSys.info,'irkaPH'));
            testCase.verifyTrue(isstruct(redSys.info.irkaPH));
        end

        function redSysMethod(testCase, sys)
            redSys = cirkaPH(sys,5);
            testCase.verifyEqual(redSys.method, @cirkaPH)
        end

        function redSysInfo(testCase, sys)
            redSys = cirkaPH(sys,5);
            testCase.verifyTrue(isfield(redSys.info, 's0'));
            testCase.verifyTrue(isfield(redSys.info, 'Rt'));
            testCase.verifyTrue(isfield(redSys.info, 's0m'));
            testCase.verifyTrue(isfield(redSys.info, 'Rtm'));
            testCase.verifyTrue(isfield(redSys.info, 'finalStopCrit'));
            testCase.verifyTrue(isfield(redSys.info, 'kIrka'));
            testCase.verifyTrue(isfield(redSys.info, 'nHighDimLU'));
            testCase.verifyTrue(isfield(redSys.info, 'originalOrder'));
            testCase.verifyTrue(isfield(redSys.info, 'modelFctOrder'));
            testCase.verifyTrue(isfield(redSys.info, 'relH2err'));
            testCase.verifyTrue(isfield(redSys.info, 'irkaPH'));
            testCase.verifyTrue(isstruct(redSys.info.irkaPH));
        end
    end

    methods (Test, TestTags = {'Options'})
        % Test optionalal computation arguments (opts struct)

        function structurePreservation(test,sys)
            % Check arnoldiPH argument structurePreservation

            opts.irkaPH.arnoldiPH.constAlg = 'phs';
            opts.modelFctPH.structurePreservation = 'specialInverse';
            opts.irkaPH.arnoldiPH.structurePreservation = 'scaling';
            [redSys,V] = cirkaPH(sys,5,opts);
            test.verifyEqual(rank(V), 5);
            test.verifyLessThan(norm(eye(size(V'*V)) - V'*V), 1e-14);
            test.verifyEqual(redSys.dim, 5);
            test.verifyEqual(redSys.parameters.irkaPH.arnoldiPH.constAlg, opts.irkaPH.arnoldiPH.constAlg);
            test.verifyEqual(redSys.parameters.modelFctPH.structurePreservation, opts.modelFctPH.structurePreservation);
            test.verifyEqual(redSys.parameters.irkaPH.arnoldiPH.structurePreservation, opts.irkaPH.arnoldiPH.structurePreservation);

            opts.modelFctPH.structurePreservation = 'specialInverse';
            opts.irkaPH.arnoldiPH.structurePreservation = 'specialInverse';
            [redSys,V] = cirkaPH(sys,5,opts);
            test.verifyEqual(rank(V), 5);
            test.verifyLessThan(norm(eye(size(V'*V)) - V'*V), 1e-14);
            test.verifyEqual(redSys.dim, 5);
            test.verifyEqual(redSys.parameters.irkaPH.arnoldiPH.constAlg, opts.irkaPH.arnoldiPH.constAlg);
            test.verifyEqual(redSys.parameters.modelFctPH.structurePreservation, opts.modelFctPH.structurePreservation);
            test.verifyEqual(redSys.parameters.irkaPH.arnoldiPH.structurePreservation, opts.irkaPH.arnoldiPH.structurePreservation);

            opts.modelFctPH.structurePreservation = 'Cholesky';
            opts.irkaPH.arnoldiPH.structurePreservation = 'Cholesky';
            [redSys,V] = cirkaPH(sys,5,opts);
            test.verifyEqual(rank(V), 5);
            test.verifyEqual(redSys.dim, 5);
            test.verifyEqual(redSys.parameters.irkaPH.arnoldiPH.constAlg, opts.irkaPH.arnoldiPH.constAlg);
            test.verifyEqual(redSys.parameters.modelFctPH.structurePreservation, opts.modelFctPH.structurePreservation);
            test.verifyEqual(redSys.parameters.irkaPH.arnoldiPH.structurePreservation, opts.irkaPH.arnoldiPH.structurePreservation);

            opts.modelFctPH.structurePreservation = 'Cholesky+';
            opts.irkaPH.arnoldiPH.structurePreservation = 'Cholesky+';
            [redSys,V] = cirkaPH(sys,5,opts);
            test.verifyEqual(rank(V), 5);
            test.verifyEqual(redSys.dim, 5);
            test.verifyEqual(redSys.parameters.irkaPH.arnoldiPH.constAlg, opts.irkaPH.arnoldiPH.constAlg);
            test.verifyEqual(redSys.parameters.modelFctPH.structurePreservation, opts.modelFctPH.structurePreservation);
            test.verifyEqual(redSys.parameters.irkaPH.arnoldiPH.structurePreservation, opts.irkaPH.arnoldiPH.structurePreservation);

            opts.modelFctPH.structurePreservation = 'QV';
            opts.irkaPH.arnoldiPH.structurePreservation = 'QV';
            [redSys,V] = cirkaPH(sys,5,opts);
            test.verifyEqual(rank(V), 5);
            test.verifyLessThan(norm(eye(size(V'*V)) - V'*V), 1e-14);
            test.verifyEqual(redSys.dim, 5);
            test.verifyEqual(redSys.parameters.irkaPH.arnoldiPH.constAlg, opts.irkaPH.arnoldiPH.constAlg);
            test.verifyEqual(redSys.parameters.modelFctPH.structurePreservation, opts.modelFctPH.structurePreservation);
            test.verifyEqual(redSys.parameters.irkaPH.arnoldiPH.structurePreservation, opts.irkaPH.arnoldiPH.structurePreservation);
        end

        function shiftInitialization(testCase, sys)
            Opts.initShifts = 'zeros';
            testSys = cirkaPH(sys, 5, Opts);
            testCase.verifyEqual(testSys.parameters.initShifts, Opts.initShifts);

            Opts.initShifts = 'eig_circle';
            testSys = cirkaPH(sys, 5, Opts);
            testCase.verifyEqual(testSys.parameters.initShifts, Opts.initShifts);

            Opts.initShifts = 'linear';
            testSys = cirkaPH(sys, 5, Opts);
            testCase.verifyEqual(testSys.parameters.initShifts, Opts.initShifts);

            Opts.initShifts = 'diag';
            testSys = cirkaPH(sys, 5, Opts);
            testCase.verifyEqual(testSys.parameters.initShifts, Opts.initShifts);
        end

        function stoppingCriteria(testCase,sys)
            % Test stopping criteria of cirkaPH

            % default
            opts = struct();
            redSys = cirkaPH(sys,5,opts);
            testCase.verifyEqual(redSys.dim, 5);
            testCase.verifyEqual(redSys.parameters.stopCrit,'combAny');

            opts.stopCrit = 's0';
            redSys = cirkaPH(sys,5,opts);
            testCase.verifyEqual(redSys.dim, 5);
            testCase.verifyEqual(redSys.parameters.stopCrit, 's0');

            opts.stopCrit = 'sysr';
            redSys = cirkaPH(sys,5,opts);
            testCase.verifyEqual(redSys.dim, 5);
            testCase.verifyEqual(redSys.parameters.stopCrit, 'sysr');

            opts.stopCrit = 'combAll';
            redSys = cirkaPH(sys,5,opts);
            testCase.verifyEqual(redSys.dim, 5);
            testCase.verifyEqual(redSys.parameters.stopCrit, 'combAll');

            opts.stopCrit = 'combAny';
            redSys = cirkaPH(sys,5,opts);
            testCase.verifyEqual(redSys.dim, 5);
            testCase.verifyEqual(redSys.parameters.stopCrit, 'combAny');

            opts.stopCrit = 's0+tanDir'; % not meaningful for SISO but should work anyway
            opts.degTol = 10;
            redSys = cirkaPH(sys,5,opts);
            testCase.verifyEqual(redSys.dim, 5);
            testCase.verifyEqual(redSys.parameters.stopCrit, 's0+tanDir');
            testCase.verifyEqual(redSys.parameters.degTol, opts.degTol);
        end

        function stoppingCriteriaMIMO(testCase,sysMIMO)
            % Test stopping criteria of cirkaPH (MIMO case)

            % default
            opts = struct();
            redSys = cirkaPH(sysMIMO,5,opts);
            testCase.verifyEqual(redSys.dim, 5);
            testCase.verifyEqual(redSys.parameters.stopCrit,'combAny');

            opts.stopCrit = 's0';
            redSys = cirkaPH(sysMIMO,5,opts);
            testCase.verifyEqual(redSys.dim, 5);
            testCase.verifyEqual(redSys.parameters.stopCrit, 's0');

            opts.stopCrit = 'sysr';
            redSys = cirkaPH(sysMIMO,5,opts);
            testCase.verifyEqual(redSys.dim, 5);
            testCase.verifyEqual(redSys.parameters.stopCrit, 'sysr');

            opts.stopCrit = 'combAll';
            redSys = cirkaPH(sysMIMO,5,opts);
            testCase.verifyEqual(redSys.dim, 5);
            testCase.verifyEqual(redSys.parameters.stopCrit, 'combAll');

            opts.stopCrit = 'combAny';
            redSys = cirkaPH(sysMIMO,5,opts);
            testCase.verifyEqual(redSys.dim, 5);
            testCase.verifyEqual(redSys.parameters.stopCrit, 'combAny');

            opts.stopCrit = 's0+tanDir';
            opts.degTol = 10;
            redSys = cirkaPH(sysMIMO,5,opts);
            testCase.verifyEqual(redSys.dim, 5);
            testCase.verifyEqual(redSys.parameters.stopCrit, 's0+tanDir');
            testCase.verifyEqual(redSys.parameters.degTol, opts.degTol);
        end

        function restrictMaxIter(testCase, sys)
            % Test option maxIter

            opts.maxIter = 2;
            redSys = cirkaPH(sys,5,opts);
            testCase.verifyLessThanOrEqual(length(redSys.info.kIrka),opts.maxIter);
        end

        function restrictIrkaMaxIter(testCase, sys)
            % Restrict irkaPH maxIter --> should result in much more cirkaPH
            % outer iterations
            Opts.irkaPH.maxIter = 2;
            redSys_default = cirkaPH(sys,5);    % no restriction for comparison
            redSys = cirkaPH(sys,5,Opts);       % restriction
            testCase.verifyGreaterThan(length(redSys.info.kIrka), length(redSys_default.info.kIrka));
            testCase.verifyTrue(all(redSys.info.kIrka <= 2));
            testCase.verifyGreaterThan(redSys.info.modelFctOrder,redSys_default.info.modelFctOrder);
        end

        function decreaseTolerance(testCase, sys)
            % Test option tol by decreasing it and verifying result

            opts.maxIter = 200; % make sure to allow enough iterations to converge
            opts.tol = 1e-6; % default is 1e-3
            redSys = cirkaPH(sys,20,opts);
            testCase.verifyLessThanOrEqual(min(redSys.info.finalStopCrit),opts.tol)
        end

    end

    methods (Test, TestTags = {'Options','verbose'})
        function verboseOff(test)
            sys = setup_MassSpringDamperSystem(40,2,1,1);

            opts.verbose = false;
            opts.suppressWarn = true;
            warning('on','all')
            disp(" --- Starting non-verbose mode --- ");
            [redSys,V] = cirkaPH(sys,5,opts);
            disp(" --- Finished computation --- ");
            warning('off','all')
        end
        function maximumConsoleOutput(test)
            sys = setup_MassSpringDamperSystem(40,2,1,1);

            opts.verbose = true;
            opts.suppressWarn = false;
            opts.irkaPH.summary = true;
            opts.irkaPH.verbose = true;
            warning('on','all')
            disp(" --- Starting all-output mode --- ");
            [redSys,V] = cirkaPH(sys,5,opts);
            disp(" --- Finished computation --- ");
            warning('off','all')
        end
    end

    methods (Test, TestTags = {'Functionality'})
        % Simulate more realistic scenarios with larger original models
        % Compare results with other functions (irkaPH)

        function accuracy(testCase)
            % make sure reduced system approximates the original
            sys = setup_MassSpringDamperSystem(200,2,1,1); % order 200
            redSys = cirkaPH(sys,30);
            testCase.verifyLessThanOrEqual(norm(redSys-sys)/norm(sys), 1e-4);
        end

        function cirkaPHvsIrkaPH(testCase)
            % Compare results of cirkaPH to those of irkaPH

            sys = setup_MassSpringDamperSystem(200,2,1,1); % order 200

            redSys_irka = irkaPH(sys,30);
            redSys_cirka = cirkaPH(sys,30);
            err_irka = norm(redSys_irka-sys)/norm(sys);
            err_cirka = norm(redSys_cirka-sys)/norm(sys);
            testCase.verifyLessThan(abs(err_irka - err_cirka), testCase.tolerance);
        end

        function cirkaPHvsIrkaPHMIMO(testCase)
            % Compare results of cirkaPH to those of irkaPH (MIMO case)

            sysMIMO = setup_MassSpringDamperSystem(300,2,1,1,'MIMO'); % order 300

            opts.irkaph.constAlg = 'phs';
            opts.constAlg = 'phs';
            opts.tol = 1e-10;
            opts.irkaph.tol = opts.tol;
            redSys_irka = irkaPH(sysMIMO,30,opts);
            redSys_cirka = cirkaPH(sysMIMO,30,opts);
            err_irka = norm(redSys_irka-sysMIMO)/norm(sysMIMO);
            err_cirka = norm(redSys_cirka-sysMIMO)/norm(sysMIMO);
            testCase.verifyLessThanOrEqual(abs(err_irka-err_cirka), testCase.tolerance);
        end

    end

    methods (Test, TestTags = {'Input','Options'})

        function surrogateInitialization(testCase, sys, shifts)
            % Test reduction with custom initial surrogate model shifts

            shifts_surrogate = [shifts,shifts/2];
            opts.s0m = shifts_surrogate;
            redSys = cirkaPH(sys,shifts,opts);
            testCase.verifyEqual(redSys.parameters.s0m, shifts_surrogate);
            testCase.verifyNotEqual(redSys.info.s0m, shifts_surrogate);
        end

        function surrogateInitializationAndTangent(test, sysMIMO, shifts)
            % Test reduction with custom initial surrogate model shifts and
            % tangent directions

            shifts_surrogate = [shifts,shifts/2];
            %             tangent = rand(size(sysMIMO.G,2),length(shifts));
            tangent = ones(size(sysMIMO.G,2),length(shifts));
            tangent_surrogate = [tangent, 2*tangent];
            opts.irkaph.constAlg = 'phs';
            opts.s0m = shifts_surrogate;
            opts.Rtm = tangent_surrogate;
            redSys = cirkaPH(sysMIMO,shifts,opts);
            test.verifyEqual(redSys.parameters.s0m, shifts_surrogate);
        end

    end

    methods (Test, TestTags = {'Auxiliary Functions'})
        function modelFctPHTest(testCase,sys)
            shifts1 = [1e-3, 1e-4, 1e-1, 1e1, 1e3];
            shifts2 = [1,2,3,4,5];
            shifts3 = [shifts1,shifts2];
            shifts4 = [1, 2, 1e-3, 10, 20, 30];
            Opts.structurePreservation = 'specialInverse';

            [sysm1,s0mTot1,RtmTot1,V1,W1,nLU1] = modelFctPH(sys,shifts1,Opts);
            [sysm2,s0mTot2,RtmTot2,V2,W2,nLU2] = modelFctPH(sys,shifts2,s0mTot1,V1,W1,Opts);
            [sysm3,s0mTot3,RtmTot3,V3,W3,nLU3] = modelFctPH(sys,shifts3,s0mTot2,V2,W2,Opts);
            [sysm4,s0mTot4,RtmTot4,V4,W4,nLU4] = modelFctPH(sys,shifts4,s0mTot2,V2,W2,Opts);

            % Check reduced system dimension
            testCase.verifyEqual(sysm1.dim, length(shifts1));
            testCase.verifyEqual(sysm2.dim, length(shifts1) + length(shifts2));
            testCase.verifyEqual(sysm3.dim, length(shifts1) + length(shifts2));
            testCase.verifyEqual(sysm4.dim, length(unique([shifts1,shifts2,shifts4])));
            % Check equivalence of shift vectors
            testCase.verifyEqual(s0mTot1, shifts1);
            testCase.verifyEqual(s0mTot2, [shifts1,shifts2]);
            testCase.verifyEqual(s0mTot3, shifts3);
            testCase.verifyEqual(sort(s0mTot4), unique([shifts1,shifts2,shifts4]));
            % Check returned matrices Vi
            %             testCase.verifyEqual(V3,V2,'AbsTol',1e-14,'RelTol',1e-14);
            testCase.verifyEqual(rank(V1*V1'*V2),rank(V1));
            testCase.verifyEqual(rank(V2*V2'*V3),rank(V2));
            testCase.verifyEqual(rank(V2*V2'*V4),rank(V2));
            % Check number of nLU decompositions
            testCase.verifyNotEqual(nLU1,0);
            testCase.verifyNotEqual(nLU2,0);
            testCase.verifyEqual(nLU3, 0);
            testCase.verifyEqual(nLU4, length(unique([shifts1,shifts2,shifts4])) - length(s0mTot2));
            % Verify that reduced systems are PH
            phs.inputValidation(sysm1);
            phs.inputValidation(sysm2);
            phs.inputValidation(sysm3);
        end
    end
end