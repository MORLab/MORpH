classdef test_arnoldiPH < matlab.unittest.TestCase
    % TEST test_arnoldiPH - Tests the functionalities of arnoldiPH
    %
    % Press <F5> or enter "runtests("test_arnoldiPH")" to run this testscript
    %
    % Description:
    %   This script tests the function 'arnoldiPH' of the MORpH-toolbox. The
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

    %% Properties
    properties (TestParameter)
        sys = {setup_MassSpringDamperSystem(40,2,1,1)};   % standard MSD system for most of the tests
        sysMIMO = {setup_MassSpringDamperSystem(40,2,1,1,'MIMO')}; % MIMO MSD system for MIMO testing
        shifts = {logspace(-10,10,5)};  % standard shifts
    end

    properties
        tolerance = 5e-8;   % default tolerance
    end

    %% Fixture functions (setup test environment)
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

    %% Test methods
    methods (Test, TestTags = {'Input'})
        % check input combinations

        function standardInput(testCase, sys, shifts)
            % Verify that standard input (system + shifts) works
            % Verify that moments are matched at shifts

            redSys = arnoldiPH(sys,shifts);
            testCase.verifyEqual(redSys.dim, 5);
            % accuracy of transfer function at shifts
            sys_ss = ss(sys);
            redSys_ss = ss(redSys);
            G_original = @(s) sys_ss.C*((s*eye(size(sys_ss.A)) - sys_ss.A)\sys_ss.B);
            G_reduced = @(s) redSys_ss.C*((s*eye(size(redSys_ss.A)) - redSys_ss.A)\redSys_ss.B);
            for i = 1:length(shifts)
                delta = abs(G_original(shifts(i)) - G_reduced(shifts(i)));
                testCase.verifyLessThan(delta, 1e-10);
            end
        end

        function standardInputWithOpts(testCase, sys, shifts)
            % Pass additional opts-struct

            opts.constAlg = 'phs';
            opts.structurePreservation = 'specialInverse';
            opts.verbose = false;
            redSys = arnoldiPH(sys,shifts,opts);
            testCase.verifyEqual(redSys.dim, 5);
        end

        function imaginaryShifts(testCase, sys)
            % Test reduction at imaginary shifts
            % Verify that moments are matched at shifts

            shifts = 1i*logspace(-10,10,5);
            opts.constAlg = 'phs';
            redSys = arnoldiPH(sys,shifts,opts);
            testCase.verifyEqual(redSys.dim, 10)
            % accuracy of transfer function at imaginary shifts
            sys_ss = ss(sys);
            redSys_ss = ss(redSys);
            G_original = @(s) sys_ss.C*((s*eye(size(sys_ss.A)) - sys_ss.A)\sys_ss.B);
            G_reduced = @(s) redSys_ss.C*((s*eye(size(redSys_ss.A)) - redSys_ss.A)\redSys_ss.B);
            for i = 1:length(shifts)
                delta = abs(G_original(shifts(i)) - G_reduced(shifts(i)));
                testCase.verifyLessThan(delta, 1e-10);
            end
        end

        function complexShifts(testCase, sys)
            % Test reduction at complex shifts
            % Verify that moments are matched at shifts

            shifts = 1i*logspace(-10,10,5) + logspace(-10,10,5);
            opts.constAlg = 'phs';
            redSys = arnoldiPH(sys,shifts,opts);
            testCase.verifyEqual(redSys.dim, 10)
            % accuracy of transfer function at complex shifts
            sys_ss = ss(sys);
            redSys_ss = ss(redSys);
            G_original = @(s) sys_ss.C*((s*eye(size(sys_ss.A)) - sys_ss.A)\sys_ss.B);
            G_reduced = @(s) redSys_ss.C*((s*eye(size(redSys_ss.A)) - redSys_ss.A)\redSys_ss.B);
            for i = 1:length(shifts)
                delta = abs(G_original(shifts(i)) - G_reduced(shifts(i)));
                testCase.verifyLessThan(delta, 1e-10);
            end
        end

        function shiftsAtInfAndZero(testCase, sys)
            % Test reductionw with shifts at zero and infinity (Markov parameters)

            shifts = [zeros(1,3),Inf(1,2)];
            opts.constAlg = 'phs';

            redSys = arnoldiPH(sys,shifts,opts);
            testCase.verifyEqual(redSys.dim, 5);
        end

        function MIMOWithTangent(testCase, sysMIMO, shifts)
            % Test MIMO reduction with supplied tangent directions
            % Verify that moments are matched at shifts + tangent
            tangent = rand(2,5);
            opts.constAlg = 'phs';

            redSys = arnoldiPH(sysMIMO,shifts,tangent,opts);
            testCase.verifyEqual(redSys.dim, 5);

            % accuracy of transfer function at shifts
            sys_ss = ss(sysMIMO);
            redSys_ss = ss(redSys);
            G_original = @(s) sys_ss.C*((s*eye(size(sys_ss.A)) - sys_ss.A)\sys_ss.B);
            G_reduced = @(s) redSys_ss.C*((s*eye(size(redSys_ss.A)) - redSys_ss.A)\redSys_ss.B);

            for i = 1:length(shifts)
                delta = norm(G_original(shifts(i))*tangent(:,i) - G_reduced(shifts(i))*tangent(:,i));
                testCase.verifyLessThan(delta, testCase.tolerance);
            end
        end

        function MIMOWithTangentWithoutOpts(testCase, sysMIMO, shifts)
            % Test MIMO reduction with supplied tangent directions but without opts
            % Verify that moments are matched at shifts + tangent

            tangent = rand(2,5);

            redSys = arnoldiPH(sysMIMO,shifts,tangent);
            testCase.verifyEqual(redSys.dim, 5);

            % accuracy of transfer function at shifts
            sys_ss = ss(sysMIMO);
            redSys_ss = ss(redSys);
            G_original = @(s) sys_ss.C*((s*eye(size(sys_ss.A)) - sys_ss.A)\sys_ss.B);
            G_reduced = @(s) redSys_ss.C*((s*eye(size(redSys_ss.A)) - redSys_ss.A)\redSys_ss.B);

            for i = 1:length(shifts)
                delta = norm(G_original(shifts(i))*tangent(:,i) - G_reduced(shifts(i))*tangent(:,i));
                testCase.verifyLessThan(delta, testCase.tolerance);
            end
        end

        function MIMOWithComplexTangent(testCase, sysMIMO, shifts)
            % Test MIMO reduction with complex tangent directions
            % Note that 'sss' does not support this case!
            %   --> use opts.constAlg = 'phs'

            tangent = rand(2,5) + 1i*rand(2,5);
            opts.constAlg = 'phs';

            redSys = arnoldiPH(sysMIMO,shifts,tangent,opts);
            testCase.verifyEqual(redSys.dim,  10);

            % accuracy of transfer function at shifts
            sys_ss = ss(sysMIMO);
            redSys_ss = ss(redSys);
            G_original = @(s) sys_ss.C*((s*eye(size(sys_ss.A)) - sys_ss.A)\sys_ss.B);
            G_reduced = @(s) redSys_ss.C*((s*eye(size(redSys_ss.A)) - redSys_ss.A)\redSys_ss.B);

            for i = 1:length(shifts)
                delta = norm(G_original(shifts(i))*tangent(:,i) - G_reduced(shifts(i))*tangent(:,i));
                testCase.verifyLessThan(delta, testCase.tolerance);
            end
        end

        function MIMOWithComplexShiftsAndTangent(testCase, sysMIMO, shifts)
            % Test MIMO reduction with complex shifts and tangent directions
            % Note that 'sss' does not support this case!
            %   --> use opts.constAlg = 'phs'

            tangent = rand(2,5) + 1i*rand(2,5);
            shifts = 1i*logspace(-10,10,5) + logspace(-10,10,5);
            opts.constAlg = 'phs';

            redSys = arnoldiPH(sysMIMO,shifts,tangent,opts);
            testCase.verifyEqual(redSys.dim,  10);

            % accuracy of transfer function at shifts
            sys_ss = ss(sysMIMO);
            redSys_ss = ss(redSys);
            G_original = @(s) sys_ss.C*((s*eye(size(sys_ss.A)) - sys_ss.A)\sys_ss.B);
            G_reduced = @(s) redSys_ss.C*((s*eye(size(redSys_ss.A)) - redSys_ss.A)\redSys_ss.B);

            for i = 1:length(shifts)
                delta = norm(G_original(shifts(i))*tangent(:,i) - G_reduced(shifts(i))*tangent(:,i));
                testCase.verifyLessThan(delta, testCase.tolerance);
                % Complex conjugates should also work
                delta = norm(G_original(conj(shifts(i)))*conj(tangent(:,i)) - G_reduced(conj(shifts(i)))*conj(tangent(:,i)));
                testCase.verifyLessThan(delta, testCase.tolerance);
            end
        end

        function MIMOWithComplexTangentRepetitive(testCase, sysMIMO)
            % Test MIMO reduction with repeated shifts/tangent directions for
            % higher order moment matching

            shifts = repmat(logspace(-10,10,5) + 1i*logspace(-10,10,5),1,3);
            tangent = repmat(rand(2,5) + 1i*rand(2,5),1,3);
            opts.constAlg = 'phs';

            redSys = arnoldiPH(sysMIMO,shifts,tangent,opts);
            testCase.verifyEqual(redSys.dim,  30);
        end

        function MIMOWithoutTangent(testCase, sysMIMO, shifts)
            % Test MIMO reduction without provided tangent directions
            % --> Block Krylov approach
            % Verify accuracy at shifts

            opts.constAlg = 'phs';

            redSys = arnoldiPH(sysMIMO,shifts,opts);
            testCase.verifyEqual(redSys.dim,  10);

            % accuracy of transfer function at shifts
            sys_ss = ss(sysMIMO);
            redSys_ss = ss(redSys);
            G_original = @(s) sys_ss.C*((s*eye(size(sys_ss.A)) - sys_ss.A)\sys_ss.B);
            G_reduced = @(s) redSys_ss.C*((s*eye(size(redSys_ss.A)) - redSys_ss.A)\redSys_ss.B);
            for i = 1:length(shifts)
                delta = norm(G_original(shifts(i)) - G_reduced(shifts(i)));
                testCase.verifyLessThan(delta, 1e-10);
            end
        end

        function MIMOWithoutTangentWithoutOpts(testCase,sysMIMO)
            % Test MIMO reduction without provided tangent directions and
            % without opts
            % --> Block Krylov approach
            % Verify accuracy at shifts

            shifts = logspace(-10,10,5);

            redSys = arnoldiPH(sysMIMO,shifts);
            testCase.verifyEqual(redSys.dim,  10);

            % accuracy of transfer function at shifts
            sys_ss = ss(sysMIMO);
            redSys_ss = ss(redSys);
            G_original = @(s) sys_ss.C*((s*eye(size(sys_ss.A)) - sys_ss.A)\sys_ss.B);
            G_reduced = @(s) redSys_ss.C*((s*eye(size(redSys_ss.A)) - redSys_ss.A)\redSys_ss.B);
            for i = 1:length(shifts)
                delta = norm(G_original(shifts(i)) - G_reduced(shifts(i)));
                testCase.verifyLessThan(delta, 1e-10);
            end
        end

    end

    methods (Test, TestTags = {'Output'})
        % check if output is correct

        function standardInputWithOptionalOutput(testCase, sys, shifts)
            % Check that optional output is provided

            [redSys, V] = arnoldiPH(sys,shifts);
            testCase.verifyEqual(redSys.dim, 5);
            testCase.verifyEqual(rank(V), 5);
            [redSys, V, W] = arnoldiPH(sys,shifts);
            testCase.verifyEqual(redSys.dim, 5);
            testCase.verifyEqual(rank(W), 5);
            testCase.verifyEqual(size(W), size(V));
            [~, ~, ~, nLU] = arnoldiPH(sys, shifts);
            testCase.verifyEqual(nLU, length(shifts));
        end

        function nLU(testCase, sys)
            % Check that optional output nLU is provided and correct

            shifts1 = logspace(-10,10,5);
            shifts2 = zeros(5,1);
            shifts3 = [1i*shifts1, -1i*shifts1];

            opts.constAlg = 'phs';
            [~,~,~,nLU] = arnoldiPH(sys,shifts1,opts);
            testCase.verifyEqual(nLU, 5);
            [~,~,~,nLU] = arnoldiPH(sys,shifts2,opts);
            testCase.verifyEqual(nLU, 1);
            [~,~,~,nLU] = arnoldiPH(sys,shifts3,opts);
            testCase.verifyEqual(nLU, 5);

            opts.constAlg = 'sss';
            [~,~,~,nLU] = arnoldiPH(sys,shifts1,opts);
            testCase.verifyEqual(nLU, 5);
            [~,~,~,nLU] = arnoldiPH(sys,shifts2,opts);
            testCase.verifyEqual(nLU, 1);
            [~,~,~,nLU] = arnoldiPH(sys,shifts3,opts);
            testCase.verifyEqual(nLU, 5);
        end

        function phsRedMethod(testCase, sys, shifts)
            redSys = arnoldiPH(sys, shifts);
            testCase.verifyEqual(redSys.method, @arnoldiPH);
        end

        function phsRedParameters(testCase, sys, shifts)
            opts.structurePreservation = 'specialInverse';
            redSys = arnoldiPH(sys,shifts,opts);
            testCase.verifyEqual(redSys.parameters.shifts, shifts);
            testCase.verifyEqual(redSys.parameters.structurePreservation, opts.structurePreservation);

        end

        function phsRedInfo(testCase, sys, shifts)
            redSys = arnoldiPH(sys,shifts);
            testCase.verifyEmpty(redSys.info);
        end
    end

    methods (Test, TestTags = {'Options'})
        % check if options are applied correctly

        function verbose(testCase, sys, shifts)
            % Test non-verbose mode
            % Use command window to verify that no output was created

            warning('on','all')
            disp("Testing non-verbose mode");
            opts.verbose = false;
            opts.phs.verbose = false;
            redSys = arnoldiPH(sys,shifts,opts);
            testCase.verifyEqual(redSys.dim, 5);
            disp("Test non-verbose mode: done");
            warning('off','all')
        end

        function scaling(testCase, sys, shifts)
            % Test structure preserving option 'scaling'

            opts.constAlg = 'phs';
            opts.structurePreservation = 'scaling';
            [redSys,V] = arnoldiPH(sys,shifts,opts);
            testCase.verifyEqual(rank(V),  5);
            testCase.verifyLessThan(norm(eye(size(V'*V)) - V'*V), 1e-14);
            testCase.verifyEqual(redSys.dim, 5);
            % accuracy of transfer function at shifts
            sys_ss = ss(sys);
            redSys_ss = ss(redSys);
            G_original = @(s) sys_ss.C*((s*eye(size(sys_ss.A)) - sys_ss.A)\sys_ss.B);
            G_reduced = @(s) redSys_ss.C*((s*redSys_ss.E - redSys_ss.A)\redSys_ss.B);
            for i = 1:length(shifts)
                delta = abs(G_original(shifts(i)) - G_reduced(shifts(i)));
                testCase.verifyLessThan(delta, 1e-10);
            end

            opts.constAlg = 'sss';
            opts.structurePreservation = 'scaling';
            [redSys,V] = arnoldiPH(sys,shifts,opts);
            testCase.verifyEqual(rank(V),  5);
            testCase.verifyLessThan(norm(eye(size(V'*V)) - V'*V), 1e-14);
            testCase.verifyEqual(redSys.dim, 5);
            % accuracy of transfer function at shifts
            sys_ss = ss(sys);
            redSys_ss = ss(redSys);
            G_original = @(s) sys_ss.C*((s*eye(size(sys_ss.A)) - sys_ss.A)\sys_ss.B);
            G_reduced = @(s) redSys_ss.C*((s*redSys_ss.E - redSys_ss.A)\redSys_ss.B);
            for i = 1:length(shifts)
                delta = abs(G_original(shifts(i)) - G_reduced(shifts(i)));
                testCase.verifyLessThan(delta, 1e-10);
            end
        end

        function specialInverse(testCase, sys, shifts)
            % Test structure preserving option 'specialInverse'

            opts.constAlg = 'phs';
            opts.structurePreservation = 'specialInverse';
            [redSys,V] = arnoldiPH(sys,shifts,opts);
            testCase.verifyEqual(rank(V),  5);
            testCase.verifyLessThan(norm(eye(size(V'*V)) - V'*V), 1e-14);
            testCase.verifyEqual(redSys.dim, 5);
            % accuracy of transfer function at shifts
            sys_ss = ss(sys);
            redSys_ss = ss(redSys);
            G_original = @(s) sys_ss.C*((s*eye(size(sys_ss.A)) - sys_ss.A)\sys_ss.B);
            G_reduced = @(s) redSys_ss.C*((s*eye(size(redSys_ss.A)) - redSys_ss.A)\redSys_ss.B);
            for i = 1:length(shifts)
                delta = abs(G_original(shifts(i)) - G_reduced(shifts(i)));
                testCase.verifyLessThan(delta, 1e-10);
            end

            opts.constAlg = 'sss';
            opts.structurePreservation = 'specialInverse';
            [redSys,V] = arnoldiPH(sys,shifts,opts);
            testCase.verifyEqual(rank(V),  5);
            testCase.verifyLessThan(norm(eye(size(V'*V)) - V'*V), 1e-14);
            testCase.verifyEqual(redSys.dim, 5);
            % accuracy of transfer function at shifts
            sys_ss = ss(sys);
            redSys_ss = ss(redSys);
            G_original = @(s) sys_ss.C*((s*eye(size(sys_ss.A)) - sys_ss.A)\sys_ss.B);
            G_reduced = @(s) redSys_ss.C*((s*eye(size(redSys_ss.A)) - redSys_ss.A)\redSys_ss.B);
            for i = 1:length(shifts)
                delta = abs(G_original(shifts(i)) - G_reduced(shifts(i)));
                testCase.verifyLessThan(delta, 1e-10);
            end
        end

        function Cholesky(testCase, sys, shifts)
            % Test structure preserving option 'Cholesky'

            opts.constAlg = 'phs';
            opts.structurePreservation = 'Cholesky';
            [redSys,V] = arnoldiPH(sys,shifts,opts);
            testCase.verifyEqual(rank(V),  5);
            testCase.verifyEqual(redSys.dim, 5);
            % accuracy of transfer function at shifts
            sys_ss = ss(sys);
            redSys_ss = ss(redSys);
            G_original = @(s) sys_ss.C*((s*eye(size(sys_ss.A)) - sys_ss.A)\sys_ss.B);
            G_reduced = @(s) redSys_ss.C*((s*eye(size(redSys_ss.A)) - redSys_ss.A)\redSys_ss.B);
            for i = 1:length(shifts)
                delta = abs(G_original(shifts(i)) - G_reduced(shifts(i)));
                testCase.verifyLessThan(delta, 1e-10);
            end

            opts.constAlg = 'sss';
            opts.structurePreservation = 'Cholesky';
            [redSys,V] = arnoldiPH(sys,shifts,opts);
            testCase.verifyEqual(rank(V),  5);
            testCase.verifyEqual(redSys.dim, 5);
            % accuracy of transfer function at shifts
            sys_ss = ss(sys);
            redSys_ss = ss(redSys);
            G_original = @(s) sys_ss.C*((s*eye(size(sys_ss.A)) - sys_ss.A)\sys_ss.B);
            G_reduced = @(s) redSys_ss.C*((s*eye(size(redSys_ss.A)) - redSys_ss.A)\redSys_ss.B);
            for i = 1:length(shifts)
                delta = abs(G_original(shifts(i)) - G_reduced(shifts(i)));
                testCase.verifyLessThan(delta, 1e-10);
            end
        end

        function Cholesky_plus(testCase, sys, shifts)
            % Test structure preserving option 'Cholesky+'

            opts.constAlg = 'phs';
            opts.structurePreservation = 'Cholesky+';
            [redSys,V] = arnoldiPH(sys,shifts,opts);
            testCase.verifyEqual(rank(V),  5);
            testCase.verifyEqual(redSys.dim, 5);
            % accuracy of transfer function at shifts
            sys_ss = ss(sys);
            redSys_ss = ss(redSys);
            G_original = @(s) sys_ss.C*((s*eye(size(sys_ss.A)) - sys_ss.A)\sys_ss.B);
            G_reduced = @(s) redSys_ss.C*((s*eye(size(redSys_ss.A)) - redSys_ss.A)\redSys_ss.B);
            for i = 1:length(shifts)
                delta = abs(G_original(shifts(i)) - G_reduced(shifts(i)));
                testCase.verifyLessThan(delta, 1e-10);
            end

            opts.constAlg = 'sss';
            opts.structurePreservation = 'Cholesky+';
            [redSys,V] = arnoldiPH(sys,shifts,opts);
            testCase.verifyEqual(rank(V),  5);
            testCase.verifyEqual(redSys.dim, 5);
            % accuracy of transfer function at shifts
            sys_ss = ss(sys);
            redSys_ss = ss(redSys);
            G_original = @(s) sys_ss.C*((s*eye(size(sys_ss.A)) - sys_ss.A)\sys_ss.B);
            G_reduced = @(s) redSys_ss.C*((s*eye(size(redSys_ss.A)) - redSys_ss.A)\redSys_ss.B);
            for i = 1:length(shifts)
                delta = abs(G_original(shifts(i)) - G_reduced(shifts(i)));
                testCase.verifyLessThan(delta, 1e-10);
            end
        end

        function QV(testCase, sys, shifts)
            % Test structure preserving option 'QV'

            opts.constAlg = 'phs';
            opts.structurePreservation = 'QV';
            [redSys,V] = arnoldiPH(sys,shifts,opts);
            testCase.verifyLessThan(norm(eye(size(V'*V)) - V'*V), 1e-14);
            testCase.verifyEqual(rank(V),  5);
            testCase.verifyEqual(redSys.dim, 5);
            % accuracy of transfer function at shifts
            sys_ss = ss(sys);
            redSys_ss = ss(redSys);
            G_original = @(s) sys_ss.C*((s*eye(size(sys_ss.A)) - sys_ss.A)\sys_ss.B);
            G_reduced = @(s) redSys_ss.C*((s*redSys_ss.E - redSys_ss.A)\redSys_ss.B);
            for i = 1:length(shifts)
                delta = abs(G_original(shifts(i)) - G_reduced(shifts(i)));
                testCase.verifyLessThan(delta, 1e-10);
            end

            opts.constAlg = 'sss';
            opts.structurePreservation = 'QV';
            [redSys,V] = arnoldiPH(sys,shifts,opts);
            testCase.verifyLessThan(norm(eye(size(V'*V)) - V'*V), 1e-14);
            testCase.verifyEqual(rank(V),  5);
            testCase.verifyEqual(redSys.dim, 5);
            % accuracy of transfer function at shifts
            sys_ss = ss(sys);
            redSys_ss = ss(redSys);
            G_original = @(s) sys_ss.C*((s*eye(size(sys_ss.A)) - sys_ss.A)\sys_ss.B);
            G_reduced = @(s) redSys_ss.C*((s*redSys_ss.E - redSys_ss.A)\redSys_ss.B);
            for i = 1:length(shifts)
                delta = abs(G_original(shifts(i)) - G_reduced(shifts(i)));
                testCase.verifyLessThan(delta, 1e-10);
            end
        end

        function orthogonalizationOptions(testCase, sys, shifts)
            % Test orthogonalization options
            % Verify that matrix V has orthogonal columns (if orthogonalized)

            opts.reorthog = false;
            opts.structurePreservation = 'specialInverse'; % needed to avoid postprocessing of V

            opts.constAlg = 'phs';
            %             opts.orth = false;
            %             [redSys,V] = arnoldiPH(sys,shifts,opts);
            %             testCase.verifyEqual(redSys.parameters.constAlg,'phs');
            %             testCase.verifyFalse(redSys.parameters.orth);
            opts.orth = 'mgs';
            [redSys,V] = arnoldiPH(sys,shifts,opts);
            testCase.verifyEqual(redSys.parameters.constAlg,'phs');
            testCase.verifyEqual(redSys.parameters.orth,'mgs');
            testCase.verifyLessThan(norm(V'*V-eye(size(V'*V))), testCase.tolerance)
            opts.orth = '2mgs';
            [redSys,V] = arnoldiPH(sys,shifts,opts);
            testCase.verifyEqual(redSys.parameters.constAlg,'phs');
            testCase.verifyEqual(redSys.parameters.orth,'2mgs');
            testCase.verifyLessThan(norm(V'*V-eye(size(V'*V))), testCase.tolerance)
            opts.orth = 'dgks'; % will throw error
            testCase.verifyError(@()arnoldiPH(sys,shifts,opts),"MORpH:arnoldiPH:bad_option_combination");
            testCase.verifyEqual(redSys.parameters.constAlg,'phs');
            testCase.verifyEqual(redSys.parameters.orth,'2mgs');

            opts.constAlg = 'sss';
            %             opts.orth = false;
            %             [redSys,V] = arnoldiPH(sys,shifts,opts);
            %             testCase.verifyEqual(redSys.parameters.constAlg,'sss');
            %             testCase.verifyFalse(redSys.parameters.orth);
            opts.orth = 'mgs';
            [redSys,V] = arnoldiPH(sys,shifts,opts);
            testCase.verifyEqual(redSys.parameters.constAlg,'sss');
            testCase.verifyEqual(redSys.parameters.orth,'mgs');
            testCase.verifyLessThan(norm(V'*V-eye(size(V'*V))), testCase.tolerance)
            opts.orth = '2mgs';
            [redSys,V] = arnoldiPH(sys,shifts,opts);
            testCase.verifyEqual(redSys.parameters.constAlg,'sss');
            testCase.verifyEqual(redSys.parameters.orth,'2mgs');
            testCase.verifyLessThan(norm(V'*V-eye(size(V'*V))), testCase.tolerance)
            opts.orth = 'dgks';
            [redSys,V] = arnoldiPH(sys,shifts,opts);
            testCase.verifyEqual(redSys.parameters.constAlg,'sss');
            testCase.verifyEqual(redSys.parameters.orth,'dgks');
            testCase.verifyLessThan(norm(V'*V-eye(size(V'*V))), testCase.tolerance)
        end

        %         function LUdecomposition(testCase, sys, sysMIMO)
        %         % Test Opts.lse = {'full','sparse'}
        %
        %             Opts1.constAlg = 'phs';
        %             Opts2.constAlg = 'phs';
        %             Opts1.lse = 'full';
        %             Opts2.lse = 'sparse';
        %             Opts1.reorthog = false;
        %             Opts2.reorthog = false;
        %             Opts1.structurePreservation = 'specialInverse';
        %             Opts2.structurePreservation = 'specialInverse';
        %
        %             % SISO
        %                 % logspace shifts
        %                 shifts = logspace(-5,5,10);
        %                 shifts = [shifts,shifts];   % double shifts such that LU decompositions are perfomed
        %
        %                 [redSys1, V1, W1, nLU1] = arnoldiPH(sys,shifts,Opts1);
        %                 [redSys2, V2, W2, nLU2] = arnoldiPH(sys,shifts,Opts2);
        %
        %                 testCase.verifyLessThan(norm(redSys1-redSys2)/norm(redSys1),testCase.tolerance);
        %                 testCase.verifyLessThan((norm(redSys1-sys)-norm(redSys2-sys))/norm(sys),testCase.tolerance);
        %                 testCase.verifyLessThan(norm(V1-V2)/norm(V1),testCase.tolerance);
        %                 testCase.verifyLessThan(norm(W1-W2)/norm(W1),testCase.tolerance);
        %                 testCase.verifyEqual(nLU1, nLU2);
        %
        %                 % imaginary shifts
        %                 shifts = 1i*logspace(0,5,10);
        %                 shifts = [shifts,shifts];   % double shifts such that LU decompositions are perfomed
        %
        %                 [redSys1, V1, W1, nLU1] = arnoldiPH(sys,shifts,Opts1);
        %                 [redSys2, V2, W2, nLU2] = arnoldiPH(sys,shifts,Opts2);
        %
        %                 testCase.verifyLessThan(norm(redSys1-redSys2)/norm(redSys1),testCase.tolerance);
        %                 testCase.verifyLessThan((norm(redSys1-sys)-norm(redSys2-sys))/norm(sys),testCase.tolerance);
        %                 testCase.verifyLessThan(norm(V1-V2)/norm(V1),testCase.tolerance);
        %                 testCase.verifyLessThan(norm(W1-W2)/norm(W1),testCase.tolerance);
        %                 testCase.verifyEqual(nLU1, nLU2);
        %
        %                 % complex shifts
        %                 shifts = logspace(-5,5,10) + 1i*logspace(0,5,10);
        %                 shifts = [shifts,shifts];   % double shifts such that LU decompositions are perfomed
        %
        %                 [redSys1, V1, W1, nLU1] = arnoldiPH(sys,shifts,Opts1);
        %                 [redSys2, V2, W2, nLU2] = arnoldiPH(sys,shifts,Opts2);
        %
        %                 testCase.verifyLessThan(norm(redSys1-redSys2)/norm(redSys1),testCase.tolerance);
        %                 testCase.verifyLessThan((norm(redSys1-sys)-norm(redSys2-sys))/norm(sys),testCase.tolerance);
        %                 testCase.verifyLessThan(norm(V1-V2)/norm(V1),testCase.tolerance);
        %                 testCase.verifyLessThan(norm(W1-W2)/norm(W1),testCase.tolerance);
        %                 testCase.verifyEqual(nLU1, nLU2);
        %
        %             % MIMO
        %                 % logspace shifts
        %                 shifts = logspace(-5,5,10);
        %                 tangent = rand(size(sysMIMO.G,2),length(shifts));
        %                 shifts = [shifts,shifts];   % double shifts such that LU decompositions are perfomed
        %                 tangent = [tangent,tangent];
        %
        %                 [redSys1, V1, W1, nLU1] = arnoldiPH(sysMIMO,shifts,tangent,Opts1);
        %                 [redSys2, V2, W2, nLU2] = arnoldiPH(sysMIMO,shifts,tangent,Opts2);
        %
        %                 testCase.verifyLessThan(norm(redSys1-redSys2)/norm(redSys1),testCase.tolerance);
        %                 testCase.verifyLessThan((norm(redSys1-sysMIMO)-norm(redSys2-sysMIMO))/norm(sysMIMO),testCase.tolerance);
        %                 testCase.verifyLessThan(norm(V1-V2)/norm(V1),testCase.tolerance);
        %                 testCase.verifyLessThan(norm(W1-W2)/norm(W1),testCase.tolerance);
        %                 testCase.verifyEqual(nLU1, nLU2);
        %
        %                 % imaginary shifts
        %                 shifts = 1i*logspace(0,5,10);
        %                 shifts = [shifts,shifts];   % double shifts such that LU decompositions are perfomed
        %
        %                 [redSys1, V1, W1, nLU1] = arnoldiPH(sysMIMO,shifts,tangent,Opts1);
        %                 [redSys2, V2, W2, nLU2] = arnoldiPH(sysMIMO,shifts,tangent,Opts2);
        %
        %                 testCase.verifyLessThan(norm(redSys1-redSys2)/norm(redSys1),testCase.tolerance);
        %                 testCase.verifyLessThan((norm(redSys1-sysMIMO)-norm(redSys2-sysMIMO))/norm(sysMIMO),testCase.tolerance);
        %                 testCase.verifyLessThan(norm(V1-V2)/norm(V1),testCase.tolerance);
        %                 testCase.verifyLessThan(norm(W1-W2)/norm(W1),testCase.tolerance);
        %                 testCase.verifyEqual(nLU1, nLU2);
        %
        %                 % complex shifts
        %                 shifts = logspace(-5,5,10) + 1i*logspace(0,5,10);
        %                 shifts = [shifts,shifts];   % double shifts such that LU decompositions are perfomed
        %
        %                 [redSys1, V1, W1, nLU1] = arnoldiPH(sysMIMO,shifts,tangent,Opts1);
        %                 [redSys2, V2, W2, nLU2] = arnoldiPH(sysMIMO,shifts,tangent,Opts2);
        %
        %                 testCase.verifyLessThan(norm(redSys1-redSys2)/norm(redSys1),testCase.tolerance);
        %                 testCase.verifyLessThan((norm(redSys1-sysMIMO)-norm(redSys2-sysMIMO))/norm(sysMIMO),testCase.tolerance);
        %                 testCase.verifyLessThan(norm(V1-V2)/norm(V1),testCase.tolerance);
        %                     testCase.verifyLessThan(norm(W1-W2)/norm(W1),testCase.tolerance);
        %                 testCase.verifyEqual(nLU1, nLU2);
        %
        %         end

    end

    methods (Test, TestTags = {'Functionality'})
        % Verify correct outputs for more realistic scenarios

        function matchesMoments(testCase)
            % Verify that model reduction of a more complex system (order >=
            % 100) still guarantees moment matching
            % Verify 1st order moment matching and 2nd order moment matching

            sys = setup_MassSpringDamperSystem(200,10,10,10);
            G = @(s) (sys.G+sys.P)'*sys.Q*((s*sys.E-(sys.J-sys.R)*sys.Q)\(sys.G-sys.P));
            G_der = @(s) (sys.G+sys.P)'*sys.Q*((s*sys.E-(sys.J-sys.R)*sys.Q)\((s*sys.E-(sys.J-sys.R)*sys.Q)\(sys.E*(sys.G-sys.P))));
            %             G_der = @(s) -(sys.G+sys.P)'*sys.Q*(s*sys.E-(sys.J-sys.R)*sys.Q)^(-2)*sys.E*(sys.G-sys.P);

            shifts = repmat(logspace(-5,5,10),1,2);
            opts.constAlg = 'phs';
            opts.lse = 'sparse';

            redSys = arnoldiPH(sys,shifts,opts);

            % 1st order moment matching
            G_red = @(s) (redSys.G+redSys.P)'*redSys.Q*((s*redSys.E-(redSys.J-redSys.R)*redSys.Q)\(redSys.G-redSys.P));

            for i = 1:length(shifts)
                testCase.verifyLessThan(abs(G(shifts(i))-G_red(shifts(i))),testCase.tolerance);
            end

            % 2nd order moment matching
            G_red_der = @(s) (redSys.G+redSys.P)'*redSys.Q*((s*redSys.E-(redSys.J-redSys.R)*redSys.Q)\((s*redSys.E-(redSys.J-redSys.R)*redSys.Q)\(redSys.E*(redSys.G-redSys.P))));
            %             G_red_der = @(s) -(redSys.G+redSys.P)'*redSys.Q*(s*redSys.E-(redSys.J-redSys.R)*redSys.Q)^(-2)*redSys.E*(redSys.G-redSys.P);

            for i = 1:length(shifts)
                testCase.verifyLessThan(abs(G_der(shifts(i))-G_red_der(shifts(i))),testCase.tolerance);
            end
        end

    end

end