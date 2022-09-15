classdef test_phsRed < matlab.unittest.TestCase
    % TEST test_phsRed - Tests the functionalities of the phsRed class
    %
    % Press <F5> or enter "runtests("test_phsRed")" to run this testscript
    %
    % Description:
    %   This test verifies the functionalities of the phsRed class.
    %   These scenarios are considered:
    %       - Input: Test if all input combinations are processed correctly
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
        J = {[0 -1; 1 0]};
        R = {[1 0; 0 0]};
        Q = {eye(2)};
        G = {[1; 0]};
        E = {eye(2)};
        P = {[0;0]};
        S = {0};
        N = {0};
    end

    %% TESTS

    methods (Test, TestTags = {'Input'})
        function minimalInput(testCase,J,R,Q,G)
            % Minimal input for creating a reasonable phs-object representing a
            % real reduced PH-system model

            testSys = phsRed(J,R,Q,G);
            testCase.verifyProperties(testSys);
            testCase.verifyEqualSystemMatrices(testSys,J,R,Q,G,testCase.standardE,testCase.standardP,testCase.standardS,testCase.standardN);
            testCase.verifyEqual(testSys.method,[]);
            testCase.verifyEqual(testSys.parameters,[]);
            testCase.verifyFalse(testSys.isMIMO);
            testCase.verifyFalse(testSys.isImplicit);
            testCase.verifyFalse(testSys.isDAE);
        end

        function additionalProperties(testCase,J,R,Q,G)
            % Provide additional properties method, parameters and info

            method = "test";
            params.test = true;
            info.someInfo = 'info';
            testSys = phsRed(J,R,Q,G);
            testSys.method = method;
            testSys.parameters = params;
            testSys.info = info;
            testCase.verifyProperties(testSys);
            testCase.verifyEqual(testSys.method,method);
            testCase.verifyEqual(testSys.info,info);

            method = 'test'; %also test with char array
            testSys = phsRed(J,R,Q,G);
            testSys.method = method;
            testSys.parameters = params;
            testSys.info = info;
            testCase.verifyProperties(testSys);
            testCase.verifyEqual(testSys.method,method);
            testCase.verifyTrue(testSys.parameters.test);
            testCase.verifyEqualSystemMatrices(testSys,J,R,Q,G,testCase.standardE,testCase.standardP,testCase.standardS,testCase.standardN);
            testCase.verifyEqual(testSys.info,info);
        end

        function maximumInput(testCase,J,R,Q,G,E,P,S,N)
            % All intended input parameters used

            method = 'test';
            params.test = true;
            info.someInfo = 'info';
            opts.inputValidation = false;
            testSys = phsRed(J,R,Q,G,E,P,S,N,opts);
            testSys.method = method;
            testSys.parameters = params;
            testSys.info = info;
            testCase.verifyProperties(testSys);
            testCase.verifyEqual(testSys.method,method)
            testCase.verifyEqual(testSys.parameters,params);
            testCase.verifyEqual(testSys.info, info);
            testCase.verifyEqualSystemMatrices(testSys,J,R,Q,G,E,P,S,N);
            phs.inputValidation(testSys);   % should not error
        end

        function onlyMethod(testCase,J,R,Q,G)
            % only method provided, parameters missing

            method = "test";
            testSys = phsRed(J,R,Q,G);
            testSys.method = method;
            testCase.verifyProperties(testSys);
            testCase.verifyEqualSystemMatrices(testSys,J,R,Q,G,testCase.standardE,testCase.standardP,testCase.standardS,testCase.standardN);
            testCase.verifyEqual(testSys.method,method);
            testCase.verifyEqual(testSys.parameters,[]);
            testCase.verifyEqual(testSys.info,[]);
        end

        function onlyParameters(testCase,J,R,Q,G)
            % only parameters provided, method missing

            params.test = true;
            testSys = phsRed(J,R,Q,G);
            testSys.parameters = params;
            testCase.verifyProperties(testSys);
            testCase.verifyEqualSystemMatrices(testSys,J,R,Q,G,testCase.standardE,testCase.standardP,testCase.standardS,testCase.standardN);
            testCase.verifyEqual(testSys.method,[]);
            testCase.verifyEqual(testSys.parameters,params);
        end

        function onlyE(testCase,J,R,Q,G,E)
            % Optional input matrix E provided

            testSys = phsRed(J,R,Q,G,E);
            testCase.verifyProperties(testSys);
            testCase.verifyEqualSystemMatrices(testSys,J,R,Q,G,E,testCase.standardP,testCase.standardS,testCase.standardN);
        end

        function onlyEP(testCase,J,R,Q,G,E,P)
            % Optional input matrices E and P provided

            testSys = phsRed(J,R,Q,G,E,P);
            testCase.verifyProperties(testSys);
            testCase.verifyEqualSystemMatrices(testSys,J,R,Q,G,E,P,testCase.standardS,testCase.standardN);
        end

        function makeSparse(testCase, J, R, Q, G, E, P, S, N)
            % Test phs.makeSparse function
            Opts.inputValidation = false;
            testSys = phsRed(J, R, Q, G, E, P, S, N, Opts);
            testCase.verifyEqualSystemMatrices(testSys, J, R, Q, G, E, P, S, N);
            testSys = testSys.makeSparse;
            testCase.verifyTrue(issparse(testSys.J));
            testCase.verifyEqual(testSys.J, sparse(J));
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
            testCase.verifyEqual(testSys.J,full(J));
            testCase.verifyEqual(testSys.R,full(R));
            testCase.verifyEqual(testSys.Q,full(Q));
            testCase.verifyEqual(testSys.G,full(G));
            testCase.verifyEqual(testSys.P,full(P));
            testCase.verifyEqual(testSys.E,full(E));
            testCase.verifyEqual(testSys.S,full(S));
            testCase.verifyEqual(testSys.N,full(N));
        end
    end
end