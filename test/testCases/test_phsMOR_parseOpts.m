classdef test_phsMOR_parseOpts < matlab.unittest.TestCase
    % TEST test_phsMOR_parseOpts - Tests the functionalities of
    %           phsMOR_parseOpts function
    %
    % Press <F5> or enter "runtests("test_phsMOR_parseOpts")" to run this testscript
    %
    % Description:
    %   < Describe what this test is intended to do here. Include, e.g.,
    %     considered test scenarios (grouped by tags)                    >
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

    %% TESTS
    methods (Test, TestTags = {'Input'})
        % Check input combinations

        function allCorrect(testCase)
            % All passed options are compliant with admissible values

            opts = struct();
            opts.value1 = "value1";
            opts.value2 = 2;
            opts.value3 = [1 1; 2 2];

            optsAdmissible = struct();
            optsAdmissible.value1 = {"blabla","value1",1};
            optsAdmissible.value2 = {2,3,[1 1]};
            optsAdmissible.value3 = {"3",3,[1 1; 2 2]};

            opts_verified = phsMOR_parseOpts(opts,optsAdmissible);
            testCase.verifyEqual(opts_verified,opts);
        end

        function optsNotStruct(testCase)
            % the passed opts is not a structure --> error

            opts = "notAStruct";

            optsAdmissible = struct();
            optsAdmissible.value1 = {"blabla","value1",1};
            optsAdmissible.value2 = {2,3,[1 1]};
            optsAdmissible.value3 = {"3",3,[1 1; 2 2]};

            testCase.verifyError(@() phsMOR_parseOpts(opts,optsAdmissible),...
                'MORpH:phsMOR_parseOpts:badPattern');
        end

        function oneMissing(testCase)
            % One value (value3) is missing in the passed opts struct

            opts = struct();
            opts.value1 = "value1";
            opts.value2 = 2;

            optsAdmissible = struct();
            optsAdmissible.value1 = {"blabla","value1",1};
            optsAdmissible.value2 = {2,3,[1 1]};
            optsAdmissible.value3 = {"3",3,[1 1; 2 2]};

            opts_verified = phsMOR_parseOpts(opts,optsAdmissible);

            opts_expected.value1 = "value1";
            opts_expected.value2 = 2;
            opts_expected.value3 = "3";
            testCase.verifyEqual(opts_verified,opts_expected);
        end

        function admissibleNotSpecified(testCase)
            % One additional value (value1) provided which is not verified

            opts = struct();
            opts.value1 = "value1";
            opts.value2 = 2;

            optsAdmissible = struct();
            optsAdmissible.value1 = {"blabla","value1",1};

            opts_verified = phsMOR_parseOpts(opts,optsAdmissible);

            opts_expected.value1 = "value1";
            opts_expected.value2 = 2;
            testCase.verifyEqual(opts_verified,opts_expected);
        end

        function oneWrong(testCase)
            % One wrong value (value2) --> error

            opts = struct();
            opts.value1 = "value1";
            opts.value2 = 5;    % wrong --> can only be one of {2,3,[1 1]}
            opts.value3 = [1 1; 2 2];

            optsAdmissible = struct();
            optsAdmissible.value1 = {"blabla","value1",1};
            optsAdmissible.value2 = {2,3,[1 1]};
            optsAdmissible.value3 = {"3",3,[1 1; 2 2]};

            testCase.verifyError(@() phsMOR_parseOpts(opts,optsAdmissible),...
                "MORpH:phsMOR_parseOpts:input_not_admissible");
        end

        function arbitraryValuesPossible(testCase)
            % optsAdmissible only specifies datatype and default value, passed
            % opts struct is compliant

            opts = struct();
            opts.value1 = 2;
            optsAdmissible = struct();
            optsAdmissible.value1 = 1.0; % only datatype and default value

            opts_verified = phsMOR_parseOpts(opts,optsAdmissible);
            testCase.verifyEqual(opts,opts_verified);
        end

        function arbitraryValuesPossibleWrongDataType(testCase)
            % passed value (opts.value1) has wrong datatype --> error

            opts = struct();
            opts.value1 = "2";
            optsAdmissible = struct();
            optsAdmissible.value1 = 1.0;

            testCase.verifyError(@() phsMOR_parseOpts(opts,optsAdmissible), ...
                "MORpH:phsMOR_parseOpts:input_not_admissible");
        end

        function arbitraryValuesPossibleNotProvided(testCase)
            % optsAdmissible specifies value needed which is not provided by
            % opts --> default value should be chosen

            opts = struct();
            opts.value2 = "2";
            optsAdmissible = struct();
            optsAdmissible.value1 = 1.0;

            opts_verified = phsMOR_parseOpts(opts,optsAdmissible);

            opts_expected.value1 = 1.0;
            opts_expected.value2 = "2";
            testCase.verifyEqual(opts_expected,opts_verified);
        end

        function nestedCorrect(testCase)
            % Nested structs are provided

            opts = struct();
            opts.value1 = "value1";
            opts.value2.value1 = 2.1;
            opts.value2.value2 = "2.2";

            optsAdmissible = struct();
            optsAdmissible.value1 = {"blabla","value1",1};
            optsAdmissible.value2.value1 = {2.1,"2.1"};
            optsAdmissible.value2.value2 = {2.2,"2.2"};

            opts_verified = phsMOR_parseOpts(opts,optsAdmissible);
            testCase.verifyEqual(opts_verified,opts);
        end

        function nestedIncorrect(testCase)
            % Nested opts structure has wrong value

            opts = struct();
            opts.value1 = "value1";
            opts.value2.value1 = 2.1;
            opts.value2.value2 = "wrong";   %wrong

            optsAdmissible = struct();
            optsAdmissible.value1 = {"blabla","value1",1};
            optsAdmissible.value2.value1 = {2.1,"2.1"};
            optsAdmissible.value2.value2 = {2.2,"2.2"};

            testCase.verifyError(@() phsMOR_parseOpts(opts,optsAdmissible), ...
                "MORpH:phsMOR_parseOpts:input_not_admissible");
        end

        function nestedOneMissing(testCase)
            % Nested opts structure does provide value (opts.value2.value2)
            % --> default option chosen

            opts = struct();
            opts.value1 = "value1";
            opts.value2.value1 = 2.1;

            optsAdmissible = struct();
            optsAdmissible.value1 = {"blabla","value1",1};
            optsAdmissible.value2.value1 = {2.1,"2.1"};
            optsAdmissible.value2.value2 = {2.2,"2.2"}; % missing in opts

            opts_verified = phsMOR_parseOpts(opts,optsAdmissible);

            optsExpected = struct();
            optsExpected.value1 = "value1";
            optsExpected.value2.value1 = 2.1;
            optsExpected.value2.value2 = 2.2;
            testCase.verifyEqual(opts_verified,optsExpected);
        end

        function nestedStructureMissing(testCase)
            % Nested opts structure does provide value (opts.value2)
            % --> default option chosen

            opts = struct();
            opts.value1 = "value1";

            optsAdmissible = struct();
            optsAdmissible.value1 = {"blabla","value1",1};
            optsAdmissible.value2.value1 = {2.1,"2.1"}; % missing in opts
            optsAdmissible.value2.value2 = {2.2,"2.2"}; % missing in opts

            opts_verified = phsMOR_parseOpts(opts,optsAdmissible);

            optsExpected = struct();
            optsExpected.value1 = "value1";
            optsExpected.value2.value1 = 2.1;
            optsExpected.value2.value2 = 2.2;
            testCase.verifyEqual(opts_verified,optsExpected);
        end

        function additionalNested(testCase)
            % opts defines one struct value not included in optsAdmissible

            opts = struct();
            opts.value1 = "value1";
            opts.value2.value1 = 2.1;   % missing in opts.admissible
            opts.value2.value2 = "2.2"; % missing in opts.admissible

            optsAdmissible = struct();
            optsAdmissible.value1 = {"blabla","value1",1};

            opts_verified = phsMOR_parseOpts(opts,optsAdmissible);
            testCase.verifyEqual(opts_verified,opts);
        end

        function nestedStructureWrong(testCase)
            % Nested structure pattern is not respected --> error

            opts = struct();
            opts.value1 = "value1";
            opts.value2 = 2;    % should be nested to be compliant

            optsAdmissible = struct();
            optsAdmissible.value1 = {"blabla","value1",1};
            optsAdmissible.value2.value1 = {2.1,"2.1"};
            optsAdmissible.value2.value2 = {2.2,"2.2"}; % missing in opts

            testCase.verifyError(@() phsMOR_parseOpts(opts,optsAdmissible),...
                'MORpH:phsMOR_parseOpts:badPattern');
        end

        function signature(testCase)
            optsAdmissible.value1 = {'hi'};
            opts.value1 = 4;
            try
                phsMOR_parseOpts(opts, optsAdmissible);
            catch exception
                testCase.verifyEqual(exception.identifier, 'MORpH:phsMOR_parseOpts:input_not_admissible');
            end
        end

        function fillOnly(testCase)
            % Test fill-only mode
            opts = struct();
            opts.value1 = "value1";
            opts.value2 = "wrong";
            opts.value5 = 7;

            optsAdmissible = struct();
            optsAdmissible.value1 = {"blabla","value1",1};
            optsAdmissible.value2 = {2,3,[1 1]};
            optsAdmissible.value3 = {"3",3,[1 1; 2 2]};
            optsAdmissible.value4 = 4;
            optsAdmissible.value5 = 5;

            opts_verified = phsMOR_parseOpts(opts,optsAdmissible,'fillOnly');

            opts_expected.value1 = "value1";
            opts_expected.value2 = "wrong";
            opts_expected.value3 = "3";
            opts_expected.value4 = 4;
            opts_expected.value5 = 7;
            testCase.verifyEqual(opts_verified,opts_expected);
        end

    end
end