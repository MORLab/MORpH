function opts = phsMOR_parseOpts(opts, optsAdmissible, mode)
% PHSMOR_PARSEOPTS - Returns a struct of options whose values were compared
% to a structure with admissible values.
%
% Syntax:
%   opts = PHSMOR_PARSEOPTS(opts,optsAdmissible)
%   opts = PHSMOR_PARSEOPTS(opts,optsAdmissible,mode)
%
% Description:
%   Pass a structure whose values need to be verified with respect to a
%   given set of admissible values specified by another structure
%   optsAdmissible. optsAdmissible must have the same field names for
%   fields which need to be compared.
%
%   For optsAdmissible, use cells to specify a finite set of admissible
%   values or assign a single value to admit all values of the same
%   datatype if you only want to provide a default value.
%
%   If any value does not follow the specified rules, phsMOR_parseOpts will
%   throw an error.
%
%   Note that "string" and 'char' values can be used interchangeably.
%
% Input Arguments:
%       *Required Input Arguments:*
%       - opts:             Structure whose values need to be verified
%       - optsAdmissible:   Structure with same field names and fields of
%                           cells with admissible values or default values.
%                           First cell value will be used as default.
%       *Optional Input Arguments:*
%       - mode: one of the following:
%           - 'fillOnly':   Do not check correctness of provided values but
%                           fill in missing values (which are part of
%                           optsAdmissible but not opts) only.
%
% Output Arguments:
%       - opts:             Structure with verified values (or default values)
%
% Examples:
%       opts.value1 = "Hello";
%       opts.value2 = 2;
%       opts.value3 = 'value3';
%       opts.value5 = 2.5;
%       optsAdmissible.value1 = {"Goodbye","Hello"};
%       optsAdmissible.value2 = {1,3,4};
%       optsAdmissible.value4 = {4,'value4'};
%       optsAdmissible.value5 = 1.0; %can be any other value of same type
%
%       opts_verified = phsMOR_parseOpts(opts,optsAdmissible)
%       % will error because of non-admissible inputs
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

%% Input parsing
if nargin < 3
    mode = 'check';
end

if ~isstruct(opts)
    % throw error if structure is incorrect
    error('MORpH:phsMOR_parseOpts:badPattern',...
        'The passed value for opts does not follow the admissible structure pattern.\nCheck for structure type and correct nesting!');
end

%% Iterate through admissible values
fnames = fieldnames(optsAdmissible);
for i = 1:length(fnames)
    name = fnames{i};

    % Add missing entries
    if ~isfield(opts, name)
        if iscell(optsAdmissible.(name))
            admissibleValues = optsAdmissible.(name);
            opts.(name) = admissibleValues{1};
        elseif isstruct(optsAdmissible.(name))
            opts.(name) = phsMOR_parseOpts(struct(), optsAdmissible.(name));
        else
            opts.(name) = optsAdmissible.(name);
        end
        continue    % Skip checking since default value should be correct
    end

    % Skip checking in fill-only mode
    if isequal(mode, 'fillOnly'), continue; end

    % Check values
    if iscell(optsAdmissible.(name))
        % Finite set of admissible values
        if ~checkValue(opts.(name), optsAdmissible.(name))
            throwAsCaller(inputException(name, opts.(name), "(Value is not in set of admissible values)"))
        end
    elseif isstruct(optsAdmissible.(name))
        % Nested definition! --> recursive call
        opts.(name) = phsMOR_parseOpts(opts.(name), optsAdmissible.(name));
    else
        % Only one default value provided (check data type)
        if ~checkType(opts.(name), optsAdmissible.(name))
            throwAsCaller(inputException(name, opts.(name), "(Wrong data type)"));
            % if numel(opts.(name)) ~= numel(optsAdmissible.(name)) && ~ischar(optsAdmissible.(name)) && ~isstring(optsAdmissible.(name))
            %   throwAsCaller(inputException(name, opts.(name), "(Wrong number of elements)"));
            % end
        end
    end
end


%% =============== AUXILIARY FUNCTIONS =====================
    function valueOK = checkValue(value, admissibleValueCell)
        % Checks if value is in provided cell array
        valueOK = any(cellfun(@(x) isequal(value,x), admissibleValueCell));
    end

    function typeOK = checkType(value, admissibleValue)
        % Checks if value has same data type as admissible value
        % (string and char are treated equal)
        typeOK = isa(value,class(admissibleValue)) ...
            || (ischar(value) && isstring(admissibleValue)) ...
            || (isstring(value) && ischar(admissibleValue));
    end

    function exception = inputException(name, value, hint)
        % Returns exception with message text according to wrong value
        if nargin < 3
            hint = [];
        end
        if ischar(value)
            valueTerm = strcat("('", value, "') ");
        elseif isstring(value)
            valueTerm = strcat("(""", value, """) ");
        elseif isa(value, 'double') && numel(value) == 1
            valueTerm = strcat("(", string(value), ") ");
        else
            valueTerm = [];
        end
        msg = strcat("The given value ", valueTerm, "for ", name, " is not admissible! ", hint);
        exception = MException("MORpH:phsMOR_parseOpts:input_not_admissible",msg);
    end

end
