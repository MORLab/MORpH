classdef phsRed < phs
    % PHSRED (class) - offers same functionalities as phs class and adds some
    %                  features for reduced models
    %
    % Syntax:
    %   redSys = PHSRED(J,R,Q,G,E,P,S,N,Opts)
    %   (see phs-class)
    %
    % Description:
    %       This class extends the phs-class. It is useful for reduction
    %       algorithms: It stores all the attributes and offers all functions
    %       provided by the phs class. In addition, it stores the method and
    %       parameters that were used for deriving the reduced model (if
    %       information is provided by the reduction algorithm)
    %
    % Input Arguments:
    %       *Required Input Arguments:*
    %       - J,R,Q,G:  system matrices
    %
    %       *Optional Input Arguments:*
    %       - E,P,S,N:  system matrices
    %       - Opts:     see class 'phs'
    %
    % Additional properties (compared to phs):
    %       - method:       method used for reduction (e.g. function handle)
    %       - parameters:	structure with reduction parameters (depends on method)
    %       - info:         Additional information about the reduced system
    %
    % Output Arguments:
    %       - redSys:   phsRed object
    %
    % See Also:
    %       phs
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

    properties (Access = public)
        parameters
        method
        info
    end

    methods
        function redSys = phsRed(J, R, Q, G, varargin)

            % Use empty constructor and set system properties manually to
            % enable full matrices
            redSys = redSys@phs();
            if nargin == 0, return, end     % empty object

            [J, R, Q, G, E, P, S, N, Opts] = phs.parseInputs(J, R, Q, G, varargin{:});
            redSys.Opts = Opts;
            redSys.J = full(J);
            redSys.R = full(R);
            redSys.Q = full(Q);
            redSys.G = full(G);
            redSys.E = full(E);
            redSys.P = full(P);
            redSys.S = full(S);
            redSys.N = full(N);

            if redSys.Opts.inputValidation
                phs.inputValidation(redSys);
            else
                if Opts.verbose
                    warning('MORpH:phs:noInputValidation',...
                        'No input validation for correctness of port-Hamiltonian system!')
                end
            end

        end %of function phsRed (constructor)
    end %methods public

end