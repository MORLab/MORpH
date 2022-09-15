function sys_frd = frd(sys, varargin)
% FRD - Converts a PH system model (phs object) to a Frequency Response
%       Data model (frd).
%
% Syntax:
%   sys_frd = sys.FRD
%   sys_frd = FRD(sys)
%   sys_frd = FRD(sys, freqs)
%   sys_frd = FRD(sys, resp, freqs)
%
% Description:
%       FRD computes the Frequency Response Data model of the phs-object
%       sys. You can provide the frequencies where the response should be
%       evaluated or you can even provide the response and the frequencies.
%       If not provided, frd will call freqresp(sys) to compute the
%       frequencies/response data.
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:      phs-object
%       *Optional Input Arguments:*
%       - freqs:    Vector of real-valued frequencies
%       - resp:     Vector of (complex) response belonging to freqs. Note
%                   that you must pass resp as second argument and freqs as
%                   third if you want to provide response data.
%
% Output Arguments:
%       - sys_frd:  Frequency Response Data model
%
% Examples:
%       The following code creates a simple mass-spring-damper model,
%       computes the frequency response at ten logarithmically distributed
%       frequency points and converts this data to an frd object. It then
%       calls the bode function.
%
%       sys = setup_MassSpringDamperSystem(20,1,1,0.1)
%       freqs = logspace(-5,5,10);
%       resp = sys.freqresp(freqs);
%
%       sys_frd = frd(sys, resp, freqs)
%       bode(sys_frd);
%
% See Also:
%       frd, phs, DynamicSystem/freqresp, phs/freqresp
%
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

switch length(varargin)
    case 0
        [response, frequencies] = freqresp(sys);
    case 1  % frequencies provided
        frequencies = varargin{1};
        response = freqresp(sys, frequencies);
    otherwise   % at least frequencies and response provided
        response = varargin{1};
        frequencies = varargin{2};
end
sys_frd = frd(response, frequencies);

end