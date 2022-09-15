%% Installation of MORpH toolbox
%
% This file runs a setup utility to install the MORpH toolbox and (if desired) 
% relevant third-party toolboxes. You can decide which toolboxes to install 
% during execution via user dialogs. The script will detect already installed toolboxes
% automatically.
%
% Toolboxes will be installed here: /src/thirdParty/external
% 
% Potentially installed toolboxes are:
%   - M-M.E.S.S., see <a href="https://github.com/mpimd-csc/mmess">https://github.com/mpimd-csc/mmess</a>
%   - CVX, see <a href="https://github.com/cvxr/CVX">https://github.com/cvxr/CVX</a>
%   - Manopt, see <a href="https://github.com/NicolasBoumal/manopt">https://github.com/NicolasBoumal/manopt</a>
%   - GRANSO, see <a href="https://gitlab.com/timmitchell/GRANSO/">https://gitlab.com/timmitchell/GRANSO/</a>
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

%% Print welcome text
home
header("MORPH INSTALLATION", 'title')

fprintf(['Welcome to the setup utility of the MORpH toolbox! \n' ...
    'This script will help you to set up the MATLAB path and install optional third-party software.\n\n'])


% Find path of MORpH toolbox from location of this file
[morph_path, filename, extension] = fileparts(mfilename('fullpath'));

%% Add MORpH toolbox to MATLAB path
if ~exist("phs", 'class')
    addpath(genpath(morph_path))
    rmpath(genpath(strcat(morph_path,"\.git")))
    savepath
    printStatus('INFO', 'Added MORpH toolbox to MATLAB path.\n');
else
    printStatus('OK', 'MORpH toolbox is already on path.\n');
end

%% Control System Toolbox
installation_control = checkToolboxInstallation("Control System Toolbox");
if installation_control
    printStatus('OK', "Installation of Control System Toolbox found\n");
else
    header("CONTROL SYSTEM TOOLBOX")
    printStatus('INFO', "Could not find an installation of Control System Toolbox");
    fprintf(strcat("It seems that you have not installed the MATLAB Control System Toolbox. (If you have, you can ignore this message)\n",...
                   "The MORpH toolbox uses some of the functionalities of the Control System Toolbox, including e.g. \n",...
                   "    - bode plots\n",...
                   "    - norm commands\n",...
                   "    - lsim\n",...
                   "We try to provide as many alternatives as possible such that most of the features of the MORpH toolbox can also be used\n",...
                   "without the Control System Toolbox. Nevertheless, we strongly recommend installing it (e.g. by using the 'Add-Ons Manager'\n",...
                   "which you can find in the MATLAB 'Home' section) to get the 'full experience'"));
end

%% Third-party Software Download and Installation
cd(strcat(morph_path,"\src\thirdParty\external\"))

% Third-party settings
thirdParty.mmess = struct('fullname','M.E.S.S.','download_cmd','git clone https://github.com/mpimd-csc/mmess.git','install_cmd','mmess\mess_path.m',...
    'descr','solving large matrix equations','usedBy','balPH, lyapOpt, prbt, sfmor and ss2phs');
thirdParty.cvx = struct('fullname','CVX','download_cmd','git clone https://github.com/cvxr/CVX.git','install_cmd','CVX\cvx_setup.m',...
    'descr','solving convex optimization problems','usedBy','prlPasEnf and ss2phs');
thirdParty.manopt = struct('fullname','Manopt','download_cmd','git clone https://github.com/NicolasBoumal/manopt.git','install_cmd','manopt\importmanopt.m',...
    'descr','solving optimization problems on manifolds','usedBy','lyapOpt and prOpt');
thirdParty.granso = struct('fullname','GRANSO','download_cmd','git clone https://gitlab.com/timmitchell/GRANSO.git','install_cmd',[],...
    'descr','solving constrained nonsmooth optimization problems','usedBy','sobmor and ihaPH');

fn = fieldnames(thirdParty);
for i=1:numel(fn)
    try
        thirdPartyCheck(thirdParty.(fn{i}).fullname);
        printStatus('OK', strcat(thirdParty.(fn{i}).fullname," is already installed.",{newline}));
    catch
        intro = strcat("Do you want to download and install the ", thirdParty.(fn{i}).fullname," toolbox?",{newline},...
                       "It helps with ",thirdParty.(fn{i}).descr," and is used by the functions ", thirdParty.(fn{i}).usedBy,".",{newline},...
                       "Please note that by downloading this software, you agree to comply with its specific license!",{newline});
        download_MESS = askConsent(intro);
        if download_MESS
            [status,cmdout] = system(thirdParty.(fn{i}).download_cmd);
            if status
                printStatus('ERROR',strcat("Download of ",thirdParty.(fn{i}).fullname," toolbox failed. ",...
                    'The download requires a <a href="https://git-scm.com/">GIT</a> installation and internet connection.',{newline}));
            else
                if ~isempty(thirdParty.(fn{i}).install_cmd)
                    try
                        run(thirdParty.(fn{i}).install_cmd);
                        printStatus('OK', strcat("Finished installation of ",thirdParty.(fn{i}).fullname," toolbox.",{newline}));
                    catch
                        printStatus('ERROR',strcat("Installation of ",thirdParty.(fn{i}).fullname," toolbox failed.",{newline},...
                            "Please install the toolbox manually.",{newline}));
                    end
                else
                    newpath = strcat(morph_path,'\src\thirdParty\external\',thirdParty.(fn{i}).fullname);
                    addpath(genpath(newpath))
                    rmpath(genpath(strcat(newpath,"\.git")))
                    printStatus('OK', strcat("Finished installation of ",thirdParty.(fn{i}).fullname," toolbox.",{newline}));
                end
            end
        else
            printStatus('SKIPPED', strcat("Did not install ",thirdParty.(fn{i}).fullname," toolbox.",{newline}));
        end
    end
end
savepath;
cd(morph_path)

%% Finish
fprintf("\nSetup finished - enjoy working with the MORpH toolbox! \n")

%% ================= AUXILIARY FUNCTIONS ================================
function consent = askConsent(question)
    reply = input(strcat(question, "\n[Y/N] >>> "), 's');
    if strcmp(reply, 'Y') || strcmp(reply, 'y')
        consent = true;
    else
        consent = false;
    end
end

function printStatus(status, message)
    fprintf(strcat("[", status, "]    ", message, "\n"));
end

function installationFound = checkToolboxInstallation(name)
    addons = matlab.addons.installedAddons;
    installationFound =  max(ismember(addons.Name, name) .* addons.Enabled);
end

function header(text, opt)
    if nargin <= 1
        opt = 'normal';
    end
    switch opt
        case 'normal'
            fprintf(strcat(text, "\n"));
            fprintf(strcat(repmat('-', 1, strlength(text)), "\n"));
        case 'title'
            fprintf(strcat(repmat('=', 1, strlength(text)+4), "\n"));
            fprintf(strcat("| ", text, " |\n"));
            fprintf(strcat(repmat('=', 1, strlength(text)+4), "\n"));
        otherwise
            error('unknown header option')
    end
end
