%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%__author__         = "Tim Kayser"
%__date__           = "17.11.2021"
%__version__        = "1.0"
%__maintainer__     = "Tim Kayser"
%__email__          = "kaysert@purdue.edu"
%__status__         = "running, WIP"
%__copyright__      = "Purdue Space Program - Liquids"
%__credits__        = ["NASA CEA", "Purdue CEA MATLAB wrapper ", 
%                       "Tim Kayser"]
%__license__        = "GPL"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Notes for improvements
% - Multiply Isp by .81 for realistic stuff
% - Multiply C* by .9 for realistic stuff
function [isp,cstar,Ae_At_ex,gamma,combustionTemp,cp,conductivity,enthalpy,rho] = PSP_1DOF_CEA_function_wrapper(p_c,p_ex, O_F, filename, debug)

%CEA_function_wrapper => This function takes in a chamber pressure, an
%external pressure and an oxidiser fuel ratio and runs the CEA tool from
%NASA to acquire all the data the CEA tool provides. THis particular
%version returns the specific impulse, as well as the characteristic
%velocity of the engine as outputs. The input of this function is in psi,
%the output is in KMS, taken directly from the CEA tool. This can be easily
%adpated though.
%
%The debug parameter returns a bunch of input and output parameters and
%also checks for negative/surprisingly deviating values
%
%   INPUTS:     UNIT    data type
%       p_c     psi     float OR int
%       p_ex    psi     float OR int
%       O_F     -       float OR int
%       debug   -       boolean
%
%   OUTPUTS:    UNIT    data type
%       Isp     s       float    
%       cstar   m/s     float
%
%   EXAMPLE CALL:
%       [isp, cstar] = CEA_function_wrapper(500, 14.7, 2.6, true)

%% ______ GENRERAL COMMENTS
% According to the "Kalt Bendall" relation for hypersonic
% flow separation the separation pressure
% => p_sep = 0.667 * (P_c / p_a)^(-0.2) * pa 
%      = 5.36 - for a psi of 300 as an upper estimation.


% => THis function assumes the temperature of 90.170 °K for LOX and
%       298.150 °K for RP-1

%% ______ define CEA .inp file parameters
% create input variable - container map
inp = containers.Map;

inp('type') = 'eq';             % Sets the type of CEA calculation (equilibrium (eq) or frozen (fr))
inp('p') = p_c;                 % Chamber pressure
inp('pip') = p_c / p_ex;        % chamber pressure -> P_c / P_ex
inp('p_unit') = 'psi';          % Chamber pressure units
inp('fuel') = {'C2H5OH(L)', 'C8H18(L),isooct'};         % Set Fuel
inp('fuel_wt%') = [98, 2]; % set percentage of fuel/water mixture
inp('ox') = {'O2(L)'};          % Set Oxidizer
inp('o/f') = O_F;               % Set O/F ratio

% set the name of the generated output file according to the parameters
% chosen
% namestring = strcat('PSP_BM1_estimates_pip_', ...
%     num2str(int8(p_c / p_ex)), '_p_c_',...
%     num2str(p_c), '_O_F_', ...
%     num2str(O_F));
namestring = filename;

inp('file_name') = namestring;    % Input/output file name

%% ______ debug output
% creating console output if the debug flag is set:
if debug
    warning on all;
    % check inputs for unviable options:
    % p_c,p_ex, O_F,
    % Check chmaber pressure p_c
    if p_c < 50
        warning('WARNING. \n a chamber pressure of %s.',p_c, 'psi seems to be quite low');
    elseif p_c > 1000
        warning('WARNING. \n I thought we were not conidering such high chamber pressures, I like it');
    end
    
    % check exit pressure p_ex
    if p_ex < 0 || p_c < 0 || O_F < 0
        error('ERROR. \n a negative input for the parameters makes no physical sense')
    end 
    
    % check oxidiser fuel ratio O_F
    if O_F < 1
        warning('WARNING. \na O/F of %s.',O_F, 'seems very low. The optimum is somewhere around 2.6 according to literature');
    elseif O_F > 5
        warning('WARNING. \na O/F of %s.',O_F, 'seems very high. The optimum is somewhere around 2.6 according to literature');
    end 

    % output all parameters
    disp('______INPUT______');
    disp('running CEA with the following parameters:');
    disp('  p_c:');
    disp(p_c);
    disp('chamber pressure unit:');
    disp('   psi');
    disp('pressure ratio:');
    disp(p_c/p_ex);
    disp('fuel & oxidiser:');
    disp('   Jet-A & LOX');
    disp('O/F:');
    disp(O_F);
    disp('type:');
    disp('   eq');
    disp('______OUTPUT______');
    disp('Isp unit:');
    disp('   m/s');
    disp('C* unit:');
    disp('   m/s');
end


%% ______ Run the CEA program
CEA_RUN = true;
CEA_SAVE_FILE = 'cea.mat';  % File name for .mat output file

if CEA_RUN
    data = cea_rocket_run(inp);     % Call the CEA MATLAB code
    save(CEA_SAVE_FILE, 'data');
else
    load(CEA_SAVE_FILE);
end

data_eq = data('eq');

%% ______ extract the desired output data

cstar = squeeze(data_eq('cstar'));      % IMPORTANT
isp = squeeze(data_eq('isp'));
Ae_At = data_eq('ae/at');
gamma = data_eq('gammas');
combustionTemp = data_eq('t');
cp = squeeze(data_eq('cp'));            % for film cooling analysis
conductivity = squeeze(data_eq('k'));   % for film cooling analysis
enthalpy = squeeze(data_eq('h'));       % for film cooling analysis
rho = squeeze(data_eq('rho'));          % for film cooling analysis

% _______ (some) OTHER OUTPUT OPTIONS
% prandtl = squeeze(data_eq('prandtl'));% only needed for heat transfer
% cp = squeeze(data_eq('cp'));          % only needed for heat transfer
% gammas = squeeze(data_eq('gammas'));  % only needed for heat transfer
% ae_at = squeeze(data_eq('ae/at'));    % not wanted at the moment
% ivac = squeeze(data_eq('ivac'));      % not wanted at the moment

%% ______ set outputs
isp = isp(end);
cstar = cstar(end);
Ae_At_ex = Ae_At(end);
gamma = gamma(end);
combustionTemp = combustionTemp(1);
cp = cp(1);
conductivity = conductivity(1);
enthalpy = enthalpy(1);
rho = rho(1);

end

