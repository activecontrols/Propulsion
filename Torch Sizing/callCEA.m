%% CEA Function
% Author: Yash Amin
%   
% This function will return the entire CEA output to make runtimes faster
% and have a universal CEA run function
%
%% Inputs:
% type: Frozen or Equilibrium -> 'fr' or 'eq'
% Pc: Chamber Pressure -> Can be single value or array of values
% PcUnits: Chamber Pressure units -> 'bar' OR 'psi'
% mixture_type: Tells CEA the amount of fuel and oxidizer i.e. Mixture
% Ratio or Equivalance Ratio -> 'o/f' OR 'phi' Can be single value or array of values
% exit_cond_type: Exit Conditions Type -> 'pip' OR 'sub' OR 'sup'
% exit_cond: Exit Condition -> Pressure Ratio Or Area Ratio Can be single value or array of values
% fuel: Name of Fuel -> must be from propellant list
% fuel_temp_units: Fuel Temperature Units -> 'K', 'C', 'R', or 'F'
% fuel_t: Fuel Temperature
% fuel_wt: Fuel Weight Percentage -> can be in decimal form or percent form
% ox: Name of ox -> must be from propellant list
% ox_temp_units: ox Temperature Units -> 'K', 'C', 'R', or 'F'
% ox_wt: ox Weight Percentage -> can be in decimal form or percent form



function [data_return] = callCEA(type, Pc, PcUnits, mixture_type, mix_val, exit_cond_type, exit_cond, fuel, fuel_temp_units, fuel_t, fuel_wt, ox, ox_temp_units, ox_t, ox_wt)

fclose all;
p = genpath('CEA Calculation Files');
addpath(p)


if ~strcmp(type,'fr') && ~strcmp(type,'eq')
    msg = 'Type not valid, choose between fr or eq';
    error(msg);
end

if ~strcmp(PcUnits,'bar') && ~strcmp(PcUnits,'psi')
    msg = 'Chamber Pressure Unit not valid, choose between bar or psi';
    error(msg);
end

if ~strcmp(mixture_type,'o/f') && ~strcmp(mixture_type,'phi')
    msg = 'Mixture type not valid, choose between o/f or phi';
    error(msg);
end

if ~strcmp(exit_cond_type,'pip') && ~strcmp(exit_cond_type,'sub') && ~strcmp(exit_cond_type,'sup')
    msg = 'Exit Condition type not valid, choose between pip or sub or sup';
    error(msg);
end


if strcmp(type,'fr') 
    pick = 1;
elseif strcmp(type,'eq') 
    pick = 0;
end

run = true;
CEA_RUN = run;
CEA_SAVE_FILE = 'cea.mat';

inp = containers.Map;
inp('type') = type;              % Sets the type of CEA calculation


inp('p') = Pc;                % Chamber pressure
inp('p_unit') = PcUnits;              % Chamber pressure units

if strcmp(mixture_type,'o/f')

    inp('o/f') = mix_val;               % Mixture ratio

elseif strcmp(mixture_type,'phi')

    inp('phi') = mix_val;               % Equivalance Ratio

end

if strcmp(exit_cond_type,'pip')

    inp('pip') = exit_cond;                   % Pressure ratios

elseif strcmp(exit_cond_type,'sub')

    inp('sub') = exit_cond;               % Subsonic area ratios

elseif strcmp(exit_cond_type,'sup')

    inp('sup') = exit_cond;               % Supersonic area ratios

end


inp('fuel') = fuel;             % Fuel name from thermo.inp

inp('fuel_t_units') = fuel_temp_units;

inp('fuel_wt%') = fuel_wt;

inp('fuel_t') = fuel_t;                % Fuel inlet temperature

inp('ox') = ox;

inp('ox_t_units') = ox_temp_units;

inp('ox_wt%') = ox_wt;

inp('ox_t') = ox_t;

% inp('file_name') = 'CEArun_.inp';    % Input/output file name
if CEA_RUN
    data = cea_rocket_run(inp);     % Call the CEA MATLAB code
    save(CEA_SAVE_FILE, 'data');
else
    load(CEA_SAVE_FILE);
end

% The output data structure, called 'data' in this case, is also a MATLAB
% map. 'data' contains a single entry for each of the CEA calculation types
% listed ('eq' and 'fr'). For instance, if only 'fr' is listed, then 'data'
% will only contain a single entry under data('fr').




if strcmp(type,'fr')

    data_return = data('fr');

elseif strcmp(type,'eq')

    data_return = data('eq');

end



fclose all;


end