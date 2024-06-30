% @function:      get_film_cooling_pcnt
% @description:   determine optimal film cooling percent
% @contributors:  Liam Schenk
%                 Cole Cicora
%                 Nicholas Faber
% @inputs:
%          @variable:    @unit:  @info:
%          of_ratio      -       oxidizer to fuel ratio
%          m_dot         lbm/s   total propellant mass flow rate
%          chamber_p     psi     chamber pressure
%          status        BOOL    give script status update
% @outputs:
%          @variable:    @unit:  @info:
%          best_fc_pcnt  -       optimal film cooling percent of mass flow rate
%          opt_fc_diam   in      optimal film cooling hole diameter
%          opt_fc_numh   -       optimal number of film cooling holes
% @example:
%          [fc_pcnt,diam,num_holes] =
%          get_film_cooling_pcnt(1.2,13.2,300,true);
function [best_fc_pcnt, opt_fc_diam, opt_fc_numh] = get_film_cooling_pcnt(of_ratio, m_dot, chamber_p, status)

% Set function paths:
if ismac || isunix
    p_sep = '/';
elseif ispc
    p_sep = '\';
end
function_path = append(pwd, p_sep, 'bin');
cea_path = append(pwd, p_sep, 'PSP_CEA_function_wrapper');
main_path = pwd;
addpath(function_path);
addpath(cea_path)
name_string = 'film_cooling_CEA_file';
input_name = append(name_string, '.inp');
output_name = append(name_string, '.out');

% Define constants/conversion factors:
c = conversions;
bit_sizes = readmatrix('bit_sizes.csv'); % in; diam
gravity = 32.174; % ft/s^2
gravity = gravity * c.F2IN; % in^2/s
gas_const = 8.314; % J/(K * mol)
external_p = 14.696; % psi; sea level
boiling_temp = 352.15; % K; sea level
heat_combustion = 12.5; % BTU/lb
enthalpy_vap = 38.6; % kJ/mol
initial_temp = 76.964; % degF

% Chamber wall properties:
rho_wall = 0.052; % lbm/in^3
k_wall = 0.34591544112; % W/(m * K)
k_wall = k_wall * c.WPMK2BTUSFI; % BTU/(degF * in * s)
cp_wall = 0.38; % BTU/(lbm * degF)
alpha_wall = k_wall / (rho_wall * cp_wall); % in^2/s

% Define additional run parameters:
% @note: the discharge coefficient below is a CONSERVATIVE estimate
%        also, the fuel dP upper limit was provided by fluids
rho_fuel = 0.030911956; % lbm/in^3
dis_coef_fc = 0.3; % N/A
mdot_fuel = m_dot / (of_ratio + 1); % lbm/s
dp_fuel = 0.2 * chamber_p; % psi
cp_fuel = 0.7635; % BTU/(lb * degF)

% Set iteration bounds:
% @changeme: if more coarse/fine solution set desired
% @note: the % film fuel flow rate was chosen from a ~ short ~ lit review
%        5-15% is often used in addition to ablation/regen
%        20-35% if primary cooling method
%        set to 10~11% for CMS
fc_pcnt_range = 0.10:0.002:0.11;
fc_hole_range = 20:5:70;

% Define storage variables:
num_outputs = length(fc_pcnt_range) * length(fc_hole_range);
pcnt_cache = zeros(1, num_outputs);
hole_cache = zeros(num_outputs, 2);
jet_cache = zeros(1, num_outputs);
flux_cache = zeros(1, num_outputs);

% Determine wall temp associated with film cooling %:
iter = 1;
count = 1;
if status
    f = waitbar(0,'Initializing rocket','Name','Current rocket progress','Units','normalized','Position',[0.38 0.3 0.24 0.08]);
end
for fc_pcnt = fc_pcnt_range
    % Determine hole geometry:
    mdot_fc = mdot_fuel * fc_pcnt; % lbm/s
    area_fc = mdot_fc / (sqrt(2 * dp_fuel * rho_fuel * gravity) * dis_coef_fc); % in^2

    for num_fc_holes = fc_hole_range
        % Define real hole size:
        diam_fc = 2 * sqrt(area_fc / (pi * num_fc_holes));       % in
        if diam_fc < (0.95 * bit_sizes(1)) || diam_fc > (1.05 * bit_sizes(end))
            num_outputs = num_outputs - 1;
            pcnt_cache(:,end) = [];
            hole_cache(end,:) = [];
            jet_cache(:,end) = [];
            flux_cache(:,end) = [];
            if status
                waitbar(count / num_outputs,f,sprintf("Setting engine parameters, %.2f%% Complete", 100*count / num_outputs));
            end
            continue
        end
        temp_matrix = repmat(diam_fc, [1 length(bit_sizes)]);
        [~,indexMin] = min(abs(temp_matrix - bit_sizes));
        real_diam_fc = bit_sizes(indexMin);                      % in
        real_area_fc = pi * num_fc_holes * (real_diam_fc / 2)^2; % in^2
        
        % Fetch chamber properties:
        mdot_fc_temp = sqrt(2 * dp_fuel * rho_fuel * gravity) * dis_coef_fc * real_area_fc; % lbm/s
        mdot_fuel_temp = mdot_fuel - mdot_fc_temp;                                          % lbm/s
        of_ratio_temp = (m_dot / mdot_fuel_temp) - 1;                                       % n/a
        inj_vel = mdot_fc_temp / (rho_fuel * real_area_fc);                                 % in/s

        % Solve core flow chamber properties:
        % @note: run at higher of ratio
        [~,~,~,~,chamber_temp,cp_gas,k_gas,enthalpy_gas,rho_gas] = PSP_1DOF_CEA_function_wrapper(chamber_p,external_p,of_ratio_temp,name_string,false);
        delete(append(cea_path, p_sep, input_name));
        delete(append(main_path, p_sep, output_name));
        delete(append(main_path, p_sep, input_name));

        % Convert CEA result units:
        % @variable:     @unit: ~ output from CEA ~
        % chamber_temp:  K
        % cp_gas:        j/(kg*K)
        % k_gas:         W/(m*K)
        % enthalpy_gas:  J/kg
        % rho_gas:       kg/m^3
        chamber_temp = (chamber_temp - 273.15)*1.8 + 32; % degF
        cp_gas = cp_gas * c.JKGK2BTULF;                  % BTU/(lb*degF)
        k_gas = k_gas * c.WPMK2BTUSFI;                   % BTU/(in*s*degF)
        enthalpy_gas = -1 * enthalpy_gas * c.JPK2BTUL;   % BTU/lb
        rho_gas = rho_gas * c.KG2LB / c.CM2CI;           % lb/in^3

        % Perform RT/Weber analysis:
        surf_tension = (dp_fuel * real_diam_fc) / 4; % lb/in
        surf_temp = ((1 / boiling_temp) + (gas_const / (1000*enthalpy_vap))*log(external_p / chamber_p))^(-1); % boiling temp; K
        visc_temp = surf_temp; % K
        surf_temp = (surf_temp - 273.15)*1.8 + 32; % degF
        dynamic_visc = exp(-7.37146 + (2770.25)/(74.6787 + visc_temp)); % mPa*s; vogel eq
        dynamic_visc = dynamic_visc * c.MPS2LIS; % lb/(in * s)
        ohnesorge = dynamic_visc / (rho_fuel * real_diam_fc * surf_tension)^0.5;
        sauter_diam = 1.436 * real_diam_fc * (1+(3*ohnesorge)) ^ (1/6); % in

        % Determine drop lifetime:
        spalding = (cp_gas * (chamber_temp - surf_temp) + (heat_combustion / of_ratio_temp)) / enthalpy_gas; % n/a
        k_const = (8 * k_gas * log(1 + spalding)) / (rho_fuel * cp_gas); % proportionality constant; in^2/s
        drop_lifetime = (sauter_diam^2) / k_const; % s
        k_fuel = 0.0003 * visc_temp^1.8815; % mW/(m * K)
        k_fuel = k_fuel * c.WPMK2BTUSFI * 0.001; % BTU/(degF * in * s)
        prandtl = dynamic_visc * cp_fuel / k_fuel; % n/a

        % Evaluate aerodynamic effects:
        weber = rho_fuel * inj_vel^2 * real_diam_fc / surf_tension;            % n/a
        reynolds = rho_fuel * inj_vel * real_diam_fc / dynamic_visc;           % n/a
        compact_length = 4 * sauter_diam * sqrt(weber) + 0.393701;             % in; accounts for initial section
        if reynolds <= 80
            drag_coef = 27 * reynolds^(-0.84);
        elseif reynolds >= 10000
            drag_coef = 2;
        else
            drag_coef = 0.271 * reynolds^(0.217);
        end
        drag_force = drag_coef * pi * rho_gas * (sauter_diam * inj_vel)^2 / 8; % lbf
        drop_mass = pi * rho_fuel * (sauter_diam^3) / 6;                       % lb
        drop_accel = (drop_mass*gravity - drag_force) / drop_mass;             % in/s^2

        % Determine jet penetration:
        drop_penetration = inj_vel * drop_lifetime + 0.5 * drop_accel * drop_lifetime^2; % in
        jet_penetration = compact_length + drop_penetration;                             % in
        jet_lifetime = drop_lifetime + (compact_length / inj_vel);                       % s

        % Re-run CEA:
        % @note: run at higher of ratio
        [~,c_star,area_ratio,~,chamber_temp,~,k_gas,~,rho_gas] = PSP_1DOF_CEA_function_wrapper(chamber_p,external_p,of_ratio,name_string,false);
        delete(append(cea_path, p_sep, input_name));
        delete(append(main_path, p_sep, output_name));
        delete(append(main_path, p_sep, input_name));

        % Convert CEA result units:
        % @variable:     @unit: ~ output from CEA ~
        % chamber_temp:  K
        % k_gas:         W/(m*K)
        % rho_gas:       kg/m^3
        % c_star:        m/s
        chamber_temp = (chamber_temp - 273.15)*1.8 + 32; % degF
        k_gas = k_gas * c.WPMK2BTUSFI;                   % BTU/(in*s*degF)
        rho_gas = rho_gas * c.KG2LB / c.CM2CI;           % lb/in^3
        c_star = c_star * c.M2IN;                        % in/s

        % Determine heat transfer coefficient:
        % @note: this time with film cooling
        R_gas = chamber_p / (rho_gas * chamber_temp);
        area_throat = m_dot / (chamber_p * k_gas * ((sqrt((2 / (k_gas + 1)) ^ ((k_gas + 1) / (k_gas - 1))) / sqrt(R_gas * k_gas * chamber_temp)))); % in^2
        diam_throat = 2 * sqrt(area_throat / pi); % in
        radius_curv = diam_throat * 1.5 / 2;      % in
        wall_temp_guess = chamber_temp / 2;       % degF
        
        % Iterate over wall temperature:
        while 1
            correction_factor = 1 / (0.5 * (wall_temp_guess / chamber_temp) + 0.5)^0.68;
            heat_trans_coef = ((0.026 * cp_fuel * dynamic_visc^0.2)/(diam_throat * prandtl^0.6)) * (chamber_p * gravity / c_star)^0.8 * (diam_throat / radius_curv)^0.1;
            heat_trans_coef = heat_trans_coef * correction_factor * area_ratio^0.9;
            last_wall_temp = (1 - exp(heat_trans_coef^2 * alpha_wall * jet_lifetime / k_wall^2)*erfc(heat_trans_coef * sqrt(jet_lifetime * alpha_wall) / k_wall)) * (chamber_temp - initial_temp) + initial_temp;
            if abs((last_wall_temp - wall_temp_guess) / last_wall_temp) < 0.001
                break;
            else
                wall_temp_guess = last_wall_temp; % degF
            end
        end
        fc_heat_flux = heat_trans_coef * (chamber_temp - last_wall_temp); % BTU/(s * in^2)

        % Compile results:
        % @changeme: efficiency?
        pcnt_cache(iter) = fc_pcnt;
        hole_cache(iter,1) = num_fc_holes;
        hole_cache(iter,2) = real_diam_fc;
        jet_cache(iter) = jet_penetration;
        flux_cache(iter) = fc_heat_flux;

        % Provide progress, if requested:
        if status
            waitbar(count / num_outputs,f,sprintf("Setting engine parameters, %.2f%% Complete", 100*count / num_outputs));
        end
        iter = iter + 1;
        count = count + 1;
    end
end

% Parse results:
[~,flux_idx] = sort(flux_cache,'ascend');
[~,jet_idx] = sort(jet_cache,'descend');
temp_flux = flux_idx(1);
temp_jet = jet_idx(1);
best_idx = 1;
for idx = 2:num_outputs
    if flux_idx(idx) < temp_flux && jet_idx(idx) < temp_jet
        best_idx = idx;
        temp_flux = flux_idx(idx);
        temp_jet = jet_idx(idx);
    end
end
best_fc_pcnt = pcnt_cache(best_idx);
opt_fc_numh = hole_cache(best_idx,1);
opt_fc_diam = hole_cache(best_idx,2);

% Show results:
cd(main_path)
if status
    close(f)
end

% Clear:
delete(append(pwd, p_sep, 'cea.mat'));
end
