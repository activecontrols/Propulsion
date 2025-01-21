%%Author: Rahaf Adi
function output =otherguessingmachine(At_d,exit_area_d, mdot, OF_array)
%% Throat area guess
mdot=mdot*0.45359237;
CEA_input_name = help;
output = [];
P_e = 14.7;
counter = 0;
At_input=((At_d/2)^2)*pi();
exit_area=((exit_area_d/2)^2)*pi();
D_t = At_d * 0.0254;
R_of_curve = 1.5 * D_t / 2;
%OF_array = 1:0.5:3;
area_ratio=exit_area/At_input;
m=1;
heatflux_factor = 1;
PC=zeros(1,m);
h = zeros(1,m);
Tg = zeros(1,m);
Tg_tot = zeros(1,m);
h_tot = zeros(1,m);
for OF=1:length(OF_array)
    convergence = 0;
    Pc_Max = 5000; % [psi]
    Pc_Min = 0; % [psi]
    pc_At_guess= (Pc_Max+Pc_Min)/2 ;
    counter=0;
    while ~(convergence)
        % CEA call,
        [CEAdata] = callCEA('fr', pc_At_guess, 'psi', 'o/f', OF_array(OF), 'sup', area_ratio, 'H2', 'K', 297, 100, 'O2', 'K', 297, 100);
        pc_At_guess_SI = pc_At_guess* 6894.76;
        cstar = squeeze(CEAdata('cstar'));
        At_guess = (cstar(1) * mdot / pc_At_guess_SI) * 1550; % [m^2] to [in^2]
       
        %%if pc increases at decreases
        if abs(At_guess - At_input) > (0.00001) && counter < 100
            if At_guess > At_input
                Pc_Min= pc_At_guess;
            else
                Pc_Max = pc_At_guess;
            end
            pc_At_guess = (Pc_Max + Pc_Min) / 2; % Bisection step
            counter = counter + 1;
        else

            convergence = 1;
            PC(1,m)= pc_At_guess;
            h = squeeze(CEAdata('h'));
            h_tot(1,m) = -h(1);
            T = squeeze(CEAdata('t')); 
            Tg_tot(1,m) = T(1);
            Tg = T(2);
            M = squeeze(CEAdata('mach'));
            cp = squeeze(CEAdata('cp')); 
            cp_g_tot = cp(1);
            mu = squeeze(CEAdata('visc'));
            mu_g_tot = mu(1); 
            Pr = squeeze(CEAdata('prandtl')); 
            Pr_g_tot = Pr(1);
            gamma = squeeze(CEAdata('gammas'));
            sigma = (.5 * Tg / Tg_tot(1,m) * (1 + (gamma(2) - 1) / 2 * M(2) ^ 2) + .5) ^ -.68 * (1 + (gamma(2) - 1) / 2 * M(2) ^ 2) ^ -.12; % film coefficient correction factor [N/A] (Huzel & Huang 86).
            h_g(1,m) = heatflux_factor * (0.026 / D_t ^ 0.2) * (mu_g_tot ^ 0.2 * cp_g_tot / Pr_g_tot ^ 0.6) * (PC(1,m) / cstar(1)) ^ 0.8 * (D_t / R_of_curve) ^ 0.1 * (1 / area_ratio) ^ .9 * sigma; % gas film coefficient [W/m^2-K] - bartz equation (Huzel & Huang 86).
            m=m+1;
        end
        abs(At_guess - At_input)
    end
    %% chamber Outputs (plots)
   
    output = [PC; cp_g_tot; Tg; h_g];

end
 
 