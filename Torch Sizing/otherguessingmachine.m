%%Author: Rahaf Adi
function PC=otherguessingmachine(At_d,exit_area_d, mdot, OF_array)
%%exit area=0.1925in dimatere=0.116
%%At=0.1925 under exi
%%mdot_int0.01077734592
%%0.0033 lbs per s (make sure this is in metric)
%% Throat area guess
%At_d=0.5
%exit_area_d=0.116
%mdot=0.0107777
%mdot conversion
mdot=mdot*0.45359237
counter = 0;
At_input=((At_d/2)^2)*pi();
exit_area=((exit_area_d/2)^2)*pi();
%OF_array = 1:0.5:3;
area_ratio=exit_area/At_input;
m=1
PC=zeros(1,m)
for OF=1:length(OF_array)
    convergence = 0;
    Pc_Max = 5000; % [psi]
    Pc_Min = 0; % [psi]
    pc_At_guess= (Pc_Max+Pc_Min)/2 ;
    counter=0;
    while ~(convergence)
        % CEA call,
        [CEAdata] = callCEA('fr', pc_At_guess, 'psi', 'o/f', OF_array(OF), 'sup', 0, 'H2', 'K', 297, 100, 'O2', 'K', 297, 100);
        pc_At_guess_SI = pc_At_guess* 6894.76
        cstar = squeeze(CEAdata('cstar'))
        At_guess = (cstar * mdot / pc_At_guess_SI) * 1550; % [m^2] to [in^2]
       
        %%if pc increases at decreases
        if abs(At_guess - At_input) > (0.00001) && counter < 100
            if At_guess > At_input
                Pc_Min= pc_At_guess;
            else
                Pc_Max = pc_At_guess;
            end
            pc_At_guess = (Pc_Max + Pc_Min) / 2 % Bisection step
            counter = counter + 1
        else
            convergence = 1;
            PC(1,m)= pc_At_guess;
            mach(1,m)=CEAdata('mach')
            m=m+1;
        end
        abs(At_guess - At_input)
    end
    %% chamber Outputs (plots)
    %M = squeeze(CEAdata('mach'));
    %M=M(2);
    %T = squeeze(CEAdata('t'));
    %plot(OF,T)
    %gamma = squeeze(CEAdata('gammas'));
    %fprintf('Chamber Pressure (Pc): %f psi\n', pc_At_guess);
    %fprintf('Mach Number (M): %f\n', M);
    %fprintf('Temperature (T): %f K\n', T);
    %fprintf('Gamma: %f\n', gamma);
end
 
 