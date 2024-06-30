%
% Chug Stability Code v1.7
% This program solves the characteristic equation in NASA TN D-3080 (1965)
% and NASA TN D-7026 (1970). The code extends the 1965 formulation by
% considering the complex impedance of the feedsystems (as in the 1970
% formulation), but it also considers the fuel and oxidizer flow
% correlations used in the 1965 formulation. The final solution is a plot
% with the pressure drop ratios on the axes and a curve that defines the
% neutral stability curve.
% The characteristic equation from NASA TN D-3080 uses the linearized
% orifice impedance only as the injector impedance. The
% characteristic equation from NASA TN D-7026 and SP-194 are more general
% cases with the injector impedance defined as a complex feedsytem
% impedance. In this model, the exponential coefficients from the 1965
% formulation used in correlating the data are kept in the final solution,
% but the solution still generally considers the complex feedsystem.
%
% For certain cases the general formulation is trivial (e.g. for equal
% time lags, infinite manifold compliance, and zero inertance). This
% limiting solution should be examined with the feed system off, or for
% near limiting cases using a very small frequency resolution. Newton
% Rhapson is used for the frequency solution in the single time lag
% formulations, otherwise the results are parametrically obtained. All
% analytical solutions were derived using Maple.
% See Matthew Casiano PhD Dissertation for details.
%
% The General Form of the Characteristic Equation
% X = 1+(1+MR)*dcdMR/cstar
% Y = 1-MR*(1+MR)*dcdMR/cstar;
% H(s)=1-1/(1+theta*s)*(pc/m)*(X/Zo*exp(-s*tau_o)+Y/Zf*exp(-s*tau_f));
%
% Script v1.7 Update 3/25/10, M. Casiano
%
% To do list:
% 1) Build a GUI for the single time lag methodologies as a separate tab
% 2) Allow user input save file name, also include units in text files
% 3) Other long term things- variable time lag, pc vs MR plot, incorporate
% CEA to minimize inputs and allow throttling automatically.
function CSCv17
clear % Clear all variables from workspace
global fileinp
fileinp=0; % Initialize flag indicating inputs not taken from input file.
% Redisplay above help text so it can be read in compiled version
disp(' ')
disp('Chug Stability Code v1.7')
disp('This program solves the characteristic equation in NASA TN D-3080 (1965)')
disp('and NASA TN D-7026 (1970). The code extends the 1965 formulation by')
disp('considering the complex impedance of the feedsystems (as in the 1970')
disp('formulation), but it also considers the fuel and oxidizer flow')
disp('correlations used in the 1965 formulation. The final solution is a plot')
disp('with the pressure drop ratios on the axes and a curve that defines the')
disp('neutral stability curve.')
disp(' ')
disp('The characteristic equation from NASA TN D-3080 uses the linearized')
disp('orifice impedance only as the injector impedance. The')
disp('characteristic equation from NASA TN D-7026 and SP-194 are a more general')
disp('case with the injector impedance defined as a complex feedsytem')
disp('impedance. In this model, the exponential coefficients from the 1965')
disp('formulation used in correlating the data are kept in the final solution,')
disp('but the solution still generally considers the complex feedsystem.')
disp(' ')
disp('For certain cases the general formulation is trivial (e.g. for equal')
disp('time lags, infinite manifold compliance, and zero inertance). This')
disp('limiting solution should be examined with the feed system off, or for')
disp('near limiting cases using a very small frequency resolution. Newton')
disp('Rhapson is used for the frequency solution in the single time lag')
disp('formulations, otherwise the results are parametrically obtained. All')
disp('analytical solutions were derived using Maple.')
disp('See Matthew Casiano PhD Dissertation for details.')
disp(' ')
disp('The General Form of the Characteristic Equation')
disp('X = 1+(1+MR)*dcdMR/cstar')
disp('Y = 1-MR*(1+MR)*dcdMR/cstar')
disp('H(s)=1-1/(1+theta*s)*(pc/m)*(X/Zo*exp(-s*tau_o)+Y/Zf*exp(-s*tau_f))')
disp(' ')
disp('Script v1.7 Update 3/25/10, M. Casiano')
disp(' ')
disp(' ')
disp('To do list:')
disp('1) Build a GUI for the single time lag methodologies as a separate tab')
disp('2) Allow user input save file name, also include units in text files')
disp('3) Other long term things- variable time lag, pc vs MR plot, incorporate')
disp(' CEA to minimize inputs and allow throttling automatically.')
disp(' ')
ChugGUI % Call GUI