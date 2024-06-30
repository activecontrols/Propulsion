%{
Purdue Space Program - Liquids
Rocket 3 1DoF - Conversions
Talal Zaim
Liam Schenk

Conversion factors between various units
%}

function c = conversions()

    %% Lengths
    %NM  :Nautical Mile (nm)
    %MI  :Miles (mi)
    %KM  :kilometers (km)
    %M   :meters (m)
    %LTS :light seconds (ls)
    %FT  :feet (ft)
    %IN  :inches (in)
    
    c.NM2KM = 1.8520009; % 1 nautical mile to kilometers
    c.KM2NM = 1 / c.NM2KM; % 1 kilometer to nautical miles
    c.MI2KM = 1.609344016; % 1 mile to kilometers
    c.KM2MI = 1 / c.MI2KM; % 1 kilometer to miles
    c.M2FT = 3.28084; %1 meter to feet
    c.FT2M = 1 / c.M2FT; % 1 foot to meters
    c.FT2IN = 12; % 1 foot to inches
    c.IN2FT = 1 / c.FT2IN; % 1 inch to feet
    c.M2KM = .001; %1 meter to kilometers
    c.KM2M = 1 / c.M2KM; %1 kilometer to meters
    c.KM2FT = 3280.84; % 1 kilometer to feet
    c.FT2KM = 1 / c.KM2FT; % 1 foot to kilometers
    c.IN2M = c.IN2FT * c.FT2M; % 1 inch to meters
    c.M2IN = 1 / c.IN2M; % 1 meter to inches
    c.KM2LTS = 0.00000333564095; % 1 kilometer to lightseconds
    c.LTS2KM = 1 / c.KM2LTS; % 1 lightsecond to kilometers
    c.KM2LTY = .0000000000001057001; % 1 kilometer to lightyears
    c.LTY2KM = 1 / c.KM2LTY; % 1 lightyear to kilometers
    
            % Film Cooling Integration
    %F  :Feet (ft)
    c.F2M  = 0.3048;            % Meters per foot
    c.M2F  = 1 / c.F2M;         % Feet per meter
    c.F2IN = 12.0;              % Inches per foot
    c.IN2F = 1 / c.F2IN;        % Feet per inches
    c.MI2F = 5280;              % Feet per mile
    c.F2MI = 1 / c.MI2F;        % Miles per foot
    c.M2MI  = c.M2F * c.F2MI;   % Miles per meter
    c.MI2M  = 1 / c.M2MI;       % Meters per mile
    c.F2KM = c.F2M * c.M2KM;    % Feet per kilometer
    c.KM2F = 1 / c.F2KM;        % Kilometers per foot
    
    %% Temperatures
    %K  :Kelvin (K)
    %C  :Celcius (C)
    %F  :Fahrenheit (F)
    %R  :Rankine (R)
    c.K2R = 1.8; % 1 Kelvin to Rankine
    c.R2K = 1/ c.K2R; %1 Rankine to Kelvin
    % c.C2F = (Fahrenheit - 32) * 5 / 9; % 1 Celcius to Fahrenheit
    % c.F2C = (Celcius * 9 / 5) + 32; % 1 Fahrenheit to Celcius
    
    %% Pressures
    %PSI : pounds per square inch (psi)
    %PA  : Pascal (Pa)
    %ATM : Atmospheres (Atm)
    %BARS: Bars (Bar)
    c.PSI2PA = 6894.76; % 1 psi to Pascals
    c.PA2PSI = 1 / c.PSI2PA; % 1 Pascal to psi
    c.ATM2BARS = 1.01325; % 1 atmosphere to bars
    c.BARS2ATM = 1 / c.ATM2BARS; % 1 bar to atmospheres
    c.ATM2PA = 101325; % 1 atmosphere to pascals
    c.PA2ATM = 1 / c.ATM2PA; % 1 pascal to atmosphers
    c.ATM2PSI = 14.69595; % 1 atmosphere to psi
    c.PSI2ATM = 1 / c.ATM2PSI; % 1 psi to atmospheres
    c.BARS2PA = 100000; % 1 bar to pascals
    c.PA2BARS = 1 / c.BARS2PA; % 1 pascal to bars
    c.BARS2PSI = 14.50377; % 1 bar to pascals
    c.PSI2BARS = 1 / c.BARS2PSI; % 1 psi to bars
    
            % Film Cooling Integration
    %B  :Bars (Bar)
    c.ATM2B = 1.01325;               % Bars per atmosphere
    c.B2ATM = 1 / c.ATM2B;           % Atmospheres per Bar
    c.B2PA = 100000;                 % Pascal per Bar
    c.PA2B = 1 / c.B2PA;             % Bars per Pascal
    c.B2PSI = 14.50377;              % PSI per Bar
    c.PSI2B = 1 / c.B2PSI;           % Bar per PSI
    
    
    %% Mass
    %KG   :Kilograms (kg)
    %SLUGS:Slugs (slug)
    %LBM  :Pound-mass (lbs)
    c.KG2SLUGS = 0.0685218; % 1 kilogram to slugs
    c.SLUGS2KG = 1 / c.KG2SLUGS; %1 slug to kilograms
    c.LBM2KG = 0.45359; % 1 pound to kilograms
    c.KG2LBM = 1 / c.LBM2KG; % 1 kilogram to pounds
    c.SLUGS2LBM = 32.17405; % 1 slug to pounds
    c.LBM2SLUGS = 1 / c.SLUGS2LBM; %1 pound to slugs
    c.KG2G = 1000; % 1 kilogram to grams
    c.G2KG = 1 / c.KG2G; % 1 gram to kilograms
    
    % Film Cooling Integration
    %LB     :Pound-mass (lbs)
    %S      :Slugs (slug)
    %G      :Grams (gram)
    c.LB2KG = 0.45359237;       % Kilograms per pound
    c.KG2LB = 1 / c.LB2KG;      % Pounds per kilogram
    c.KG2S = 0.0685218;         % Slugs per kilogram
    c.S2KG = 1 / c.KG2S;        % Kilograms per slug
    c.S2LB = 32.17405;          % Pounds per slug
    c.LB2S = 1 / c.S2LB;        % Slugs per pound
    c.KG2G = 1000;              % Grams per kilogram
    c.G2KG = 1 / c.KG2G;        % Kilograms per gram
    
    %% Forces
    %N   :Newtons (N)
    %LBF :Pound Force (lbf)
    c.N2LBF = 0.2248089; % 1 Newton to pound forces
    c.LBF2N = 1 / c.N2LBF; % 1 pound force to newtons
    
    %% Angles
    %D  :Degrees 
    %RAD:Radians (rad)
    c.D2RAD = 0.01745329; % 1 degree to radians
    c.RAD2D = 1 / c.D2RAD; % 1 radian to degrees
    
    % Film Cooling Integration
    % DEG:  Degrees
    c.RAD2DEG = 180 / pi;       % Degrees per radian
    c.DEG2RAD = pi / 180;       % Radians per degree
    
    %% Area
    %SQM :square meter (m^2)
    %SQFT:square foot (ft^2)
    c.SQM2SQFT = 10.76391; %1 square meter to square feet
    c.SQFT2SQM = 1 / c.SQM2SQFT; %1 square foot to square meters
    c.SM2SF = 10.76391;         % Square feet per square meter
    
    % Film Cooling Integration
    % SM: square meter [m^2]
    % SF: square foot [ft^2]
    % SI: square inch [in^2]
    c.SF2SM = 1 / c.SM2SF;      % Square meter per square foot
    c.SF2SI = 144;              % Square inches per square foot
    c.SI2SF = 1 / c.SF2SI;      % Square feet per square inch
    
    %% Time
    %DAY :Day
    %HR  :Hour (hr)
    %MIN :minutes (min)
    %SEC :seconds (sec)
    c.DAY2HR = 24; % 1 day to hours
    c.HR2DAY = 1 / c.DAY2HR; % 1 hour to days
    c.HR2MIN = 60; % 1 hour to minutes
    c.MIN2HR = 1 / c.HR2MIN; % 1 minute to hours
    c.MIN2SEC = 60; % 1 minute to seconds
    c.SEC2MIN = 1 / c.MIN2SEC; % 1 minute to seconds
    
    %% Energy
    %J  :joule (J)
    %KJ :kilojoule (KJ)
    %WH :Watt Hour (watt hr)
    c.J2KJ = 0.001; % 1 joule to kilojoules
    c.KJ2J = 1 / c.J2KJ; % 1 kilojoule to joules
    c.KJ2WH = 0.2777778; % 1 kilojoule to watt hours
    c.WH2KJ = 1 / c.KJ2WH; % 1 watt hour to kilojoules
    c.J2WH = 0.0002777778; % 1 joule to watt hours
    c.WH2J = 1 / c.J2WH; % 1 watt hour to joules
    
    %% Speed
    %MPH :miles per hour (mph)
    %KPH :Kilometers per hour (kph)
    c.MPH2KPH = 1.609344; % 1 mile per hour to kilometers per hour
    c.KPH2MPH = 1 / c.MPH2KPH; % 1 kilometer per hour to miles per hour
    
    % Film Cooling Integration
    % MPS: meters per second [m/s]
    % FPS: feet per second [ft/s]
    c.MPS2FPS = 3.281;          % Feet per second per mps
    c.FPS2MPS = 1 / c.MPS2FPS;  % Meters per second per fps
    c.MPH2FPS = 1.467;          % Feet per second per miph
    c.FPS2MPH = 1 / c.MPH2FPS;  % Miles per hour per fps
    
    %% Volume
    % CF: cubic feet [ft^3]
    % CI: cubic inches [in^3]
    % CM: cubic meters [m^3]
    c.CF2CI = 1728;             % Cubic inches per cubic feet
    c.CI2CF = 1 / c.CF2CI;      % Cubic feet per cubic inches
    c.CM2CI = 61024;            % Cubic inches per cubic meter
    c.CI2CM = 1 / c.CM2CI;      % Cubic meters per cubic inch
    c.CM2CF = 35.315;           % Cubic feet per cubic meter
    c.CF2CM = 1 / c.CM2CF;      % Cubic meters per cubic feet
    
    %% Specific Heat Conversions
    % BTULF: BTU per pound degree Fahrenheit [BTU/(lb * degF)]
    % JKGK: Joules per kilogram Kelvin [J/(kg * K)]
    c.BTULF2JKGK = 4186.7982;        % Joules per kilogram Kelvin per BTU/(lb * degF)
    c.JKGK2BTULF = 1 / c.BTULF2JKGK; % BTU per pound Fahrenheit per J/kgK
    
    %% THERMAL CONDUCTIVITY CONVERSIONS:
    % BTUSFI: BTU per second degree Fahrenheit inch [BTU/(s * in * degF)]
    % WPMK: Watts per meter Kelvin [W/(m * K)]
    c.WPMK2BTUSFI = 0.27733887194;     % BTU per second Fahrenheit inch per wpmk
    c.BTUSFI2WPMK = 1 / c.WPMK2BTUSFI; % Watts per meter kelvin per btusfi
    
    %% ENTHALPY CONVERSIONS:
    % JPK: joules per kilogram [J/kg]
    % BTUL: BTU per pound [BTU/lb]
    c.JPK2BTUL = 0.0004299226;   % BTU per pound per joules per kilogram
    c.BTUL2JPK = 1 / c.JPK2BTUL; % Joules per kilogram per BTU per pound
    
    %% VISCOSITY CONVERSIONS:
    % PS: pascal seconds [Pa*s]
    % LIS: pound per second inch [lb/(s * in)]
    % MPS: mPa seconds [mPa*s]
    c.PS2LIS = 0.05599741459;     % Pound per inch seconds per pascal seconds
    c.LIS2PS = 1 / c.PS2LIS;      % Pascal seconds per pound per inch seconds
    c.MPS2LIS = 0.00005599741459; % Pound per inch seconds per mpascal seconds
    c.LIS2MPS = 1 / c.MPS2LIS;    % mPascal seconds per pound per inch seconds

end