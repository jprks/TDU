function varargout = tdu % main function
clear all;
clc;
%% Universal Constants
universal_gas_constant= 8.3144598; % [J/(mol*K)]
standard_gravity = 9.81; % [m/s^2]
unitless='[-]';

%% Gas Properties of Propellant
propellant_name = 'Nitrogen';
specific_heat_ratio = 1.4; % 1.4 for Nitrogen = 1.67 for Xenon
molar_mass = .028; % .014006 for Nitrogen - 0.131293 for Xenon
molar_mass_units = '[kg/mol]';

%% Nozzle Geometry
exit_radius = .014160; % radius at nozzle exit
throat_radius = .01; % radius at nozzle throat
length_units = '[m]';
conical_half_angle = 15; % half angle of conical nozzle, 15 degrees is optimal
angle_units = '[deg]';

%% Chamber Conditions
chamber_temperature = 273; % 273 K = 0C
temperature_units = '[K]';
pressure_units = '[Pa]';

%% Initial Tank Conditions
spherical_tank_radius = 0.025;
length = 0.05; % The length of the tank towards the nozzle
width = 0.1;
height = 0.1;
length_units = '[m]';
% spherical_tank_volume = 4/3*(pi*tank_radius^3);
% torus_tank_volume = 1/4*pi^2*(inner radius+outer_radius)*(outer_radius-inner_radius)^2;
% rectangular_tank_volume = length*width*height;
tank_volume_units = '[m^3]';
tank_radius_units = '[m]';
% Using estimated propellant masses we can calculate a target chamber
% pressure
tank_radius = spherical_tank_radius;
theoretical_propellant_mass = 0.460;
propellant_mass_units = '[Kg]';
initial_system_mass = 2.66;
system_mass_units = '[kg]';
% theoretical_chamber_pressure = (3*theoretical_propellant_mass*universal_gas_constant*chamber_temperature)/(4*pi*molar_mass*tank_radius^3);
% theoretical_chamber_pressure_psi = theoretical_chamber_pressure/6895.75729;
chamber_pressure_psi_units = '[Psi]';

%% Math (the fun part)

%           NOZZLE NOMENCLATURE
%   ********************************
%
%              /-
%             /
%     ===\---/          c = chamber
%     (c) (t)  (e)      t = throat
%     ===/---\          e = exit
%             \
%              \-
%   ********************************

specific_gas_constant = universal_gas_constant/molar_mass;
specific_gas_constant_units = '[J/kg K]';
exit_area = pi*exit_radius^2;
throat_area = pi*throat_radius^2;
area_units = '[m^2]';
length = (exit_radius-throat_radius)/tan(deg2rad(conical_half_angle)); % [m]
length_units = '[m]';



%% Exit conditions
%   it's much easier to define the flow characteristics at the exit in
%   terms of ratios for now. We'll resolve them into actual useful values
%   later.
%   these are ratios between the throat parameter (numerator) and the exit
%   parameter (denominator).
%   For example, area_ratio = throat_area/exit_area

%(already know this one)    area_ratio = exit_mach(((specific_heat_ratio+1)/2)/(1+(specific_heat_ratio-1)/2*exit_mach^2))^((specific_heat_ratio+1)/(2*(specific_heat_ratio-1)));
exit_mach = solve_mach(exit_area,throat_area,specific_heat_ratio);

temperature_ratio = 1+((specific_heat_ratio-1)/2)*exit_mach^2;
pressure_ratio = temperature_ratio^((specific_heat_ratio-1)/specific_heat_ratio);

exit_temperature = chamber_temperature/temperature_ratio;

exit_velocity = sqrt((2*specific_heat_ratio*specific_gas_constant*(chamber_temperature))/(specific_heat_ratio-1)*(1-1/(1+(specific_heat_ratio-1)/2*exit_mach^2)));
velocity_units = '[m/s]';

%% Final Tank Conditions

system_delta_velocity = 122;
system_delta_velocity_units = '[m/s]';
specific_impulse_ve = exit_velocity/standard_gravity;
final_system_mass = initial_system_mass/(exp(system_delta_velocity/(specific_impulse_ve*standard_gravity)));
propellant_mass = initial_system_mass - final_system_mass;
mass_ratio = final_system_mass/initial_system_mass; % Dimensionless
chamber_pressure = (3*propellant_mass*universal_gas_constant*chamber_temperature)/(4*pi*molar_mass*tank_radius^3);
exit_pressure = chamber_pressure/pressure_ratio;

%% Throat conditions

% We are designing a nozzle that "chokes" the fluid flow at the throat,
% so we assume that at the throat, Mach number is isentropically
% brought to the boundary between supersonic and subsonic flow (Mach=1)
throat_temperature = chamber_temperature*(2/(specific_heat_ratio+1));
throat_pressure = chamber_pressure*(2/(specific_heat_ratio+1))^(specific_heat_ratio/(specific_heat_ratio-1));
throat_velocity = sqrt(specific_heat_ratio*specific_gas_constant*throat_temperature);
% use nozzle areas to find Mach number at the exit via the Newton-Raphson method
%   since
%       mass_flow_rate = density*speed*cross_sectional_area
%   and from the ideal gas law, PV=RT
throat_density = throat_pressure/(specific_gas_constant*throat_temperature); %[kg/m^3]
throat_density_units = '[kg/m^3]';
%   and at the throat Mach = 1 so flow rate (speed of flow) is the fluid's speed of sound at the throat
throat_flowrate = sqrt(specific_heat_ratio*specific_gas_constant*throat_temperature);
mass_flowrate = throat_density*throat_flowrate*throat_area;
%     mass_flowrate = (specific_heat_ratio*chamber_pressure*exit_area*exit_mach)/(specific_heat_ratio*specific_gas_constant*chamber_temperature)^(1/2)*(1+(specific_heat_ratio-1)*(exit_mach)^2/2)^((2-specific_heat_ratio)/2);
mass_flowrate_units = '[kg/s]';

%% Thrust and Specific Impulse

thrust = mass_flowrate*exit_velocity + exit_pressure*exit_area;
force_units = '[N]';
specific_impulse_thrust = thrust/(mass_flowrate*standard_gravity);
%     specific_impulse = exit_velocity/standard_gravity;
isp_units='[s]';

%% System Calculations

system_velocity = 122;
initial_mass = 2.66;
initial_mass_units = '[kg]';
%     calculated_propellant_mass = initial_mass*(exp(system_velocity/(specific_impulse*standard_gravity)-1));
%     mass_ratio = exp(system_velocity/(specific_impulse*standard_gravity));
mass_ratio_units = '[-]';
system_velocity = standard_gravity*specific_impulse_ve*log(mass_ratio);
system_velocity_units = '[m/s]';
burn_time = (propellant_mass*system_delta_velocity)./(thrust);
burn_time_units = '[s]';

%% format & display outputs

linedivider='------------';
result =  {'Propellant','',propellant_name;
    linedivider,'','';
    'Specific heat ratio', specific_heat_ratio, unitless;
    'Molar mass', molar_mass, molar_mass_units;
    'Specific gas constant',specific_gas_constant,specific_gas_constant_units;
    linedivider,'','';
    'Chamber temperature', chamber_temperature, temperature_units;
    'Chamber pressure', chamber_pressure, pressure_units;
%     'Chamber pressure (Psi)',theoretical_chamber_pressure_psi,chamber_pressure_psi_units;
%     'Theoretical Chamber Pressure', theoretical_chamber_pressure,chamber_pressure_psi_units;
    'Exit radius', exit_radius, length_units;
    'Throat radius', throat_radius, length_units;
    'Nozzle Length',length,length_units;
    'Half-angle',conical_half_angle, angle_units;
    linedivider,'','';
    'Propellant Mass',propellant_mass,propellant_mass_units;
    'Theoretical Propellant Mass',theoretical_propellant_mass,propellant_mass_units;
    'Tank Radius',tank_radius,tank_radius_units;
    'Mass Ratio',mass_ratio,mass_ratio_units;
    linedivider,'','';
    'Length',length,length_units;
    'Exit area',exit_area,area_units;
    'Throat area',throat_area,area_units;
    'Throat temperature',throat_temperature,temperature_units;
    'Throat pressure',throat_pressure,pressure_units;
    'Mass flow rate',mass_flowrate,mass_flowrate_units;
    'Throat flowrate',throat_flowrate,mass_flowrate_units;
    linedivider,'','';
    'Exit temperature',exit_temperature,temperature_units;
    'Exit pressure',exit_pressure,pressure_units;
    linedivider,'','';
    'Exhaust velocity',exit_velocity,velocity_units;
    'Thrust',thrust,force_units;
    'Specific impulse (Ve)',specific_impulse_ve,isp_units;
    'Specific impulse (Thrust)',specific_impulse_thrust,isp_units
    'System Velocity',system_velocity,system_velocity_units;
    'Burn Time',burn_time,burn_time_units;
    linedivider,'','';
    'Exit Mach',exit_mach,unitless;
    'A/At',exit_area/throat_area,unitless;
    'T/Tc',exit_temperature/chamber_temperature,unitless;
    'P/Pc',exit_pressure/chamber_pressure,unitless;
    'v/at',exit_velocity/throat_velocity,unitless;
    };

display(result);
end

function Mach = solve_mach(A,At,k)
%   solve Mach number from area ratio by Newton-Raphson Method. (assume
%   supersonic)
%   https://www.grc.nasa.gov/WWW/winddocs/utilities/b4wind_guide/mach.html
P = 2/(k+1);
Q = 1-P;
R = (A/At).^((2*Q)/P);
a = Q.^(1/P);
r = (R-1)/(2*a);
X = 1/((1+r)+sqrt(r*(r+2)));  % initial guess
diff = 1;  % initalize termination criteria
while abs(diff) > .0001
    F = (P*X+Q).^(1/P)-R*X;
    dF = (P*X+Q).^((1/P)-1)-R;
    Xnew = X - F/dF;
    diff = Xnew - X;
    X = Xnew;
end
Mach = 1/sqrt(X);
end

function display(result)
[n,~]=size(result);
for i = 1:n
    fprintf('\n%24s\t%12f\t%s',result{i,:});
end
end
