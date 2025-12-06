% ENGR492 Intro to Aircraft Design
% Assignment 4
% Jasper Palmer

clc;
clear;

% DATUM: Nose of a/c
m_f = 7600; % Fuel weight, kg
m_e = 10500; % Empty weight, kg
m_p = 4500; % Payload/crew weight, kg
m_0 = m_e + m_f + m_p; % MTOW, kg

WS = 3800; % Wing loading, N/m^2
TW = 0.35; % Trust to weight
x_LE = 7.25; % Root leading edge distance from datum, m
x_ng = 3.5; % Nose gear distance from datum, m
b = 21; % Wingspan, m
x_mg = 10.56; % Main gear distance from datum, m
A = 7.6; % Aspect ratio
T = 4.78; % gear track, m
C_r = 4.4; % Root chord, m
B = 7.06; % Main gear wheelbase, m
lambda = 0.25; % Taper ratio, []
phi_ot = 50; % Overturn angle, deg
Lambda_LE = 26; % Leading edge sweep, deg
g = 9.81;

% Requirements
alpha_to = 15; % Fuselage angle of attack at takeoff, deg
H_c = 0.3; % Fuselage ground clearance at rotation
a_landing = -3.5; % Landing decelleration, m/s^2
a_to = 4.5; % Takeoff acceleration, m/s^2
V_taxi = 20*0.5144; % Maximum ground speed, m/s
xbar_ac = 0.25;

%% CG Calculations

H_cg = T/(2*tand(phi_ot));

MAC = C_r*2/3*(1+lambda+lambda^2)/(1+lambda);
Ybar = b/6*(1+2*lambda)/(1+lambda);
x_mac = x_LE + Ybar*tand(Lambda_LE);

% Most aft CG
x_ac = x_mac + xbar_ac*MAC;
x_cg_aft = x_ac - 0.001;
h_aft = (x_cg_aft - x_mac)/MAC;

% Most fore CG
x_cg_fore = x_ng + B*0.8 + abs(a_landing)*H_cg/g;
h_fore = (x_cg_fore - x_mac)/MAC;


%% Gear Calculations

% Minimum gear height for desired clearance
AB = 14.338-10.560; % distance from main gear to fuselage upsweep, m
CD = H_c/cosd(alpha_to);
BD = AB*tand(alpha_to);
H_g = CD + BD;

% Overturn angle

r = 2*H_cg*V_taxi^2/(T*g);


%% Output

table([MAC; x_mac; x_cg_fore; x_cg_aft; h_fore*100; h_aft*100; H_g; r], 'RowNames',{'MAC, m','x_mac, m','x_cg_fore, m','x_cg_aft, m','h_fore, %','h_aft, %','H_g, m','r, m'},'VariableNames',{'Result'})