% ENGR 493 Intro to Aircraft Design
% Project Phase 2: Wing Design
% Jasper Palmer

%% Project Requirements
% - Fix phase 1
% - Size wing for stall conditions
% - Optimization studies to maximize (L/D)cr
% - Report wing parameters (span, Sref, A, airfoil, MAC, t/c, LE sweep, TE
%   sweep, taper, tip geometry

% Todos
% - table to show iterations (iterate wing -> iterate W0 -> iterate wing -> ...
% - Up to module 2, p.42 if using aero-solver, otherwise to p.43
% - Requires at least 2-3 iterations

%% Functions

Atmosphere = table2array(readtable("Atmosphere.csv"));
Vf = @(M, A) M*interp1(Atmosphere(:,1),Atmosphere(:,4), A);
nuf = @(A) interp1(Atmosphere(:,1),Atmosphere(:,3),A);
qf = @(V, rho) 0.5*rho*V^2;
Mf = @(V, A) interp1(Atmosphere(:,1),Atmosphere(:,4),A)/V;

%% Variables

W_0g = 40000; % Guess W0, kg
C_Lmaxg = 2.4; % Guess C_Lmax
err = 1; % Error between guess MTOW and calculated MTOW for satisfied solver

K_rf = 1.06; % reserve fuel ratio

g = 9.81; % gravitational constant, m/s^2

AIM120C = 348/2.2; % AIM-120C weight, kg
AIM9X = 186/2.2; % AIM-9X weight, kg
MK83JDAM = 1000/2.2;% JDAM weight, kg
W_p = 2500/2.2; % Payload weight, kg
W_c = 250; % Crew weight, kg
Alt_cr = 10000; % cruise altitude, m
Alt_lol = 3000; % Loiter-to-land altitude, m
M_cr = 0.85; % cruise mach
V_cr = Vf(M_cr, Alt_cr); % cruise speed, m/s
t_lol = 30*60; % Loiter-to-landing endurance (holding time)
C_lol = 21*g*1e-6; % Loter-to-landing consumption, 1/s
C_cr = 42*g*1e-6;
rho_landing = 1.1287; % Tropical day sea level air density, kg/m^3
rho_0 = 1.225; % Sea-level standard day air density, kg/m^3
rhof = @(A) rho_0*interp1(Atmosphere(:,1),Atmosphere(:,2),A);
rho_lol = rhof(Alt_lol);
nu_lol = nuf(Alt_lol);
V_lol = 128; % Loiter-to-land speed (assuming flying at maximum civil airspeed below 10000'), m/s
q_lol = qf(V_lol, rho_lol);
rho_cr = rhof(Alt_cr);
nu_cr = nuf(Alt_cr);
q_cr = qf(V_cr, rho_cr);
M_lol = Mf(V_lol, Alt_lol);

A_top = 22; % Fuselage top area based on CAD, m^2
A_side = 14; % Fuselage side area based on CAD, m^2
d = 1.7; % Diameter of fuselage based on CAD, m
L_fus = 18; % Length of fuselage, m
A_max = 3; % Maximum cross-sectional area, m^2

% DOGFIGHT
% Profile: T/O -> Climb -> Cruise for combat radius -> Accelerate + descent
% -> Loiter for dogfight -> Climb -> Cruise for combat radius -> Descent
% -> Loiter -> Landing

dogfight.C = 34*g*1e-6; % consumption at max thrust and 10000', 1/s
dogfight.E = 5*60; % dogfight endurance, s
dogfight.R = 1000*1852; % combat radius, m
dogfight.W_d = 6*AIM120C+2*AIM9X; % drop weight, kg
dogfight.M_loi = 0.75; % Best turn mach, rough estimation
dogfight.V_loi = Vf(dogfight.M_loi, 3048); % loiter speed, m/s
dogfight.rho_loi = rhof(3048);
dogfight.nu_loi = nuf(3048);
dogfight.q_loi = qf(dogfight.V_loi, dogfight.rho_loi);

% STRIKE
% Profile: T/O -> Climb -> Cruise for combat radius -> Descent -> Loiter
% for 50nm at intermediate thrust -> Climb -> Cruise -> Descent -> Loiter
% -> Landing

strike.C = 53*g*1e-6; % consumption at max non-afterburner thrust and sea level, 1/s
strike.M_loi = 0.9; % strike loiter mach
strike.V_loi = strike.M_loi*340.3; % loiter speed, m/s
strike.E = 50*1850/strike.V_loi; % endurance, s
strike.R = 1000*1852; % combat radius, m
strike.W_d = 4*MK83JDAM+2*AIM9X; % drop weight, kg
strike.rho_loi = rhof(0);
strike.nu_loi = nuf(0);

%% Solution

% DOGFIGHT
% rho_landing, V_end, V_WoD, V_TO, DeltaV_thrust, C_Lmax)

tc_root = 0.139; % Thickness/chord of airfoil at tip
tc_tip = 0.06; % thickness/chord of airfoil at tip
xc_m = 0.385; % Normalized x location of maximum thickness, average between root and tip

deltay_root = .04363/.00706; % t(x=0.06*C)/t(x=0.0015*C) from airfoil
deltay_tip = .01888/.00944;

eps = 9e9;
it = 0;

dx = 10;

b = linspace(5, 20, dx);
lambda = linspace(0, 1, dx);
Lambda_LE = linspace(0, 60, dx);
% range = linspace(1000*1852, 2000*1852, dx);
% endurance = linspace(3*60, 8*60, dx);
% W_d = linspace(4*AIM120C+2*AIM9X, 8*AIM120C+4*AIM9X, dx);

for i = 1:length(b)
    for j = 1:length(lambda)
        for k = 1:length(Lambda_LE)
            % for w = 1:length(range)
                % for x = 1:length(endurance)
                    % for y = 1:length(W_d)

                        while eps > err && it < 10000
                            

                            S_ref = performance.wingloading(W_0g, rho_landing, C_Lmaxg); % Wing loading
                            [~, ~, ~, ~, ~, Ld_loi, ~, ~] = aerodynamics.wing(W_0g, dogfight.V_loi, dogfight.nu_loi, dogfight.q_loi, S_ref, b(i), lambda(j), Lambda_LE(k), dogfight.M_loi, tc_root, tc_tip, A_top, A_side, xc_m, deltay_root, deltay_tip, L_fus, A_max, d); % Loiter aerodynamics
                            [A, MAC, Lambda_TE, Lambda_025c, S_canard, Ld_cr, C_Lmax, alpha_stall] = aerodynamics.wing(W_0g, V_cr, nu_cr, q_cr, S_ref, b(i), lambda(j), Lambda_LE(k), M_cr, tc_root, tc_tip, A_top, A_side, xc_m, deltay_root, deltay_tip, L_fus, A_max, d); % Cruise aerodynamics
                            [~, ~, ~, ~, ~, Ld_lol, ~, ~] = aerodynamics.wing(W_0g, V_lol, nu_lol, q_lol, S_ref, b(i), lambda(j), Lambda_LE(k), M_lol, tc_root, tc_tip, A_top, A_side, xc_m, deltay_root, deltay_tip, L_fus, A_max, d); % Loiter-to-landing aerodynamics
                            Ld_loi
                            Ld_lol
                            Ld_cr
                            [W_0, W_f, W_fi] = weight.weights(W_0g, W_p, dogfight.W_d, W_c, K_rf, Ld_loi, dogfight.C, dogfight.M_loi, dogfight.E, Ld_cr, C_cr, M_cr, V_cr, dogfight.R, C_lol, Ld_lol, t_lol); % Weight sizing
                            eps = abs(W_0g-W_0);
                            W_0g = W_0;

                        end

                    % end
                % end
            % end
        end
    end
end

% score = score.score(Ld_cr, AoA_stall);

% Dogfight_Output = table(b, S_ref, A, MAC, tc, Lambda_LE, Lambda_TE, lambda, Ld_loi, Ld_cr, Ld_lol, W_0, W_f, 'VariableNames',{'b','S_ref','A','MAC','tc','Λ_LE','Λ_TE','λ','(L/D)_loi','(L/D)_cr','(L/D)_lol','W_0','W_f'})
% figure
% plot(W_fi)

% STRIKE

%{
WS_stall = wingloading();
S_ref = WS_stall*W0g;
[b, A, MAC, tc, Lambda_LE, Lambda_TE, lambda, Ld_loi, Ld_cr, Ld_lol] = wing();
[W_0, W_f, W_fi] = weights(W_0g, W_p, strike.W_d, W_c, W_to, W_de, W_la, K_rf, Ld_loi, strike.C, strike.M_loi, strike.E, Ld_cr, C_cr, M_cr, V_cr, strike.R, C_lol, Ld_lol, t_lol);
Strike_Output = table(b, S_ref, A, MAC, tc, Lambda_LE, Lambda_TE, lambda, Ld_loi, Ld_cr, Ld_lol, W_0, W_f, 'VariableNames',{'b','S_ref','A','MAC','tc','Λ_LE','Λ_TE','λ','(L/D)_loi','(L/D)_cr','(L/D)_lol','W_0','W_f'})
figure
plot(W_fi)
%}
