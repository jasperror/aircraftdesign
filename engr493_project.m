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

%% Variables

W_0g = 40000; % Guess W0, kg
C_Lmaxg = 2.4; % Guess C_Lmax
err = 0.001; % Maximum root-mean-square of relative error of various guess parameters

K_rf = 1.06; % reserve fuel ratio

AIM120C = 348/2.2; % AIM-120C weight, kg
AIM9X = 186/2.2; % AIM-9X weight, kg
MK83JDAM = 1000/2.2;% JDAM weight, kg
W_p = 2500/2.2; % Payload weight, kg
W_c = 250; % Crew weight, kg
Alt_cr = 10000; % cruise altitude, m
Alt_lol = 3000; % Loiter-to-land altitude, m
M_cr = 0.85; % cruise mach
V_cr = atmosphere.V(M_cr, Alt_cr); % cruise speed, m/s
t_lol = 30*60; % Loiter-to-landing endurance (holding time)
C_lol = 21*const.g*1e-6; % Loter-to-landing consumption, 1/s
C_cr = 42*const.g*1e-6;
rho_landing = 1.1287; % Tropical day sea level air density, kg/m^3

rho_lol = atmosphere.rho(Alt_lol);
nu_lol = atmosphere.nu(Alt_lol);

V_lol = 128; % Loiter-to-land speed (assuming flying at maximum civil airspeed below 10000'), m/s
q_lol = atmosphere.q(V_lol, rho_lol);
rho_cr = atmosphere.rho(Alt_cr);
nu_cr = atmosphere.nu(Alt_cr);
q_cr = atmosphere.q(V_cr, rho_cr);
M_lol = atmosphere.M(V_lol, Alt_lol);

A_top = 22; % Fuselage top area based on CAD, m^2
A_side = 14; % Fuselage side area based on CAD, m^2
d = 1.7; % Diameter of fuselage based on CAD, m
L_fus = 18; % Length of fuselage, m
A_max = 3; % Maximum cross-sectional area, m^2

% DOGFIGHT
% Profile: T/O -> Climb -> Cruise for combat radius -> Accelerate + descent
% -> Loiter for dogfight -> Climb -> Cruise for combat radius -> Descent
% -> Loiter -> Landing

dogfight.C = 34*const.g*1e-6; % consumption at max thrust and 10000', 1/s
dogfight.E = 5*60; % dogfight endurance, s
dogfight.R = 1000*1852; % combat radius, m
dogfight.W_d = 6*AIM120C+2*AIM9X; % drop weight, kg
dogfight.M_loi = 0.75; % Best turn mach, rough estimation
dogfight.V_loi = atmosphere.V(dogfight.M_loi, 3048); % loiter speed, m/s
dogfight.rho_loi = atmosphere.rho(3048);
dogfight.nu_loi = atmosphere.nu(3048);
dogfight.q_loi = atmosphere.q(dogfight.V_loi, dogfight.rho_loi);

% STRIKE
% Profile: T/O -> Climb -> Cruise for combat radius -> Descent -> Loiter
% for 50nm at intermediate thrust -> Climb -> Cruise -> Descent -> Loiter
% -> Landing

strike.C = 53*const.g*1e-6; % consumption at max non-afterburner thrust and sea level, 1/s
strike.M_loi = 0.9; % strike loiter mach
strike.V_loi = strike.M_loi*340.3; % loiter speed, m/s
strike.E = 50*1850/strike.V_loi; % endurance, s
strike.R = 1000*1852; % combat radius, m
strike.W_d = 4*MK83JDAM+2*AIM9X; % drop weight, kg
strike.rho_loi = atmosphere.rho(0);
strike.nu_loi = atmosphere.nu(0);

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
C_Lmax_0
W_0l

for i = 1:length(b)
    for j = 1:length(lambda)
        for k = 1:length(Lambda_LE)
            % for w = 1:length(range)
            % for x = 1:length(endurance)
            % for y = 1:length(W_d)
            % dash cr
            % dash sl

%{
M,... Mach number
A_max,... Maximum cross-sectional area, m^2
l,... Body length, m
Lambda_LE,... Leading edge sweep, RAD
Lambda_025c,... Quarter-chord sweep, DEG
S_ref... Reference area, m^2
%}

            M_cr_guess = 0.004*aerodynamics.Sweeptf(Lambda_LE(), 0, 25, A_guess, lambda())+0.64;
            C_D0_cr_guess = aerodynamics.C_D0_f(M_cr_guess, atmosphere.V(M_cr_guess, h_cr), L_fus, atmosphere.nu(h_cr), aerodynamics.Swetexp(tc_root, tc_tip, S_ref_guess, lambda, b, A_top, A_side), S_ref_guess) + aerodynamics.C_Dw_f(M_cr_guess, A_max, L_fus, Lambda_LE(), aerodynamics.Sweeptf(Lambda_LE(), 0, 25, A_guess, lambda()), S_ref_guess);
            C_D0_maxV_cr_guess = aerodynamics.C_D0_f(M_dash_cr(), atmosphere.V(M_dash_cr(), 20e3*0.3048), L_fus, atmosphere.nu(20e3*0.3048), aerodynamics.Swetexp(tc_root, tc_tip, S_ref_guess, lambda, b, A_top, A_side), S_ref_guess) + aerodynamics.C_Dw_f(M_dash_cr(), A_max, L_fus, Lambda_LE(), aerodynamics.Sweeptf(Lambda_LE(), 0, 25, A_guess, lambda()), S_ref_guess);
            C_D0_maxV_sl_guess = aerodynamics.C_D0_f(M_dash_sl(), atmosphere.V(M_dash_sl(), 0*0.3048), L_fus, atmosphere.nu(0e3*0.3048), aerodynamics.Swetexp(tc_root, tc_tip, S_ref_guess, lambda, b, A_top, A_side), S_ref_guess);
            C_Lmax_0_guess = aerodynamics.maxlift();
            %{
                C_L_max,... Max lift coefficient for wing
                alpha_C_L_max... AoA of max lift, rad
                ] = maxlift(...
                deltay,... t(0.06*C)-t(0.0015*C) from airfoil
                alpha_ZL,... Zero lift AoA, rad
                C_l_max,... Max lift coefficient for airfoil
                a,... C_L_alpha (slope of CL vs alpha curve for wing), rad^-1
                Lambda_LE,... Leading edge sweep, rad
                A,... Aspect ratio
                S_ref,... Reference area, m^2
                lambda... Taper ratio
            %}
            W_0l_guess = ;
            W_0cr_guess = ;
            W_0midmission_guess = ;

            while rmsre > err && it < 1e6
                % Estimations for Performance
                % Performance Calculation
                [S_ref, T_0, WS_0 , TW_0, V_maxturn] = performance.fighter(psidot(), N_z/1.5, M_dash_cr(), M_dash_sl(), A_guess, e_guess, W_0g, C_Lmax_0_guess, W_0l_guess, C_D0_cr_guess, W_0cr_guess, h_cr, W_0midmission_guess, d, lambda(), b(), aerodynamics.Sweeptf(Lambda_LE(), 0, 25, A_guess, lambda()), C_D0_maxV_cr_guess, C_D0_maxV_sl_guess)
                % Aerodynamic Calculations
                [~, ~, ~, ~, ~, Ld_loi, ~, ~] = aerodynamics.wing(W_0g, dogfight.V_loi, dogfight.nu_loi, dogfight.q_loi, S_ref, b(i), lambda(j), Lambda_LE(k), dogfight.M_loi, tc_root, tc_tip, A_top, A_side, xc_m, deltay_root, deltay_tip, L_fus, A_max, d); % Loiter aerodynamics
                [A, MAC, Lambda_TE, Lambda_025c, S_canard, Ld_cr, C_Lmax, alpha_stall] = aerodynamics.wing(W_0g, V_cr, nu_cr, q_cr, S_ref, b(i), lambda(j), Lambda_LE(k), M_cr, tc_root, tc_tip, A_top, A_side, xc_m, deltay_root, deltay_tip, L_fus, A_max, d); % Cruise aerodynamics
                [~, ~, ~, ~, ~, Ld_lol, ~, ~] = aerodynamics.wing(W_0g, V_lol, nu_lol, q_lol, S_ref, b(i), lambda(j), Lambda_LE(k), M_lol, tc_root, tc_tip, A_top, A_side, xc_m, deltay_root, deltay_tip, L_fus, A_max, d); % Loiter-to-landing aerodynamics
                Ld_loi
                Ld_lol
                Ld_cr
                % Weight Calculations
                [W_0, W_f, W_fi] = weight.weights(W_0g, W_p, dogfight.W_d, W_c, K_rf, Ld_loi, dogfight.C, dogfight.M_loi, dogfight.E, Ld_cr, C_cr, M_cr, V_cr, dogfight.R, C_lol, Ld_lol, t_lol); % Weight sizing
                % Calculate Root Mean Square Relative Error
                rmsre = sqrt(mean([ ...
                    (W_0g-W_0)/W_0, ...
                    (A_guess-A)/A, ...
                    (e_guess-e)/e, ...
                    (C_Lmax_0_guess-C_Lmax_0)/C_Lmax_0, ...
                    (W_0l_guess-W_0l)/W_0l, ...
                    (C_D0_cr_guess-C_D0_cr)/C_D0_cr, ...
                    (W_0cr_guess-W_0cr)/W_0cr, ...
                    (W_0midmission_guess-W_0mm)/W_0mm, ...
                    (S_ref_guess-S_ref)/S_ref,...
                    ].^2));

                %{
A_guess
S_ref_guess
e_guess
W_0g
C_Lmax_0_guess?
W_0l_guess
W_0cr_guess
W_0midmission_guess

                %}
                % Update Parameters
                W_0g = W_0;
                A_guess = A;
                e_guess = e;
                W_0l_guess = W_0l;
                

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
