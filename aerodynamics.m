% ENGR 493 Intro to Aircraft Design
% Project Phase 2: Wing Design - Aerodynamics of the wing
% Jasper Palmer

% todo -- convert all sweeps to radians or degrees? same with AoA?
classdef aerodynamics
    properties
        DC_l_max_TE = 1.3; % Assuming slotted flap
        DC_l_max_LE = 0.3; % Assuming leading-edge flap
    end
    methods (Static)
        function Re = Re_f(V, L, nu) % Reynolds number
            Re = V*L/nu;
        end

        function circ = circf(Lp, rho_inf, V_inf) % Circulation (Gamma)
            circ = Lp/(rho_inf*V_inf);
        end

        % function delta = deltaf(lambda, A) % Simplified correction factor
        %
        %     A4 = [0.08, 0.02, 0.002, 0, 0.001, 0.003, 0.007, 0.016, 0.023, 0.032, 0.039];
        %     A6 = [0.11,0.03, 0.01, 0.005, 0.0051, 0.01, 0.018, 0.025, 0.035, 0.048, 0.065];
        %     A8 = [0.123, 0.038, 0.012, 0.007, 0.007, 0.015, 0.024, 0.035, 0.05, 0.065, 0.085];
        %     A10 = [0.16, 0.054, 0.022, 0.01, 0.015, 0.023, 0.033, 0.047, 0.064, 0.085, 0.105];
        %     % Ax = 4:2:10;
        %     x = 0:0.1:1;
        %
        %     if A < 4 || A > 10
        %         disp('Error: outside of delta function limits')
        %         delta = NaN;
        %         return
        %     elseif lambda < 0 }} lambda > 1
        %         disp('Error: outside of delta function limits')
        %         delta = NaN;
        %         return
        %     end
        %
        %     if A < 6
        %         lower = interp1(x, A4, lambda);
        %         upper = interp1(x, A6, lambda);
        %     elseif A < 8
        %         lower = interp1(x, A6, lambda);
        %         upper = interp1(x, A8, lambda);
        %     else
        %         lower = interp1(x, A8, lambda);
        %         upper = interp1(x, A10, lambda);
        %     end
        %
        %     delta = interp1(A, [lower, upper]);
        % end

        function e = ef(... Oswald efficiency
                A,... Aspect ratio
                Lambda_LE... Leading edge sweep, deg
                )
            el = 1.78*(1-0.045*A^0.68)-0.64;
            eu = 4.61*(1-0.045*A^0.68)*cosd(Lambda_LE)^0.15-3.1;
            e = interp1([0, 30], [el, eu], Lambda_LE);
        end

        function delta = deltaf(e)
            delta = 1-1/e;
        end

        function C_Di = C_Di_f(... Induced Drag Coefficient
                M,... Mach number
                A,... Aspect ratio
                Lambda_LE,... Leading edge sweep, deg
                C_L... Lift coefficient
                )
            % Simplified Expression: C_Di_f = @(C_L, A, delta) C_L^2/(pi*A*(1/(1+delta))); % Induced drag
            if M < 1
                K = 1/(pi*A*e);
            else
                K = A*(M^2-1)*cosd(Lambda_LE)/(4*A*sqrt(M^2-1)-2);
            end
            C_Di = K*C_L^2;
        end

        function C_f = C_f_f(... Friction Coefficient
                M,... Mach number
                V_inf,... Freestream speed, m/s
                L,... Length of object, m
                nu_inf... Kinematic viscosity, m^2/s
                )

            k = 1.015e-5; % Surface roughness, m

            if M < 1
                Recut = 38.21*(L/k)^1.053;
            else
                Recut = 44.62*(L/k)^1.053*M^1.16;
            end
            Re = min(aerodynamics.Re_f(V_inf, L, nu_inf), Recut);

            Cflam = 1.328/Re^0.5;
            Cfturb = 0.074/Re^0.2;

            % lamfrac = [0, 10; 10, 35; 25, 50; 35, 70; 0, 0; 20, 20; 5, 10; 25, 50; 50, 80; 0, 0];
            lamfrac = 0; % Zero laminar flow for military aircraft with camouglage

            C_f = Cflam*lamfrac + Cfturb*(1-lamfrac);
        end

        function FF = FF_fuselage_canopy(l, A_max) % Friction factor of fuselage/canopy
            FF = 0.9+5/l*sqrt(4/pi*A_max)+l/sqrt(4/pi*A_max)/400;
        end

        function FF = FF_tail(tc, xcm, M, Lambda_m) % Friction factor of tail
            FF = (1+0.6/xcm*tc+100*tc^4)*(1.34*M^0.18*cosd(Lambda_m)^0.28);
        end

        function FF = FF_store(l, A_max) % Friction factor of stores
            FF = 1 + 0.35/l*sqrt(4/pi*A_max);
        end

        function C_D0 = C_D0_f(...Parasitic Drag Coefficient (excl. wave drag)
                M,... Mach number
                V_inf,... Freestream speed, m/s
                L_fus,... Length of fuselage
                nu_inf,... Kinematic viscosity, m^2/s
                S_wet,... Total wetted area, m^2
                S_ref... Reference area, m^2
                )
            % Components: Fuselage, external stores, cockpit, vertical stabilizer

            % ASSUMPTIONS (FIX ME!)
            % Assuming cockpit is 20% the length of the fuselage, vertical
            % stabilizer is 30% of the length of the fuselage, assuming (big
            % assumption) that the diameter of the frontal area of the fuselage is
            % 15% of the length of the fuselage, cockpit frontal area is 5% of the
            % length of the fuselage, NACA0012 airfoil for tail using a 30 degree
            % sweep, wetted surface is a complete guess. No miscellaneous drag and
            % leakage and protuberance is 10%.

            C_fi = [aerodynamics.C_f_f(M, V_inf, L_fus, nu_inf), aerodynamics.C_f_f(M, V_inf, 3.2, nu_inf), aerodynamics.C_f_f(M, V_inf, 0.2*L_fus, nu_inf), aerodynamics.C_f_f(M, V_inf, 0.3*L_fus)];
            FF_i = [aerodynamics.FF_fuselage_canopy(0.2*L_fus, pi*(0.15*L_fus)^2/4), aerodynamics.FF_store(3.2, pi*.35^2/4), aerodynamics.FF_fuselage_canopy(0.2*L_fus, (0.05*L_fus)^2/4), aerodynamics.FF_tail(0.12, 0.3, M, 30)];
            Q_i = [1, 1.3, 1, 1.04];
            S_wet_i = [0.4*S_wet, 0.05*S_wet, 0.08*S_wet, 0.09*S_wet];
            C_D_misc = 0;
            if M < 1
                for i = 1:length(C_fi)
                    tot = tot + C_fi(i)*FF_i(i)*Q_i(i)*S_wet(i);
                end
                C_D0 = 1.1*(tot/S_ref+C_D_misc);
            else
                for i = 1:length(C_fi)
                    tot = tot + C_fi(i)*S_wet_i(i);
                end
                C_D0 = 1.1*(tot/S_ref+C_D_misc);
            end
        end

        function C_Dw = C_Dw_f(... Wave Drag Coefficient
                M,... Mach number
                A_max,... Maximum cross-sectional area, m^2
                l,... Body length, m
                Lambda_LE,... Leading edge sweep, RAD
                Lambda_025c,... Quarter-chord sweep, DEG
                S_ref,... Reference area, m^2
                E_WD) % Wave-drag efficiency factor 1.2 < E_WD < 3

            M_DD = 0.15*Lambda_025c+0.74; % Crude divergent drag mach approximation from drag maps
            M_cr = M_DD-0.08;

            Dq_SH = 9*pi/2*(A_max/l)^2; % Sears-Haack correlation
            Dq_w = E_WD*(1-0.2*(M-1.2)^0.57)*(1-pi*Lambda_LE^0.77/100)*Dq_SH;
            B = E_WD*(1-pi*Lambda_LE^0.77/100)*Dq_SH;

            if M < M_cr % Drag rise curves
                C_Dw = 0;
            elseif M < M_DD
                C_Dw = 0.002/0.08*M-0.002/0.08*M_cr;
            elseif M < 1.05
                C_Dw = (B-0.002)/(1.05-M_DD)*M+(0.002/0.08*M_DD-0.002/0.08*M_cr-(B-(0.002/0.08*M_DD-0.002/0.08*M_cr)/(1.05-M_DD)*M_DD));
            elseif M < 1.2
                C_Dw = B;
            else
                C_Dw = Dq_w/S_ref;
            end
        end
        
        function C_L = C_L_f(L, q_inf, S_ref) % Lift coefficient
            C_L = L/(q_inf*S_ref);
        end

        function alpha_i = alpha_i_f(C_L, A, delta) % Induced AoA, rad
            alpha_i = C_L/(pi*A)*(1+delta);
        end
        
        function sweepx = Sweeptf(Lambda_x, x, y, A, lambda)% Sweep angle (deg) of y% chord from sweep (rad) at x% chord, AR, and taper
            sweepx = atand(tan(Lambda_x)-4*(y-x)*(1-lambda)/(A*100*(1+lambda)));
        end

        function a = C_L_alpha_sub_f(... Subsonic lift correction, rad^1
                A,... Aspect ratio
                M,... Mach number
                a_0,... slope of C_l vs alpha curve for airfoil, rad^-1
                d,... diameter of fuselage, m
                b,... Wingspan, m
                S_exposed,... Exposed area, m^2
                S_ref,... Reference area, m^2
                Lambda_maxt) % Sweep at max thickness, deg

            % Simplified Expression: C_L_alpha = @(C_l_alpha, Lambda_05c, A, e) C_l_alpha*cosd(Lambda_05c)/... Slope of C_L vs alpha curve from slope of airfoil c_l vs alpha curve, half-chord sweep, aspect ratio, and oswald efficiency
            %     (sqrt(1+(C_l_alpha*cosd(Lambda_05c)/(pi*A*e))^2)+(C_l_alpha*cosd(Lambda_05c))/(pi*A*e));

            beta = sqrt(1-M^2);
            eta = a_0*beta/(2*pi);
            F = 1.07*(1+d/b)^2;
            a = (2*pi*A)/( 2+sqrt(4 + A^2*beta^2/eta^2)*( 1 + tand(Lambda_maxt)^2/beta^2 ) )*S_exposed/S_ref*F;

        end

        function a = C_L_alpha_sup_f(... Supersonic lift correction, rad^-1
                A,... Aspect ratio
                b,... Wingspan, m
                lambda,... Taper ratio
                Lambda_LE,... Leading edge sweep, deg
                S_ref,... Reference area, m^2
                S_exposed,... Exposed surface area, m^2
                M,... Mach number
                d... Fuselage diameter, m
                )

            beta = sqrt(M^2-1);
            factor = beta/tand(Lambda_LE);
            vv = A*tand(Lambda_LE);
            F = 1.07*(1+d/b)^2;

            if M > 1/cos(Lambda_LE)
                disp('Supersonic lift correction outside of bounds')
                disp('Mach:')
                M %#ok<NOPRT>
                disp('Lambda_LE:')
                Lambda_LE %#ok<NOPRT>
                a = NaN;
                return
            end

            % it is so unbelievably over
            % Indices; (Atan(Lambda_LE), factor, lambda)

            LHS(:,:,1) = [
                0.4, 0.4, 0.4, 0.4, 0.4, 0.5;
                0.8, 0.8, 0.8, 0.8, 0.85, 0.9;
                1.5, 1.52, 1.56, 1.62, 1.7, 1.8;
                3.1, 3.1, 3.15, 3.2, 3.3, 3.4;
                4.8, 4.7, 5.1, 4.6, 4.3, 3.7;
                6.3, 6, 5.4, 4.9, 4.4, 4;
                6.3, 6.4, 5.7, 5.2, 4.7, 4.3;
                6.3, 6.4, 6.2, 5.6, 5, 4.5
                ]; % lambda = 0

            LHS(:,:,2) = [
                0.4, 0.4, 0.4, 0.4, 0.4, 0.4;
                0.7, 0.7, 0.7, 0.7, 0.8, 0.9;
                1.5, 1.55, 1.6, 1.7, 1.8, 1.9;
                3.1, 3.2, 3.6, 3.6, 3.5, 3.3;
                4.2, 4.7, 4.6, 4.3, 3.8, 3.6;
                5.6, 5.7, 5.5, 4.9, 4.5, 4.1;
                5.8, 5.8, 5.8, 5.4, 4.8, 4.4;
                5.9, 5.9, 6, 5.7, 5.3, 4.8
                ]; % lambda = 1/5

            LHS(:,:,3) = [
                0.4, 0.4, 0.4, 0.4, 0.4, 0.45;
                0.8, 0.75, 0.7, 0.7, 0.75, 0.9;
                1.6, 1.6, 1.6, 1.7, 1.75, 1.8;
                3.2, 3.4, 3.7, 3.6, 3.4, 3.2;
                4.7, 5.1, 4.7, 4.4, 4, 3.7;
                5.4, 5.5, 5.4, 4.9, 4.5, 4.1;
                5.6, 5.6, 5.5, 5.3, 4.8, 4.4;
                5.7, 5.7, 5.7, 5.6, 5.2, 4.6
                ]; % lambda = 1/4

            LHS(:,:,4) = [
                0.4, 0.4, 0.4, 0.4, 0.4, 0.4;
                0.8, 0.8, 0.8, 0.8, 0.85, 0.9;
                1.6, 1.6, 1.65, 1.7, 1.8, 1.9;
                3.2, 3.4, 3.8, 3.7, 3.5, 3.2;
                4.7, 4.8, 4.6, 4.4, 4, 3.7;
                5.2, 5.2, 5.1, 4.9, 4.4, 4.1;
                5.5, 5.5, 5.4, 5.3, 4.8, 4.4;
                5.6, 5.6, 5.7, 5.6, 5.2, 4.7;
                ]; % lambda = 1/3

            LHS(:,:,5) = [
                0.4, 0.4, 0.4, 0.4, 0.4, 0.4;
                0.8, 0.8, 0.8, 0.8, 0.9, 0.9;
                1.5, 1.6, 1.7, 1.8, 1.85, 1.9;
                3.1, 3.4, 3.7, 3.5, 3.3, 3.1;
                4.4, 4.4, 4.4, 4.3, 3.9, 3.6;
                4.9, 4.8, 4.7, 4.6, 4.3, 4;
                5.1, 5.1, 5.1, 5, 4.7, 4.4;
                5.2, 5.2, 5.3, 5.3, 5.2, 4.6
                ]; % lambda = 1/2

            LHS(:,:,6) = [
                0.4, 0.4, 0.45, 0.5, 0.55, 0.6;
                0.7, 0.7, 0.75, 0.9, 1, 1.2;
                1.5, 1.5, 1.6, 1.7, 1.8, 2;
                3.2, 3.1, 3, 3, 2.9, 2.8;
                3.8, 2.9, 3.8, 3.5, 3.4, 3.3;
                4.1, 4.2, 4, 3.9, 3.8, 3.7;
                4.4, 4.5, 4.3, 4.2, 4.2, 4.3;
                4.6, 4.7, 4.5, 4.4, 4.6, 4.7
                ]; % lambda = 1

            RHS(:,:,1) = [
                0.5, 0.5, 0.6, 1, 1.8, 4;
                0.9, 1.1, 1.8, 1.9, 3.2, 4;
                1.8, 2.2, 2.6, 3.3, 3.8, 4;
                3.4, 3.6, 3.7, 3.8, 3.9, 4;
                3.7, 3.8, 3.8, 3.9, 3.9, 4;
                4, 4, 4, 4, 4, 4;
                4.3, 4.25, 4.2, 4.15, 4.1, 4;
                4.5, 4.4, 4.3, 4.2, 4.1, 4;
                ]; % lambda = 0

            RHS(:,:,2) = [
                0.4, 0.5, 0.7, 1.1, 2.2, 4;
                0.9, 1.1, 1.8, 2.4, 3.5, 4;
                1.9, 2.3, 2.9, 3.5, 3.8, 4;
                3.3, 3.5, 3.7, 3.8, 3.9, 4;
                3.6, 3.7, 3.8, 3.9, 3.95, 4;
                4.1, 4.2, 4.2, 4.1, 4, 4;
                4.4, 4.4, 4.3, 4.2, 4.1, 4;
                4.8, 4.6, 4.4, 4.3, 4.1, 4
                ]; % lambda = 1/5

            RHS(:,:,3) = [
                0.45, 0.5, 0.7, 1.2, 2.3, 4;
                0.9, 1.1, 2, 2.3, 3.5, 4;
                1.8, 2.3, 3, 3.5, 3.7, 4;
                3.2, 3.5, 3.6, 3.8, 3.9, 4;
                3.7, 3.8, 3.9, 4, 4, 4;
                4.1, 4.2, 4.1, 4, 4, 4;
                4.4, 4.4, 4.3, 4.2, 4.1, 4;
                4.8, 4.7, 4.4, 4.3, 4.1, 4
                ]; % lambda = 1/4

            RHS(:,:,4) = [
                0.4, 0.5, 0.7, 1.1, 2.8, 4;
                0.9, 1.2, 1.9, 2.7, 3.5, 4;
                1.9, 2.5, 3, 3.5, 3.7, 4;
                3.2, 3.5, 3.7, 3.8, 3.9, 4;
                3.7, 3.8, 3.9, 4, 4, 4;
                4.1, 4.2, 4.1, 4, 4, 4;
                4.4, 4.5, 4.3, 4.2, 4.1, 4;
                4.7, 4.6, 4.4, 4.2, 4.1, 4
                ]; % lambda = 1/3

            RHS(:,:,5) = [
                0.4, 0.5, 0.7, 1.2, 2.8, 4;
                0.9, 1.2, 1.7, 1.9, 3.4, 4;
                1.9, 2.4, 2.8, 3.4, 3.8, 4;
                3.1, 3.4, 3.6, 3.8, 3.9, 4;
                3.6, 3.8, 3.9, 4, 4, 4;
                4, 4.2, 4.2, 4.1, 4, 4;
                4.4, 4.5, 4.3, 4.2, 4.1, 4;
                4.6, 4.6, 4.4, 4.3, 4.1, 4
                ]; % lambda = 1/2

            RHS(:,:,6) = [
                0.6, 0.7, 0.9, 1.7, 2.5, 4;
                1.2, 1.4, 1.6, 2.3, 3.3, 4;
                2, 2.3, 2.7, 3.2, 3.7, 4;
                2.8, 3.2, 3.5, 3.7, 3.8, 4;
                3.3, 3.7, 3.9, 4, 4, 4;
                3.7, 4.2, 4.2, 4.1, 4, 4;
                4.3, 4.5, 4.3, 4.2, 4.1, 4;
                4.7, 4.8, 4.5, 4.3, 4.1, 4
                ]; % lambda = 1

            AtanLLEx = [0.25, 0.5, 2, 3, 4, 5, 6];
            lambdax = [0, 0.2, 0.25, 1/3, 0.5, 1];

            for i = 1:length(AtanLLEx)
                if vv > AtanLLEx(i) && vv < AtanLLEx(i+1)
                    break
                end
            end
            for j = 1:length(lambdax)
                if lambda > lambdax(j) && lambda < lambdax(j+1)
                    break
                end
            end

            if factor > 1
                factor = 1/factor;
                factorx = 1:0.2:0;

                lower = interp1(factorx, RHS(i, :, j), factor);
                upper = interp1(factorx, RHS(i+1, :, j), factor);
                ml = interp1([AtanLLEx(i),AtanLLEx(i+1)], [lower, upper], vv);

                lower = interp1(factorx, RHS(i, :, j+1), factor);
                upper = interp1(factorx, RHS(i+1, :, j+1), factor);
                mu = interp1([AtanLLEx(i),AtanLLEx(i+1)], [lower, upper], vv);

                betaCL = interp1([lambdax(j),lambdax(j+1)], [ml, mu], lambda);
                C_L_alpha = betaCL/beta;

            else
                factorx = 0:0.2:1;

                lower = interp1(factorx, LHS(i, :, j), factor);
                upper = interp1(factorx, LHS(i+1, :, j), factor);
                ml = interp1([AtanLLEx(i),AtanLLEx(i+1)], [lower, upper], vv);

                lower = interp1(factorx, LHS(i, :, j+1), factor);
                upper = interp1(factorx, LHS(i+1, :, j+1), factor);
                mu = interp1([AtanLLEx(i),AtanLLEx(i+1)], [lower, upper], vv);

                tanLCL = interp1([lambdax(j),lambdax(j+1)], [ml, mu], lambda);
                C_L_alpha = tanLCL/tand(Lambda_LE);

            end

            a = C_L_alpha*S_exposed/S_ref*F;

        end
        
        function [... Maximum lift calculations
                C_L_max,... Max lift coefficient for wing
                alpha_C_L_max... AoA of max lift, rad
                ] = maxlift(...
                deltay,... t(0.06*C)-t(0.0015*C) from airfoil
                alpha_ZL,... Zero lift AoA, rad
                C_l_max,... Max lift coefficient for airfoil
                a,... C_L_alpha (slope of CL vs alpha curve for wing), rad^-1
                Lambda_LE,... Leading edge sweep, rad
                A,... Aspect ratio
                lambda... Taper ratio
                )
            
                vals = [
                    0.9, 0.92, 0.97, 1.03, 1.1, 1.19, 1.25;
                    0.9, 0.91, 0.95, 0.98, 1.04, 1.12, 1.18;
                    0.9, 0.91, 0.93, 0.95, 0.98, 1, 1.02;
                    0.9, 0.89, 0.886, 0.883, 0.88, 0.87, 0.86;
                    0.9, 0.88, 0.87, 0.85, 0.81, 0.78, 0.73;
                    0.9, 0.87, 0.86, 0.81, 0.77, 0.69, 0.59
                ];
            
               
                dyx = 1.4:0.2:2.4;
                LLEx = 0:10:60;
               
                
                for i = 1:size(vals, 1)
                    if deltay < dyx(i+1) && deltay > dyx(i)
                        break
                    end
                end
                
                l = interp1(LLEx, vals(i,:), rad2deg(Lambda_LE));
                u = interp1(LLEx, vals(i+1,:), rad2deg(Lambda_LE));
                CLCl = interp1(dyx(i:i+1), [l, u], deltay);
                
                S_flapped_TE = 0.15*S_ref; % Approximation of slatted surface area
                S_flapped_LE = 0.08*S_ref; % Approximation of flapped surface area
                Lambda_HL_TE = aerodynamics.Sweeptf(Lambda_LE, 0, 5, A, lambda); % Sweep of trailing edge hinge line, deg
                Lambda_HL_LE = aerodynamics.Sweeptf(Lambda_LE, 0, 85, A, lambda); % Sweep of leading edge hinge line, deg
                
                DC_Lmax_TE = 0.9*aerodynamics.DC_l_max_TE*S_flapped_TE/S_ref*cosd(Lambda_HL_TE);
                DC_Lmax_LE = 0.9*aerodynamics.DC_l_max_LE*S_flapped_LE/S_ref*cosd(Lambda_HL_LE);
                DC_Lmax = DC_Lmax_TE + DC_Lmax_LE;
                
                vals = [
                    1.9, 2.1, 3.3, 4.5, 7.4, 10, 14;
                    0.2, 1, 2.4, 3.9, 5.7, 7.5, 9.8;
                    1.2, 1.7, 2.4, 3.2, 4.2, 5.3, 6.6;
                    2.3, 2, 2.1, 2.3, 2.5, 2.9, 3.2
                ];
                
                dyx = [1.2, 2, 3, 4];
                LLEx = 0:10:60;
                
                for i = 1:size(vals, 1)
                    if deltay < dyx(i+1) && deltay > dyx(i)
                        break
                    end
                end
                
                l = interp1(LLEx, vals(i,:), rad2deg(Lambda_LE));
                u = interp1(LLEx, vals(i+1,:), rad2deg(Lambda_LE));
                Dalpha_C_L_max = interp1(dyx(i:i+1), [l, u], deltay);
                
                C_L_max = C_l_max*CLCl + DC_Lmax;
                alpha_C_L_max = C_L_max/a + alpha_ZL + Dalpha_C_L_max;
        end

        function [... Wing Sizing Function
                A,... Aspect ratio
                MAC,... Mean aerodynamic chord, m
                Lambda_TE,... Trailing-edge sweep, deg
                Lambda_025c,... Quarter-chord sweep, deg
                S_canard,... Canard reference area
                Ld,... Lift to drag ratio for given conditions
                C_Lmax,... Maximum lift coefficient
                alpha_stall... Stall angle of attack, 1/rad
                ] = wing(...
                W0,... MTOW, kg
                V_inf,... Freestream speed, m/s
                nu_inf,... Kinematic viscosity, m^2/s
                q_inf,... Dynamic pressure, Pa
                S_ref,... Reference area, m^2
                b,... Wingspan, m
                lambda,... Taper ratio
                Lambda_LE,... Leading edge sweep, deg
                M,... Mach number
                tc_root,... Thickness ratio of airfoil at root
                tc_tip,... Thickness ratio of airfoil at tip
                A_top,... Fuselage top area, m^2
                A_side,... Fuselage side area, m^2
                xc_m,... Normalized location of maximum thickness of airfoil
                deltay_root,... t(x=0.06*C)-t(x=0.0015*C) from airfoil
                deltay_tip,... t(x=0.06*C)-t(x=0.0015*C) from airfoil
                d... Diameter of fuselage, m
                )

            % Should probably also check pitching moment to size canards to resize wing
            % Should probably also account for wing twist

            Lambda_maxt = Sweeptf(deg2rad(Lambda_LE), 0, xc_m, A, lambda); % Sweep angle of maximum thickness, deg
            Lambda_025c = Sweeptf(deg2rad(Lambda_LE), 0, 25, A, lambda); % Sweep angle of quarter-chord, deg
            Lambda_TE = Sweeptf(deg2rad(Lambda_LE), 0, 100, A, lambda); % Sweep angle at trailing edge, deg

            C_L = C_L_f(W0, q_inf, S_ref);

            S_canard = 0.07*S_ref; % guess at canard size

            A = b^2/(S_ref+S_canard); % Aspect ratio (not using taper formula because it isnt correct apparently)
            MAC = A/b;
            e = ef(A, Lambda_LE);
            delta = deltaf(e);

            C_Di = C_Di_f(C_L, A, delta); % Induced drag coefficient

            S_exposed = S_ref-C_root*0.1*S_ref; % Estimate of exposed area, m^2
            tc = (tc_root+tc_tip)/2;
            if tc < 0.05
                S_wet = 2.003*S_exposed;
            else
                S_wet = (1.977+0.52*tc)*S_exposed;
            end
            S_wet_fuselage = 1.7*(A_top+A_side);
            S_wet_tot = S_wet + S_wet_fuselage;
            
            Re = Ref(V_inf, MAC, nu_inf);
            [a0_root, a0_tip] = airfoil.a0(Re);
            
            if M < 1
                C_L_alpha = C_L_alpha_sub_f(A, M, (a0_root+a0_tip)/2, d, b, S_exposed, S_ref, Lambda_maxt);
            else
                C_L_alpha = C_L_alpha_sup_f(A, b, lambda, Lambda_LE, S_ref, S_exposed, M, d);
            end

            [alpha_ZL_root, alpha_ZL_tip] = airfoil.zerolift(Re);

            alpha_root = (C_L+alpha_ZL_root*C_L_alpha)/C_L_alpha; % true angle of attack, rad
            alpha_tip = (C_L+alpha_ZL_tip*C_L_alpha)/C_L_alpha;

            alpha_eff_root = alpha_root - alpha_i_f(C_L, A, delta); % effective angle of attack, rad
            alpha_eff_tip = alpha_tip - alpha_i_f(C_L, A, delta); % effective angle of attack, rad


            % Should probably account for trim drag C_Di_trim
            % Need to find C_Lmax which is a function of C_lmax
            % Should also get AoA_stall from this to use in trade studies
            % Should probably also account for twist somehow

            [C_lmax_root, C_lmax_tip] = airfoil.clmax(Re);

            [C_Lmax_root, alpha_stall_root] = maxlift(deltay_root, alpha_ZL_root, C_lmax_root, C_L_alpha, Lambda_LE, A, lambda);
            [C_Lmax_tip, alpha_stall_tip] = maxlift(deltay_tip, alpha_ZL_tip, C_lmax_tip, C_L_alpha, Lambda_LE, A, lambda);
            alpha_stall = min([alpha_stall_root, alpha_stall_tip]);
            C_Lmax = mean([C_Lmax_root, C_Lmax_tip]);
    
            [Cd_root, Cd_tip, Cl_root, Cl_tip, Cm_root, Cm_tip] = airfoil.airfoil(alpha_root, alpha_tip, alpha_eff_root, alpha_eff_tip, Re);
            
            % C_D0 = 0.004*S_wet_tot/S_ref;
            C_D0 = (Cd_root + Cd_tip)/2 + C_D0_f(M, V_inf, L_fus, nu_inf, S_wet_tot, S_ref); % Parasite drag coefficient

            C_Dw = C_Dw_f(M, A_max, l, deg2rad(Lambda_LE), Lambda_025c, S_ref, E_WD); % Wave drag coefficient

            C_D = C_D0 + C_Di + C_Dw; % Total drag coefficient

            Ld = C_L/C_D; % Lift-to-drag ratio

            % Abandoned canard sizing
            % % S_canard = ?
            % 
            % C_m = ; % Quarter-chord pitching moment of wing from airfoil
            % M_w = b*q_inf*MAC^2*C_m;
            % C_m_c = ;
            % M_c = b*q_inf*MAC_c^2*C_m_c;

        end
    end
end