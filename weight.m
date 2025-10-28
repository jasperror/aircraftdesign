% ENGR 493 Intro to Aircraft Design
% Project Phase 2: Wing Design - Weight sizing
% Jasper Palmer
classdef weight
    properties
        W_to = 0.98; % takeoff weight fraction
        W_de = 0.9925; % descent weight fraction
        W_la = 0.995; % landing weight fraction
    end
    methods (Static)
        function W_cl = W_cl(M1, M2) % Climb weight fraction, default M1=0.2
            if M1 ~= 0.1
                if M2 < 1
                    W_cl = (1.0065-0.0325*M2)/(1.0065-0.0325*M1);
                else
                    if M1 < 1
                        W_cl = (0.991-0.007*M2-0.01*M2^2)/(1.0065-0.0325*M1);
                    else
                        W_cl = (0.991-0.007*M2-0.01*M2^2)/(0.991-0.007*M1-0.01*M1^2);
                    end
                end
            else
                if M2 < 1
                    W_cl = (1.0065-0.0325*M2);
                else
                    W_cl = 0.991-0.007*M2-0.01*M2^2;
                end
            end
        end

        function W_cr = W_cr(... Cruise weight fraction
                R,... Range, m
                C,... Consumption, 1/s
                V,... Speed, m/s
                Ld) % Lift to drag ratio in cruise
            W_cr = exp((-R*C)/(V*Ld));
        end

        function W_loi = W_loi(...
                E,... Endurance, s
                C,... Consumption, 1/s
                Ld) % Lift to drag ratio in loiter
            W_loi = exp((-E*C)/(Ld));
        end

        function WeW0 = WeW0(... Empty weight fraction
                W0,... MTOW, kg
                Kt) % Technology factor, 1 for conventional designs
            WeW0 = 2.11*W0^(-0.13)*Kt;
        end

        function [... Weight Sizing function
                W_0,... Correct MTOW, kg
                W_f,... Total fuel required, kg
                W_fi]... Fuel used on each segment, kg
                = weights(W_0g,... Guess MTOW, kg
                W_p,... Payload wt, kg
                W_d,... Drop wt, kg
                W_c,... Crew wt, kg
                K_rf,... Reserve fuel ratio (1.xx)
                Ld_loi,... (L/D)_loiter
                C_loi,... Consumption in loiter, 1/s
                M_loi,... Combat loiter mach number
                E,... Combat loiter endurance, s
                Ld_cr,... (L/D)_cruise
                C_cr,... Consumption in cruise, 1/s
                M_cr,... Cruise mach number
                V_cr,... Cruise speed, m/s
                R,... Combat radius (cruise radius), m
                C_lol,... Consumption in loiter-to-landing, 1/s
                Ld_lol,... (L/D)_{Loiter to landing}
                t_lol) % Loiter to landing time (holding time), s

            it = 0;
            MP = [weight.W_to, weight.W_cl(0.1, M_cr), weight.W_cr(R, C_cr, V_cr, Ld_cr), weight.W_de,...
                weight.W_loi(E, C_loi, Ld_loi), weight.W_cl(M_loi, M_cr),...
                weight.W_cr(R, C_cr, V_cr, Ld_cr), weight.W_de, weight.W_loi(t_lol, C_lol, Ld_lol), W_la];
            while (W_0 - W_0g) > 1 && it < 10000
                it = it+1;
                W_0g = W_0;
                W_fi = zeros(length(MP));
                W = zeros(length(MP)+1);
                W(1) = W_0;

                for i = 1:numel(MP)
                    W_fi(i) = W(i)*(1-MP(i));
                    W(i+1) = W(i)-W_fi(i);
                end

                W_f = K_rf*sum(W_fi);
                W_0 = (W_c+W_p+W_d+W_f)/(1-WeW0(W_0g, 1));
            end
            % W_fstored = zeros(length(MP));
            % for i = 1:length(W_fi)
            %     W_fstored(i) = W_f-sum(W_fi(1:i));
            % end
        end
    end
end