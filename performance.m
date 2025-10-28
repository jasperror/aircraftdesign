% ENGR 493 Intro to Aircraft Design
% Project Phase 2: Wing Design - Performance (wing sizing)
% Jasper Palmer
classdef performance
    properties (Constant)
        V_app = 67;
        V_stall = performance.V_app/1.2;
        g = 9.81;
    end
    methods (Static)
        function S_ref = wingloading(W0, rho_landing, C_Lmax)
            WS_stall = 0.5*rho_landing*performance.V_stall^2*C_Lmax/performance.g;
            % Lecture note formula is wrong ^ should divide by g to get
            % correct units
            % WS_carriertakeoffland = 0.5*rho_landing*(V_end+V_WoD+DeltaV_thrust)^2*C_Lmax/(V_TO/V_stall)^2;
            % minWS = min(WS_stall, WS_carriertakeoffland);
            % S_ref = W0/minWS;
            S_ref = W0/WS_stall;
        end

        function V_end = V_end(W0)
            if W0 < 19000
                V_end = 77;
            elseif W0 < 45000
                V_end = -13/26000*W0+86;
            else
                V_end = 64;
            end
        end
    end
end