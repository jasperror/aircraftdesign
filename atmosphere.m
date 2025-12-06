classdef atmosphere

    properties (Constant)
        Atmosphere = table2array(readtable("Atmosphere.csv"));
        rho_0 = 1.225; % Sea-level air density, Pa
    end
    methods (Static)
        function V = V(M, A)
            V = M*interp1(atmosphere.Atmosphere(:,1),atmosphere.Atmosphere(:,4), A);
        end
        function nu = nu(A)
            nu = interp1(atmosphere.Atmosphere(:,1),atmosphere.Atmosphere(:,3),A);
        end
        function q = q(V, rho)
            q = 0.5*rho*V^2;
        end
        function M = M(V, A)
            M = V/interp1(atmosphere.Atmosphere(:,1),atmosphere.Atmosphere(:,4),A);
        end
        function rho = rho(A)
            rho = atmosphere.rho_0*interp1(atmosphere.Atmosphere(:,1),atmosphere.Atmosphere(:,2),A);
        end
        function sig = sig(A)
            sig = interp1(atmosphere.Atmosphere(:,1),atmosphere.Atmosphere(:,2),A);
        end
        function h = h(rho)
            h = interp1(atmosphere.Atmosphere(:,2),atmosphere.Atmosphere(:,1),rho/atmosphere.rho_0);
        end
    end
end