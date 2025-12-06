% ENGR 493 Intro to Aircraft Design
% Performance -- Wing loading (W/S) and thrust-to-weight ratio (T/W)
% Jasper Palmer

% TODO -- Should make a function that converts thrust to weight to/from design condition, this allows for substituting the performance plot of a given engine in the future
% also need to get drag for takeoff, landing, and cruise phases?

classdef performance
    properties (Constant)
        l_to = 100; % Carrier takeoff distance available, m
    end
    methods (Static)

        function V_end = V_end(W0)
            W0 = W0*2.2/1000; % kg to klbs
            V_endf = @(W) interp1(10:10:100,[173, 160, 152, 145, 140, 135, 131, 126, 122, 118],W); % % Approximation from RFP at CSV=200
            if W0 < 10
                V_end = 173; % assuming the maximum speed of the catapult is 173kts
            elseif W0 <= 100
                V_end = V_endf(W0);
            else
                V_end = V_endf(100);
                % warning("W0 above 100klbs")
            end
            V_end = V_end*0.514; % kts to m/s
        end

        function WS_0 = stall(h, C_Lmax, W_0x, V_s)
            WS_stall = 0.5*atmosphere.rho(h)*V_s^2*C_Lmax;
            WS_0 = WS_stall/W_0x;
        end

        function WS_0 = takeoff(... Wing loading for takeoff
                type,... 0 for jet (FAR takeoff), 1 for jet, 2 for propeller
                n,... Number of engines
                gr,... True/false, true = ground roll only, not 50ft obstacle (not for FAR takeoff)
                dist,... Available takeoff distance, m
                h,... Takeoff height, m
                C_Lmax,... Max lift coefficient
                W_0to,... Takeoff weight fraction
                TW_0... Thrust to weight ratio
                )
            dist = 3.28084*dist/1000; % converting takeoff distance from m to kft
            switch type
                case 0
                    TOP = (1.7673*n+20.2851)*dist+4*n-8; % first order polynomial regression fit of TOP plot for FAR data as a function of n
                    % Probably shouldnt apply to n>4 || n<1 but left them in anyways
                    % https://math.stackexchange.com/questions/680646/get-polynomial-function-from-3-points
                case 1
                    if gr
                        TOP = 410/9*dist+40;
                    else
                        TOP = 415/10.1*dist+28;
                    end
                case 2
                    if gr
                        TOP = 500/4.4*dist+30;
                    else
                        TOP = 485/5.2*dist+11;
                    end
                otherwise
                    warning("Takeoff type invalid")
            end

            TOP = 47.8803*TOP; % converting lb/ft^2 to Pa

            C_LTO = C_Lmax/(1.1)^2; % Assuming VTO=1.1Vs
            TW_to = TW_0*atmosphere.sig(h)/W_0to;

            WS_to = TOP*atmosphere.sig(h)*C_LTO*TW_to;
            WS_0 = WS_to/W_0to;
        end

        function WS_0 = Vmax(... Wing loading for maximum speed
                TW_0,... Design thrust to weight ratio
                W_0x,... Weight fraction
                h,... Altitude, m
                V_max,... Maximum speed, m/s
                C_D0,... Parasitic drag coefficient
                K... Performance coefficient
                )

            TW = TW_0*atmosphere.sig(h)/W_0x;
            q = atmosphere.q(V_max, atmosphere.rho(h));
            x = TW^2 - 4*K*C_D0;
            if x < 0
                WS = NaN;
            else
                WS = ( TW-sqrt(x) )/(2*K/q);
            end
            WS_0 = WS/W_0x;
        end

        function WS_0 = cl(... Climb wing loading
                type,... 0 for jet (at Vy), 1 for prop (at Vy)
                V_v,... Vertical velocity, ft/min
                C_D0,... Parasite drag coefficient
                TW_0,... Thrust to weight ratio
                W_0cl,... Weight fraction to climb
                h,... Height at climb
                K... Performance coefficient
                )
            V_v = V_v/196.9; % ft/min to m/s
            rho = atmosphere.rho(h);
            sig = atmosphere.sig(h);
            TW_cl = TW_0*sig/W_0cl;

            switch type
                case 0
                    A = TW_cl-C_D0*1.25^2*sqrt(K/(C_D0))-K/1.25^2/sqrt(K/(C_D0));
                    if A <= 0
                        WS_0 = NaN;
                        return
                    else
                        B = V_v/(1.25*sqrt(2/rho*sqrt(K/(C_D0))));
                    end
                case 1
                    A = TW_cl-C_D0*sqrt(K/(3*C_D0))-K/sqrt(K/(3*C_D0));
                    if A <= 0
                        WS_0 = NaN;
                        return
                    else
                        B = V_v/(sqrt(2/rho*sqrt(K/(3*C_D0))));
                    end
                otherwise
                    warning("Climb type invalid")
                    WS_0 = NaN;
                    return
            end

            WS_cl = B^2/A^2; % Derived from lecture notes and textbook
            WS_0 = WS_cl/W_0cl;

        end

        function WS_0 = ceil(... Ceiling wing loading
                type,... 0 for jet, 1 for prop
                ceiling,... Ceiling type, 0 for absolute, 1 for service, 2 for cruise
                C_D0,... Parasite drag coefficient
                TW_0,... Thrust to weight ratio
                W_0cr,... Weight fraction to ceiling
                h,... Height at ceiling
                K... Performance coefficient
                )
            switch ceiling
                case 0 % Absolute ceiling
                    V_v = 0; % ft/min
                case 1 % Service Ceiling
                    V_v = 100;
                case 2 % Cruise Ceiling
                    V_v = 300;
                otherwise
                    warning("Ceiling wing loading type invalid")
                    WS_0 = NaN;
                    return
            end

            WS_0 = performance.cl(type, V_v,C_D0,TW_0,W_0cr,h,K);
        end

        function WS_0 = cbto(... Carrier based takeoff, from §3.4
                V_WoD,... Wind speed over deck, m/s
                V_toVs,... Takeoff speed divided by stall speed
                rho,... Air density at takeoff, Pa
                W_0,... MTOW, kg
                W_0x,... Weight fraction
                TW_to,... Thrust to weight ratio at takeoff
                C_Lmax... Maximum lift coefficient
                )
            DV_T = sqrt(2*TW_to*performance.l_to); % From kinematics
            WS_to = 0.5*rho*(performance.V_end(W_0) + V_WoD + DV_T)^2*C_Lmax/V_toVs^2;
            WS_0 = WS_to/W_0x;
        end

        function WS_0 = cbl(... Carrier based landing
                W_l,.... Landing weight, kg
                rho,... Air density at landing, Pa
                C_Lmax,... Max lift coefficient
                W_0x... Landing weight fraction
                )
            W_l = W_l*2.2/1000; % kg to klbs
            if W_l <= 40
                V_app = 145;
            elseif W_l <= 70
                %V_app = -0.9*W_l+168; Relation from textbook for 40-50klbs
                V_app = -1.4*W_l+201; % Relation from RFP
            else
                V_app = NaN;
                warning("Over MLW")
            end
            h = atmosphere.h(rho);
            V_app = V_app*0.514; % kts to m/s
            V_s = V_app/1.2; % approximating stall speed from approach speed
            WS_0 = performance.stall(h, C_Lmax, W_0x, V_s);
        end

        function WS_0 = sturn(... Sustained turn wing loading
                TW_0,... Design thrust to weight ratio
                W_0x,... Weight fraction at analysis point
                ...% V_max,... Maximum speed, m/s
                K,... Performance coefficient
                n,... Design maximum load factor
                psidot,... Desired rate of turn, deg/s
                C_D0,... Parasitic drag coefficient
                h... Altitude, m
                )
            K = K; % Accounting for reduced wing efficiency in steep bank
            V = performance.V_maxturn(psidot,n);
            % if V >= V_max % Check if best turn speed is within max speed
            %     V = V_max;
            %     n = sqrt((psidot*pi/180*V/const.g)^2+1); % if not, fix it
            % end

            q = atmosphere.q(V,atmosphere.rho(h));

            TW_x = TW_0*atmosphere.sig(h)/W_0x;

            if TW_x < 2*n*sqrt(C_D0*K)
                WS_0 = NaN;
                return
            end

            WS_x = ( TW_x - sqrt(TW_x^2 - 4*n^2*K*C_D0) )/(2*K*n^2/q);
            WS_0 = WS_x/W_0x;
        end

        function WS_0 = iturn(... Instatnaneous turn wing loading
                C_Lmax,... Maximum lift coefficient at instantaneous turn
                W_0x,... Weight fraction at analysis point
                V_max,... Maximum speed, m/s
                n,... Design maximum load factor
                psidot,... Desired rate of turn, deg/s
                h... Altitude, m
                )

            V = performance.V_maxturn(psidot,n);
            if V >= V_max
                V = V_max;
                n = sqrt((psidot*pi/180*V/const.g)^2+1);
            end

            q = atmosphere.q(V,atmosphere.rho(h));
            WS_x = q*C_Lmax/n;
            WS_0 = WS_x/W_0x;

        end

        function V = V_maxturn(... Speed for maximum turn rate
                psidot,... Turn rate, deg/s
                n... Load factor
                )
            psidot = psidot*pi/180; % deg/s to rad/s
            V = const.g/psidot*sqrt(n^2-1);
        end

        function [S_ref, T_0, V_maxturn, V_cr, V_be, V_s0, WS_0 , TW_0] = fighter(...
                psidot,... Desired turn rate, deg/s
                n,... Structural load limit
                M_cr_dash,... Cruise dash mach#
                M_sl_dash,... Sea-level dash mach#
                V_crg,... % Cruise speed, m/s
                V_s,... % Maximum stall speed at landing, m/s
                K,... Performance coefficient, []
                m_0g,... Guess MTOW, kg
                C_Lmax_0,... Max C_L at SL
                W_0l,... Landing weight fraction
                C_D0_cr,... Parasite drag coefficient at cruise
                W_0cr,... Cruise weight fraction
                h_cr,... Cruise height, m
                W_0midmission,... Mid mission weight fraction
                C_D0_maxV_cr,... Parasite drag coefficient at max speed, cruise
                C_D0_maxV_SL,... Parasite drag coefficient at max speed, sea-level
                plotit... Plot it?
                )

            TW_0r = 0.1:0.001:1;
            WS_0t = zeros(7, length(TW_0r));
            V_maxturn = performance.V_maxturn(psidot, n);

            WS_0t(1, :) = performance.cbl(m_0g*W_0l, atmosphere.rho(0), C_Lmax_0, W_0l);
            WS_0t(2, :) = performance.stall(0, C_Lmax_0, W_0l, V_s);
            for i = 1:length(TW_0r)
                if i == 2
                    warning('off','all')
                end
                % Both
                WS_0t(3, i) = performance.cbto(0, 1.1, atmosphere.rho(0), m_0g, 1, TW_0r(i), C_Lmax_0);
                WS_0t(4, i) = performance.ceil(0, 2, C_D0_cr, TW_0r(i), W_0cr, h_cr, K);
                % Dogfighter
                WS_0t(5, i) = performance.sturn(TW_0r(i), W_0midmission, K, n, psidot, C_D0_cr, 20000*0.3048);
                if isnan(WS_0t(5,i)) % Checks if the minimum T/W is met for the desired turn
                    WS_0t(5,i) = WS_0t(1,i)+1;
                end
                WS_0t(6, i) = performance.Vmax(TW_0r(i), W_0midmission, 30000*0.3048, atmosphere.V(M_cr_dash, 30000*0.3048), C_D0_maxV_cr, K);
                % Strike
                WS_0t(7, i) = performance.Vmax(0.8*TW_0r(i), W_0midmission, 0, atmosphere.V(M_sl_dash, 0), C_D0_maxV_SL, K);
            end
            % if mean(abs(WS_0t(5,:))-WS_0t(1,1)-1<0.001)==1 || mean(WS_0t(5,:) == WS_0t(1,:))==1
            if isequal(WS_0t(5,:),WS_0t(1,1)+1) || isequal(WS_0t(5,:), WS_0t(1,:))
                WS_0t(5,:) = NaN;
            end
            % Find minimum T/W
            maxws = real(max(WS_0t(4:end,:),[],1));
            minws = real(min(WS_0t(1:3,:),[],1));
            % figure
            % plot(TW_0r, maxws, TW_0r, minws)
            % legend('max','min')
            % ylim([0,500])

            [~, ~, i] = zerocrossrate(maxws-minws);
            i_swap = find(i(2:end)==1, 1)+1;

            % i_swap = find(maxws-minws<0.05,1);

            minws = minws(i_swap+1:end);
            maxws = maxws(i_swap:end);
            good = [maxws, minws];
            TWmax = TW_0r(i_swap:end);
            TWmin = TW_0r(i_swap+1:end);
            TW = [TWmax, TWmin];
            itwmin = find(TW-min(TW)<0.002,1,'last'); % If T/W min is flat, find highest W/S
            tw = TW(itwmin);
            ws = good(itwmin);

            if plotit
                performance.plot(WS_0t, TW_0r, tw, ws) % Plot it
            end

            TW_0 = tw; % Thrust to weight, N/N
            WS_0 = ws; % Wing loading, N/m^2
            % W_0g
            % WS_0
            try
                S_ref = m_0g*const.g/WS_0; % Wing reference area, m^2
            catch exception
                S_ref = NaN;
                T_0 = S_ref;
                V_maxturn = S_ref;
                V_cr = S_ref;
                V_be = NaN;
                V_s0 = NaN;
                WS_0 = NaN;
                TW_0 = NaN;
                return
                % the entire module just needs to be redone to be less
                % terrible to remove this
            end
            
            T_0 = TW_0*(m_0g*const.g); % Design thrust T_0, N
            V_cr = sqrt(2/atmosphere.rho(h_cr)*WS_0*sqrt(K/C_D0_cr));
            % err = abs(V_cr-V_crg);
            % V_crg = V_cr;
            V_be = sqrt(2/atmosphere.rho(h_cr)*WS_0*sqrt(K/(3*C_D0_cr)));
            V_s0 = sqrt(2/(atmosphere.rho(0)*C_Lmax_0)*WS_0);
            % end
        end

        function plot(WS_0t, TW_0r, tw, ws)
            % close all
            figure
            hold on
            plot(WS_0t(1,:), TW_0r, WS_0t(2,:), TW_0r,WS_0t(3,:), TW_0r,WS_0t(4,:), TW_0r,WS_0t(5,:), TW_0r,WS_0t(6,:), TW_0r,WS_0t(7,:), TW_0r)
            plot(ws, tw, 'k*')
            hold off
            legend('CBL','Vs','CBTO','CEIL','STURN','VMAX CR','VMAX LO','MIN')
            xlabel('W/S, N/m^2')
            ylabel('T/W, N/N')
            xl = min(max(WS_0t(1:2,1)*1.1));
            if isnan(xl)
                xl = 4000;
            end
            xlim([0,xl])
        end
    end
end