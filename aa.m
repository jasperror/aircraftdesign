% ENGR 493 Intro to Aircraft Design
% Project Phase 2: Wing Design - Aerodynamics of the wing
% Jasper Palmer

% todo -- wing twist?
classdef aa
    properties (Constant, Hidden) % Constants/assumptions
        DC_l_max_TE = 1.3; % Assuming slotted flap
        DC_l_max_LE = 0.3; % Assuming leading-edge flap
        E_WD = 2; % Wave-drag efficiency factor 1.2 < E_WD < 3
        k = 1.015e-5; % Surface roughness, m
        xcm = 0.385; % Chord-normalized location of maximum thickmess, avg between root and tip airfoil
        xcm_tail = 0.3; % Chord-normalized location of maximum thickness in tail airfoil, []
        Lambda_maxt_tail = 30; % Sweep angle of maximum thickness of tail, deg
        tc_tail = 0.12; % Chord-normalized thickness of tail airfoil, []
        L_stores = 3.2; % Length of missiles, m
        tc = (0.139 + 0.06)/2; % Average chord-normalized thickness in wing
        deltay_root = .04363/.00706; % t(x=0.06*C)-t(x=0.0015*C) from airfoil
        deltay_tip = .01888/.00944; % t(x=0.06*C)-t(x=0.0015*C) from airfoil
    end
    properties (Access=private, Hidden) % Wing internal properties
        Lambda_maxt % Sweep at maximum thickness, deg
        C_D_misc % Miscellaneous drag coefficient, []
        C_Dw % Wave drag coefficient, []

        M_DD % Drag-divergent mach, []

        C_l_max_root % Root airfoil maximum lift coefficient
        C_l_max_tip % Tip airfoil maximum lift coefficient
        A_LHS = struct2array(load("a.mat","LHS"));
        A_RHS = struct2array(load("a.mat","RHS"));
    end
    properties (Dependent) % Outputs
        Re % Reynolds number, []
        nu % Kinematic viscosity, m^2/s
        q % Dynamic pressure, Pa
        e % Oswald's efficiency, []
        K % Performance coefficient, []
        M % Mach number, []
        M_cr % Cruise mach, []

        Lambda_025c % Quarter-chord sweep, deg
        Lambda_05c % Mid-chord sweep, deg
        Lambda_TE % Trailing edge sweep, deg

        MAC % Mean aerodynamic chord, m
        A % Aspect ratio, []

        C_Di % Induced drag coefficient, []
        C_D % Total drag coeffiicent, []
        C_L % Lift coefficient, []
        C_L_max % Maximum lift coefficient, []
        C_D0 % Parasitic drag coefficient, []
        LD % Lift-to-drag ratio, []

        S_wet_fus % Wetted fuselage area, m^2
        S_wet % Wetted wing area, m^2
        S_exposed % Exposed surface area, m^2

        C_L_alpha % Wing lift correlation, 1/rad
        alpha % Angle of attack, rad
        alpha_root % Root angle of attack, rad
        alpha_tip % Tip angle of attack, rad
        alpha_eff_root % Root effective angle of attack, rad
        alpha_eff_tip % Tip effective angle of attack, rad
        alpha_i % Induced AoA, rad
        alpha_C_L_max % Angle of attack of maximum lift, rad
        
        C_root
    end

    properties (Constant, Hidden) % Airfoil Data
        a64206_50k = struct2array(load("airfoils.mat","a64206_50k")); % Data collected from airfoiltools.com
        a64206_100k = struct2array(load("airfoils.mat","a64206_100k"));
        a64206_200k = struct2array(load("airfoils.mat","a64206_200k"));
        a64206_500k = struct2array(load("airfoils.mat","a64206_500k"));
        a64206_1000k = struct2array(load("airfoils.mat","a64206_1000k"));

        a20714_50k = struct2array(load("airfoils.mat","a20714_50k"));
        a20714_100k = struct2array(load("airfoils.mat","a20714_100k"));
        a20714_200k = struct2array(load("airfoils.mat","a20714_200k"));
        a20714_500k = struct2array(load("airfoils.mat","a20714_500k"));
        a20714_1000k = struct2array(load("airfoils.mat","a20714_1000k"));
    end
    properties (Dependent, Hidden) % Airfoil internal properties
        range % Data R# range
        a6l % Root lower data bound
        a6u % Root upper data bound
        a2l % Tip lower data bound
        a2u % Tip upper data bound
        C_m_root % Root pitching moment coefficient, []
        C_m_tip % Tip pitching moment coefficient, []
        xbar_ac_root % Root chord-normalized aerodynamic center, []
        xbar_ac_tip % Tip chord-normalized aerodynamic center, []
    end
    properties (Dependent) % Airfoil properties
        RE % R# for airfoil data, []
        alpha_ZL_root % Zero-lift AoA, rad
        alpha_ZL_tip % Zero-lift AoA, rad

        C_d0_root % Root airfoil parasitic drag coefficient, []
        C_d0_tip % Tip airfoil parasitic drag coefficient, []
        C_d0 % Airfoil parasitic drag coefficient, []

        C_l_root % Root airfoil lift coefficient, []
        C_l_tip % Tip airfoil lift coefficient, []

        C_m % Average airfoil pitching coefficient, []

        a0_root % C_L vs AoA curve slope, root, 1/rad
        a0_tip % C_L vs AoA curve slope, tip, 1/rad
        a0 % Airfoil lift correlation, 1/rad

        xbar_ac % Location of chord normalized aerodynamic center, []
    end

    properties % Inputs
        V(1,1) double {mustBePositive, mustBeFinite} = 100; % Free-stream air velocity, m/s
        h(1,1) double {mustBeFinite} = 1; % Altitude, m
        b(1,1) double {mustBePositive, mustBeFinite} = 20; % Wingspan, m
        S_ref(1,1) double {mustBePositive, mustBeFinite} = 10; % Reference area, m^2
        Lambda_LE(1,1) double {mustBeNonnegative, mustBeLessThanOrEqual(Lambda_LE, 90)} = 45; % Leading edge sweep, deg
        lambda(1,1) double {mustBeNonnegative, mustBeLessThanOrEqual(lambda,1)} = 0.5; % Taper ratio, []
        d(1,1) double {mustBePositive, mustBeFinite} = 1.5; % Fuselage diameter, m
        m(1,1) double {mustBePositive, mustBeFinite} = 20e3; % Aircraft mass, kg
        L_fus(1,1) double {mustBePositive, mustBeFinite} = 20; % Length of fuselage, m
        A_max(1,1) double {mustBePositive, mustBeFinite} = 5; % Maximum total cross-sectional area (zy plane), m^2
        A_fus_top(1,1) double {mustBePositive, mustBeFinite} = 10; % Top (yx plane) cross-sectional area of fuselage, m^2
        A_fus_side(1,1) double {mustBePositive, mustBeFinite} = 10; % Side (zx plane) cross-sectional area of fuselage, m^2
        d_canopy(1,1) double {mustBePositive, mustBeFinite} = 0.5; % Canopy frontal (yz plane) diameter, m
        MAC_tail(1,1) double {mustBePositive, mustBeFinite} = 1; % Vertical/horizontal stabilizer avg MAC, m
        L_canopy(1,1) double {mustBePositive, mustBeFinite} = 2; % Canopy z-length, m
        S_canard(1,1) double {mustBePositive, mustBeFinite} = 2; % Canard reference area, m^2
        S_flapped_TE(1,1) double {mustBeNonnegative, mustBeFinite} = 1.5; % Trailing edge device reference area, m^2
        S_flapped_LE(1,1) double {mustBeNonnegative, mustBeFinite} = 0.5; % Leading edge device reference area, m^2
    end

methods % Airfoil Data Setup
        function range = get.range(obj)
            range = obj.datasetup();
        end
        function a6l = get.a6l(obj)
            [~, a6l, ~, ~, ~, ~] = obj.datasetup();
        end
        function a6u = get.a6u(obj)
            [~, ~, a6u, ~, ~, ~] = obj.datasetup();
        end
        function a2l = get.a2l(obj)
            [~, ~, ~, a2l, ~, ~] = obj.datasetup();
        end
        function a2u = get.a2u(obj)
            [~, ~, ~, ~, a2u, ~] = obj.datasetup();
        end
        function RE = get.RE(obj)
            [~, ~, ~, ~, ~, RE] = obj.datasetup();
        end
    end
    methods % Airfoil Calculations
        function alpha_ZL_root = get.alpha_ZL_root(obj)
            [~, ~, il2] = zerocrossrate(obj.a2l(:,2));
            [~, ~, iu2] = zerocrossrate(obj.a2u(:,2));
            il2 = find(il2(2:end)==1)+1;
            iu2 = find(iu2(2:end)==1)+1;
            alpha_l2 = obj.a2l(il2, 1);
            alpha_u2 = obj.a2u(iu2, 1);
            alpha_ZL_root = deg2rad(interp1(obj.range, [alpha_l2, alpha_u2], obj.RE));
        end
        function alpha_ZL_tip = get.alpha_ZL_tip(obj)
            [~, ~, il] = zerocrossrate(obj.a6l(:,2));
            [~, ~, iu] = zerocrossrate(obj.a6u(:,2));
            il = find(il(2:end)==1)+1;
            iu = find(iu(2:end)==1)+1;
            alpha_l = obj.a6l(il, 1);
            alpha_u = obj.a6u(iu, 1);
            alpha_ZL_tip = deg2rad(interp1(obj.range, [alpha_l, alpha_u], obj.RE));
        end
        function Cd_root = get.C_d0_root(obj)
            l = interp1(obj.a2l(:,1),obj.a2l(:,3), rad2deg(obj.alpha_eff_root));
            u = interp1(obj.a2u(:,1),obj.a2u(:,3), rad2deg(obj.alpha_eff_root));

            Cd_root = interp1(obj.range,[l, u], obj.RE);
        end
        function Cd_tip = get.C_d0_tip(obj)
            l = interp1(obj.a6l(:,1),obj.a6l(:,3), rad2deg(obj.alpha_eff_tip));
            u = interp1(obj.a6u(:,1),obj.a6u(:,3), rad2deg(obj.alpha_eff_tip));
            Cd_tip = interp1(obj.range,[l, u], obj.RE);
        end
        function Cl_root = get.C_l_root(obj)
            l = interp1(obj.a2l(:,1),obj.a2l(:,2), rad2deg(obj.alpha_root));
            u = interp1(obj.a2u(:,1),obj.a2u(:,2), rad2deg(obj.alpha_root));
            Cl_root = interp1(obj.range,[l, u], obj.RE);
        end
        function Cl_tip = get.C_l_tip(obj)
            l = interp1(obj.a6l(:,1),obj.a6l(:,2), rad2deg(obj.alpha_tip));
            u = interp1(obj.a6u(:,1),obj.a6u(:,2), rad2deg(obj.alpha_tip));
            Cl_tip = interp1(obj.range,[l, u], obj.RE);
        end
        function Cm_root = get.C_m_root(obj)
            l = interp1(obj.a2l(:,1),obj.a2l(:,5), rad2deg(obj.alpha_root));
            u = interp1(obj.a2u(:,1),obj.a2u(:,5), rad2deg(obj.alpha_root));
            Cm_root = interp1(obj.range,[l, u], obj.RE);
        end
        function Cm_tip = get.C_m_tip(obj)
            l = interp1(obj.a6l(:,1),obj.a6l(:,5), rad2deg(obj.alpha_tip));
            u = interp1(obj.a6u(:,1),obj.a6u(:,5), rad2deg(obj.alpha_tip));
            Cm_tip = interp1(obj.range,[l, u], obj.RE);
        end
        function a0_root = get.a0_root(obj)
            % average the slope of Cl vs alpha curve between 0 and 5 deg
            [~, ~, il] = zerocrossrate(obj.a2l(:,1));
            [~, ~, iu] = zerocrossrate(obj.a2u(:,1)-5);
            iil = find(il(2:end)==1,1)+2;
            iiu = find(iu(2:end)==1,1)+1;
            % a0_l = mean(obj.a2l(iil:iiu,2)./obj.a2l(iil:iiu,1));
            % a0_u = mean(obj.a2u(iil:iiu,2)./obj.a2l(iil:iiu,1));
            a0_l = (obj.a2l(iiu,2)-obj.a2l(iil,2))/(obj.a2l(iiu,1)-obj.a2l(iil,1));
            a0_u = (obj.a2u(iiu,2)-obj.a2u(iil,2))/(obj.a2u(iiu,1)-obj.a2u(iil,1));
            
            a0_root = interp1(obj.range, [a0_l, a0_u], obj.RE)*180/pi;
        end
        function a0_tip = get.a0_tip(obj)
            % average the slope of Cl vs alpha curve between 0 and 5 deg
            [~, ~, il] = zerocrossrate(obj.a6l(:,1));
            [~, ~, iu] = zerocrossrate(obj.a6u(:,1)-5);
            iil = find(il(2:end)==1,1)+2;
            iiu = find(iu(2:end)==1,1)+1;
            % a0_l = mean(obj.a6l(iil:iiu,2)./obj.a6l(iil:iiu,1));
            % a0_u = mean(obj.a6u(iil:iiu,2)./obj.a6u(iil:iiu,1));
            a0_l = (obj.a2l(iiu,2)-obj.a2l(iil,2))/(obj.a2l(iiu,1)-obj.a2l(iil,1));
            a0_u = (obj.a2u(iiu,2)-obj.a2u(iil,2))/(obj.a2u(iiu,1)-obj.a2u(iil,1));

            a0_tip = interp1(obj.range, [a0_l, a0_u], obj.RE)*180/pi;
        end
        function C_lmax_root = get.C_l_max_root(obj)
            C_lmax_root = interp1(obj.range, [max(obj.a2l(:,2)), max(obj.a2u(:,2))], obj.RE);
        end
        function C_lmax_tip = get.C_l_max_tip(obj)
            C_lmax_tip = interp1(obj.range, [max(obj.a6l(:,2)), max(obj.a6u(:,2))], obj.RE);
        end
        
        function xbar_ac_root = get.xbar_ac_root(obj)
            daoa = 0.25;
            Cll = obj.a2l(:,2);
            Clu = obj.a2u(:,2);
            Cml = obj.a2l(:,5);
            Cmu = obj.a2u(:,5);
            Clla = zeros(1,length(Cll)-1);
            Clua = Clla; Cmla = Clla; Cmua = Clla;

            for i = 1:length(Cll)-1
                Clla(i) = ( Cll(i+1)-Cll(i) )/daoa;
                Clua(i) = ( Clu(i+1)-Clu(i) )/daoa;
                Cmla(i) = ( Cml(i+1)-Cml(i) )/daoa;
                Cmua(i) = ( Cmu(i+1)-Cmu(i) )/daoa;
            end

            xbar_acl = -mean(Cmla./Clla, "all") + 0.25;
            xbar_acu = -mean(Cmua./Clua, "all") + 0.25;
            xbar_ac_root = interp1(obj.range, [xbar_acl, xbar_acu], obj.RE);
        end
        function xbar_ac_tip = get.xbar_ac_tip(obj)
            daoa = 0.25;
            Cll = obj.a6l(:,2);
            Clu = obj.a6u(:,2);
            Cml = obj.a6l(:,5);
            Cmu = obj.a6u(:,5);
            Clla = zeros(1, length(Cll)-1);
            Clua = Clla; Cmla = Clla; Cmua = Clla;

            for i = 1:length(Cll)-1
                Clla(i) = ( Cll(i+1)-Cll(i) )/daoa;
                Clua(i) = ( Clu(i+1)-Clu(i) )/daoa;
                Cmla(i) = ( Cml(i+1)-Cml(i) )/daoa;
                Cmua(i) = ( Cmu(i+1)-Cmu(i) )/daoa;
            end

            xbar_acl = -mean(Cmla./Clla, "all") + 0.25;
            xbar_acu = -mean(Cmua./Clua, "all") + 0.25;

            xbar_ac_tip = interp1(obj.range, [xbar_acl, xbar_acu], obj.RE);
        end
    end
    methods % Average Airfoil
        function Cd = get.C_d0(obj)
            Cd = (obj.C_d0_root+obj.C_d0_tip)/2;
        end
        function Cm = get.C_m(obj)
            Cm = (obj.C_m_root + obj.C_m_tip)/2;
        end
        function a0 = get.a0(obj)
            a0 = (obj.a0_root + obj.a0_tip)/2;
        end
        function xbar_ac = get.xbar_ac(obj)
            xbar_ac = (obj.xbar_ac_root + obj.xbar_ac_tip)/2;
        end
    end

    methods % Wing Methods
        function alpha = get.alpha_root(obj)
            alpha = obj.C_L/obj.C_L_alpha + obj.alpha_ZL_root;
        end
        function alpha = get.alpha_tip(obj)
            alpha = obj.C_L/obj.C_L_alpha + obj.alpha_ZL_tip;
        end
        function alpha_eff_root = get.alpha_eff_root(obj)
            alpha_eff_root = obj.alpha_root - obj.alpha_i;
        end
        function alpha_eff_tip = get.alpha_eff_tip(obj)
            alpha_eff_tip = obj.alpha_tip - obj.alpha_i;
        end

        function C_L_alpha = get.C_L_alpha(obj)
            F = 1.07*(1+obj.d/obj.b)^2;
            if obj.M < 1
                beta = sqrt(1-obj.M^2);
                eta = obj.a0*beta/(2*pi);
                C_L_alpha = (2*pi*obj.A)/( 2+sqrt(4 + obj.A^2*beta^2/eta^2*( 1 + tand(obj.Lambda_maxt)^2/beta^2 )) )*obj.S_exposed/obj.S_ref*F;
            else
                if obj.M > 1/cosd(obj.Lambda_LE)
                    % warning(['Supersonic lift correction outside of bounds. Lambda_LE: ', num2str(obj.Lambda_LE), '. Mach: ', num2str(obj.M)])
                    % C_L_alpha = NaN;
                    C_L_alpha = obj.a0/(1+obj.a0/(pi*obj.A*obj.e));
                    return
                end
                beta = sqrt(obj.M^2-1);
                factor = beta/tand(obj.Lambda_LE);
                vv = obj.A*tand(obj.Lambda_LE);

                

                AtanLLEx = [0.25, 0.5, 2, 3, 4, 5, 6];
                lambdax = [0, 0.2, 0.25, 1/3, 0.5, 1];

                for i = 1:length(AtanLLEx)-1
                    if vv < AtanLLEx(1)
                        i = 1; %#ok<FXSET>
                        break
                    % elseif vv >= AtanLLEx(end)
                    %     i = size(obj.A_LHS,1); %#ok<FXSET>
                    %     break
                    end
                    if vv > AtanLLEx(i) && vv < AtanLLEx(i+1)
                        break
                    end
                end
                for j = 1:length(lambdax)-1
                    if obj.lambda <  lambdax(1)
                        j = 1;
                        break
                    % elseif obj.lambda >= lambdax(end)
                    %     j = length(lambdax)-1;
                    %     break
                    end
                    if obj.lambda > lambdax(j) && obj.lambda < lambdax(j+1)
                        break
                    end
                end

                if factor > 1
                    factor = 1/factor;
                    factorx = 1:-0.2:0;

                    lower = interp1(factorx, obj.A_RHS(i, :, j), factor);
                    upper = interp1(factorx, obj.A_RHS(i+1, :, j), factor);
                    ml = interp1([AtanLLEx(i),AtanLLEx(i+1)], [lower, upper], vv);

                    lower = interp1(factorx, obj.A_RHS(i, :, j+1), factor);
                    upper = interp1(factorx, obj.A_RHS(i+1, :, j+1), factor);
                    mu = interp1([AtanLLEx(i),AtanLLEx(i+1)], [lower, upper], vv);

                    betaCL = interp1([lambdax(j),lambdax(j+1)], [ml, mu], obj.lambda);
                    C_L_alpha = betaCL/beta;

                else
                    factorx = 0:0.2:1;

                    lower = interp1(factorx, obj.A_LHS(i, :, j), factor);
                    upper = interp1(factorx, obj.A_LHS(i+1, :, j), factor);
                    ml = interp1([AtanLLEx(i),AtanLLEx(i+1)], [lower, upper], vv);

                    lower = interp1(factorx, obj.A_LHS(i, :, j+1), factor);
                    upper = interp1(factorx, obj.A_LHS(i+1, :, j+1), factor);
                    mu = interp1([AtanLLEx(i),AtanLLEx(i+1)], [lower, upper], vv);

                    tanLCL = interp1([lambdax(j),lambdax(j+1)], [ml, mu], obj.lambda);
                    C_L_alpha = tanLCL/tand(obj.Lambda_LE);

                end

                C_L_alpha = C_L_alpha*obj.S_exposed/obj.S_ref*F; %M2p.37

            end
        end

        function alpha = get.alpha(obj)
            alpha = mean([obj.alpha_root,obj.alpha_tip]);
        end

        function Re = get.Re(obj)
            Re = obj.V*obj.MAC/obj.nu;
        end

        function nu = get.nu(obj)
            nu = atmosphere.nu(obj.h);
        end

        function MAC = get.MAC(obj)
            MAC = obj.S_ref/obj.b;
            % MAC = 2*obj.C_root/3*(1+obj.lambda+obj.lambda^2)/(1+obj.lambda);
        end
        function C_r = get.C_root(obj)
            C_r = 2*obj.S_ref/(obj.b*(1+obj.lambda));
        end
        function A = get.A(obj)
            A = obj.b^2/(obj.S_ref+obj.S_canard);
        end

        function Lambda_025c = get.Lambda_025c(obj)
            Lambda_025c = obj.Sweeptf(obj.Lambda_LE, 0, 25, obj.A, obj.lambda);
        end
        function Lambda_maxt = get.Lambda_maxt(obj)
            Lambda_maxt = obj.Sweeptf(obj.Lambda_LE, 0, obj.xcm*100, obj.A, obj.lambda);
        end
        function Lambda_TE = get.Lambda_TE(obj)
            Lambda_TE = obj.Sweeptf(obj.Lambda_LE, 0, 100, obj.A, obj.lambda);
        end
        function Lambda_05c = get.Lambda_05c(obj)
            Lambda_05c = obj.Sweeptf(obj.Lambda_LE, 0, 50, obj.A, obj.lambda);
        end

        function M = get.M(obj)
            M = atmosphere.M(obj.V,obj.h);
        end

        function e = get.e(obj)
            Dlambda = -0.357 + 0.45*exp(0.0375*obj.Lambda_025c);
            e_theo = (1+obj.A*(0.0524*(obj.lambda-Dlambda)^4-0.15*(obj.lambda-Dlambda)^3+0.1659*(obj.lambda-Dlambda)^2-0.0706*(obj.lambda-Dlambda)+0.0119))^-1;

            k_eF = 1-2*(obj.d/obj.b)^2;
            k_eD0 = 0.873;

            if obj.M > 0.8
                k_eM = -0.001521*(obj.M/0.8-1)^10.82+1;
            else
                k_eM = 1;
            end

            e = e_theo*k_eF*k_eD0*k_eM; % Scholz correlation for oswald efficiency
        end

        function C_L = get.C_L(obj)
            C_L = obj.m*const.g/(obj.q*obj.S_ref);
        end

        function LD = get.LD(obj)
            LD = obj.C_L/obj.C_D;
        end
        
        function q = get.q(obj)
            q = atmosphere.q(obj.V,atmosphere.rho(obj.h));
        end

        function C_Di = get.C_Di(obj)
            C_Di = obj.K*obj.C_L^2;
        end

        function C_D = get.C_D(obj)
            C_D = obj.C_Di + obj.C_D0;
        end
        
        function K = get.K(obj)
            f = @(M) obj.A*(M^2-1)*cosd(obj.Lambda_LE)/(4*obj.A*sqrt(M^2-1)-2); % M2,p.50
            [M_min, K_min] = fminbnd(f, 1, 1.5);
            if obj.M <= 1
                K = 1/(pi*obj.A*obj.e);
            elseif obj.M < M_min
                K = 1/(pi*obj.A*obj.e);
                K = interp1([1, M_min],[K, K_min], obj.M); % Linear interpolation over broken area on K curve
            else
                K = f(obj.M);
            end
        end

        function S_exposed = get.S_exposed(obj)
            S_exposed = obj.S_ref - obj.d*2*obj.S_ref/(obj.b*(1+obj.lambda));
        end

        function S_wet = get.S_wet(obj)
            if obj.tc < 0.05
                S_wet = 2.003*obj.S_exposed + obj.S_wet_fus;
            else
                S_wet = (1.977 + 0.52*obj.tc)*obj.S_exposed + obj.S_wet_fus;
            end
        end

        function S_wet_fus = get.S_wet_fus(obj)
            S_wet_fus = 1.7*(obj.A_fus_top + obj.A_fus_side);
        end

        function C_D_misc = get.C_D_misc(obj) %#ok<MANU>
            C_D_misc = 0; % fix me, M2 p.46
        end

        function alpha_i = get.alpha_i(obj) % Induced AoA, rad
            alpha_i = obj.C_L/(pi*obj.A)*(2-1/obj.e);
        end

        function C_D0 = get.C_D0(obj) % Parasitic drag coefficient
            % Components: Fuselage, external stores, cockpit, vertical stabilizer

            % ASSUMPTIONS (FIX ME!)
            % Assuming cockpit is 20% the length of the fuselage, vertical
            % stabilizer is 30% of the length of the fuselage, assuming (big
            % assumption) that the diameter of the frontal area of the fuselage is
            % 15% of the length of the fuselage, cockpit frontal area is 5% of the
            % length of the fuselage, NACA0012 airfoil for tail using a 30 degree
            % sweep, wetted surface is a complete guess. No miscellaneous drag and
            % leakage and protuberance is 10%.
            C_f = [
                obj.C_f(obj.L_fus);
                6*obj.C_f(obj.L_stores);
                obj.C_f(obj.L_canopy);
                3*obj.C_f(obj.MAC_tail)
                ];
            FF = [
                obj.FF_fuselage_canopy(obj.L_fus, pi/4*obj.d^2);
                6*obj.FF_store(obj.L_stores, pi/4*0.35^2);
                obj.FF_fuselage_canopy(obj.L_canopy, pi/4*obj.d_canopy^2);
                3*obj.FF_tail()
                ];
            Q = [
                1; 1.3; 1; 1.04
                ]; % Component interference coefficients, p.45
            S_wet = [
                obj.S_wet_fus, 0.05*obj.S_wet, 0.08*obj.S_wet, 0.09*obj.S_wet
                ]; %#ok<PROP>
            if obj.M < 1
                tot = sum(C_f.*FF.*Q.*S_wet,"all"); %#ok<PROP>
            else
                tot = sum(C_f.*S_wet,"all"); %#ok<PROP>
            end
            C_D0 = tot/obj.S_ref + obj.C_D_misc + obj.C_Dw + obj.C_d0;
            C_D0 = 0.004*obj.S_wet/obj.S_ref;
            C_D0 = 1.1*C_D0; % Accounting for leakage and protuberance drag
        end

        function C_Dw = get.C_Dw(obj) % Wave drag coefficient
            Dq_SH = 9*pi/2*(obj.A_max/obj.L_fus)^2; % Sears-Haack correlation
            Dq_w = obj.E_WD*(1-0.2*(obj.M-1.2)^0.57)*(1-pi*obj.Lambda_LE^0.77/100)*Dq_SH;
            B = aerodynamics.E_WD*(1-pi*obj.Lambda_LE^0.77/100)*Dq_SH;
            if obj.M < obj.M_cr
                C_Dw = 0;
            elseif obj.M < obj.M_DD
                C_Dw = 0.002/(obj.M_DD-obj.M_cr)*obj.M-0.002/(obj.M_DD-obj.M_cr)*obj.M_cr;
            elseif obj.M < 1
                C_Dw = (B/2-0.002)/(1-obj.M_DD)*obj.M-(0.5*B-0.002)/(1-obj.M_DD)*obj.M_DD+0.002;
            elseif obj.M < 1.05
                C_Dw = 0.5*B/0.05*obj.M-0.5*B/0.05+0.5*B;
            elseif obj.M < 1.2
                C_Dw = B;
            else
                C_Dw = Dq_w/obj.S_ref;
            end
        end
        
        function M_cr = get.M_cr(obj) % Cruise mach
            M_cr = obj.M_DD-0.08; % Assuming the aircraft cruises at a similar relation as airliners
        end

        function M_DD = get.M_DD(obj) % Drag-divergence mach
            M_DD = 0.004*obj.Lambda_025c + 0.72;
        end

        function C_L_max = get.C_L_max(obj)
            C_L_max_root = obj.C_L_max_f(obj.C_l_max_root, obj.deltay_root);
            C_L_max_tip = obj.C_L_max_f(obj.C_l_max_tip, obj.deltay_tip);
            C_L_max = mean([C_L_max_root, C_L_max_tip]);
        end

        function alpha_C_L_max = get.alpha_C_L_max(obj)
            alpha_C_L_max_root = obj.alpha_C_L_max_f(obj, obj.alpha_ZL_root, obj.deltay_root);
            alpha_C_L_max_tip = obj.alpha_C_L_max_f(obj, obj.alpha_ZL_tip, obj.deltay_tip);
            alpha_C_L_max = min([alpha_C_L_max_root, alpha_C_L_max_tip]);
        end
        
    end

    methods (Access=private) % Internal Functions
        function C_f = C_f(obj, L) % Friction coefficient, a function of length L, m
            if obj.M < 1
                Recut = 38.21*(L/obj.k)^1.053;
            else
                Recut = 44.62*(L/obj.k)^1.053*obj.M^1.16;
            end
            ReCf = min([Recut, obj.V*L/obj.nu]);
            cf_lam = 1.328/ReCf^0.5;
            cf_turb = 0.074/ReCf^0.2;
            lamfrac = 0; % Zero laminar flow for military aircraft with camouglage
            C_f = cf_lam*lamfrac + cf_turb*(1-lamfrac);
        end

        function FF = FF_fuselage_canopy(obj, l, A_max) %#ok<INUSD> % Friction factor of fuselage/canopy
            FF = 0.9+5/l*sqrt(4/pi*A_max)+l/sqrt(4/pi*A_max)/400;
        end

        function FF = FF_tail(obj) % Friction factor of tail
            FF = (1+0.6/obj.xcm_tail*obj.tc_tail+100*obj.tc_tail^4)*(1.34*obj.M^0.18*cosd(obj.Lambda_maxt_tail)^0.28);
        end

        function FF = FF_store(obj, l, A_max) %#ok<INUSD> % Friction factor of stores
            FF = 1 + 0.35/l*sqrt(4/pi*A_max);
        end

        function C_L_max = C_L_max_f(obj, C_l_max, deltay) % Maximum lift coefficient
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
                if deltay <= dyx(1)
                    i = 1; %#ok<FXSET>
                    break
                elseif deltay >= dyx(end)
                    i = size(vals,1); %#ok<FXSET>
                    break
                end
                if deltay <= dyx(i+1) && deltay >= dyx(i)
                    break
                end
            end
            l = interp1(LLEx, vals(i,:), obj.Lambda_LE);
            if i == size(vals, 1)
                CLCl = l;
            else
                u = interp1(LLEx, vals(i+1,:), obj.Lambda_LE);
                CLCl = interp1(dyx(i:i+1), [l, u], deltay);
            end

            Lambda_HL_TE = obj.Sweeptf(obj.Lambda_LE, 0, 5, obj.A, obj.lambda); % Sweep of trailing edge hinge line, deg
            Lambda_HL_LE = obj.Sweeptf(obj.Lambda_LE, 0, 85, obj.A, obj.lambda); % Sweep of leading edge hinge line, deg

            DC_Lmax_TE = 0.9*obj.DC_l_max_TE*obj.S_flapped_TE/obj.S_ref*cosd(Lambda_HL_TE);
            DC_Lmax_LE = 0.9*obj.DC_l_max_LE*obj.S_flapped_LE/obj.S_ref*cosd(Lambda_HL_LE);
            DC_Lmax = DC_Lmax_TE + DC_Lmax_LE;

            C_L_max = C_l_max*CLCl + DC_Lmax;
        end
        
        function alpha_C_L_max = alpha_C_L_max_f(obj, alpha_ZL, deltay) % Stall angle, rad
            vals = [
                1.9, 2.1, 3.3, 4.5, 7.4, 10, 14;
                0.2, 1, 2.4, 3.9, 5.7, 7.5, 9.8;
                1.2, 1.7, 2.4, 3.2, 4.2, 5.3, 6.6;
                2.3, 2, 2.1, 2.3, 2.5, 2.9, 3.2
                ];

            dyx = [1.2, 2, 3, 4];
            LLEx = 0:10:60;

            for i = 1:size(vals, 1)
                if deltay <= dyx(1)
                    i = 1; %#ok<FXSET>
                    break
                elseif deltay >= dyx(end)
                    i = size(vals,1); %#ok<FXSET>
                    break
                end
                if deltay <= dyx(i+1) && deltay >= dyx(i)
                    break
                end
            end
            l = interp1(LLEx, vals(i,:), obj.Lambda_LE);
            if i == size(vals, 1)
                Dalpha_C_L_max = l;
            else
                u = interp1(LLEx, vals(i+1,:), rad2deg(obj.Lambda_LE));
                Dalpha_C_L_max = interp1(dyx(i:i+1), [l, u], deltay);
            end

            alpha_C_L_max = obj.C_L_max/obj.C_L_alpha + alpha_ZL + deg2rad(Dalpha_C_L_max); % something highly suspect is going on... AoA_stall above 20deg lol
        end
        function [range, a6l, a6u, a2l, a2u, RE] = datasetup(obj) % Airfoil data setup
            RE = obj.Re;
            if RE < 50e3
                RE = 50e3;
                range = [50e3, 50.01e3];
                a6l = obj.a64206_50k;
                a6u = obj.a64206_50k;
                a2l = obj.a20714_50k;
                a2u = obj.a20714_50k;

            elseif RE < 100e3
                range = [50e3, 100e3];
                a6l = obj.a64206_50k;
                a6u = obj.a64206_100k;
                a2l = obj.a20714_50k;
                a2u = obj.a20714_100k;
                

            elseif RE < 200e3
                range = [100e3, 200e3];
                a6l = obj.a64206_100k;
                a6u = obj.a64206_200k;
                a2l = obj.a20714_100k;
                a2u = obj.a20714_200k;

            elseif RE < 500e3
                range = [200e3, 500e3];
                a6l = obj.a64206_200k;
                a6u = obj.a64206_500k;
                a2l = obj.a20714_200k;
                a2u = obj.a20714_500k;
            elseif RE < 1000e3
                range = [500e3, 1000e3];
                a6l = obj.a64206_500k;
                a6u = obj.a64206_1000k;
                a2l = obj.a20714_500k;
                a2u = obj.a20714_1000k;
            else
                range = [1000e3, 1001e3];
                a6l = obj.a64206_1000k;
                a6u = obj.a64206_1000k;
                a2l = obj.a20714_1000k;
                a2u = obj.a20714_1000k;
                RE = 1000e3;
            end
        end
    end

    methods (Static) % Normal Functions
        function sweepn = Sweeptf(Lambda_m, m, n, A, lambda)% Sweep angle (deg) of y% chord from sweep (deg) at x% chord, AR, and taper
            sweepn = atand( tand(Lambda_m) - 4*(n-m)*(1-lambda)/(A*100*(1+lambda)) );
        end
    end

end