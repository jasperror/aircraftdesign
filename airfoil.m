% ENGR 493 Intro to Aircraft Design
% Project Phase 2: Wing Design - Airfoil characteristics based on the
% NACA64-206 and the NASA SC2-0714
% Jasper Palmer
classdef airfoil % commented code? When its 2 hours until the deadline, we don't need commented code...
    properties
    end
    methods (Static)
        function loader()
            load('airfoils.mat') %#ok<LOAD> Data collected from airfoiltools.com
        end
        function [range, a6l, a6u, a2l, a2u] = datasetup(Re)
            if Re < 100e3
                range = [50e3, 100e3];
                a6l = a64206_50k;
                a6u = a64206_100k;
                a2l = a20714_50k;
                a2u = a20714_100k;

            elseif Re < 200e3
                range = [100e3, 200e3];
                a6l = a64206_100k;
                a6u = a64206_200k;
                a2l = a20714_100k;
                a2u = a20714_200k;

            elseif Re < 500e3
                range = [200e3, 500e3];
                a6l = a64206_200k;
                a6u = a64206_500k;
                a2l = a20714_200k;
                a2u = a20714_500k;
            else
                range = [500e3, 1000e3];
                a6l = a64206_500k;
                a6u = a64206_1000k;
                a2l = a20714_500k;
                a2u = a20714_1000k;
            end
        end
        function [alpha_ZL_root, alpha_ZL_tip] = zerolift(Re)
            airfoil.loader();
            [range, a6l, a6u, a2l, a2u] = airfoil.datasetup(Re);
            % alpha_ZL
            [~, ~, il6] = zerocrossrate(a6l(:,2));
            [~, ~, iu6] = zerocrossrate(a6u(:,2));
            [~, ~, il2] = zerocrossrate(a2l(:,2));
            [~, ~, iu2] = zerocrossrate(a2u(:,2));
            il6 = find(il6(2:end)==1)+1;
            iu6 = find(iu6(2:end)==1)+1;
            il2 = find(il2(2:end)==1)+1;
            iu2 = find(iu2(2:end)==1)+1;
            alpha_l6 = a6l(il6, 1);
            alpha_u6 = a6l(iu6, 1);
            alpha_l2 = a6l(il2, 1);
            alpha_u2 = a6l(iu2, 1);
            alpha_ZL_root = deg2rad(interp1(range, [alpha_l2, alpha_u2], Re));
            alpha_ZL_tip = deg2rad(interp1(range, [alpha_l6, alpha_u6], Re));
        end
        function [Cd_root, Cd_tip, Cl_root, Cl_tip, Cm_root, Cm_tip] = airfoil(alpha_root, alpha_tip, alpha_eff_root, alpha_eff_tip, Re)
            airfoil.loader();
            [range, a6l, a6u, a2l, a2u] = airfoil.datasetup(Re);
            alpha_tip = rad2deg(alpha_tip);
            alpha_root = rad2deg(alpha_root);
            % Cl
            l6 = interp1(a6l(:,1),a6l(:,2), alpha_tip);
            u6 = interp1(a6u(:,1),a6u(:,2), alpha_tip);
            l2 = interp1(a2l(:,1),a2l(:,2), alpha_root);
            u2 = interp1(a2u(:,1),a2u(:,2), alpha_root);
            Cl_root = interp1(range,[l2, u2], Re);
            Cl_tip = interp1(range, [l6, u6], Re);
            % Cd
            l6 = interp1(a6l(:,1),a6l(:,3), alpha_eff_tip);
            u6 = interp1(a6u(:,1),a6u(:,3), alpha_eff_tip);
            l2 = interp1(a2l(:,1),a2l(:,3), alpha_eff_root);
            u2 = interp1(a2u(:,1),a2u(:,3), alpha_eff_root);
            Cd_root = interp1(range,[l2, u2], Re);
            Cd_tip = interp1(range, [l6, u6], Re);
            % Cm
            l6 = interp1(a6l(:,1),a6l(:,5), alpha_tip);
            u6 = interp1(a6u(:,1),a6u(:,5), alpha_tip);
            l2 = interp1(a2l(:,1),a2l(:,5), alpha_root);
            u2 = interp1(a2u(:,1),a2u(:,5), alpha_root);
            Cm_root = interp1(range,[l2, u2], Re);
            Cm_tip = interp1(range, [l6, u6], Re);
        end
        function [a0_root, a0_tip] = a0(Re)
            airfoil.loader();
            [range, a6l, a6u, a2l, a2u] = airfoil.datasetup(Re);
            % average the slope of Cl vs alpha curve between 0 and 5 deg

            [~, ~, il6] = zerocrossrate(a6l(:,1));
            [~, ~, iu6] = zerocrossrate(a6u(:,1)-5);
            [~, ~, il2] = zerocrossrate(a2l(:,1));
            [~, ~, iu2] = zerocrossrate(a2u(:,1)-5);
            i6l = find(il6(2:end)==1)+1;
            i6u = find(iu6(2:end)==1)+1;
            i2l = find(il2(2:end)==1)+1;
            i2u = find(iu2(2:end)==1)+1;

            a0_6l = mean(a6l(i6l:i6u,2)./a6l(i6l:i6u,1));
            a0_6u = mean(a6u(i6l:i6u,2)./a6u(i6l:i6u,1));
            a0_2l = mean(a2l(i2l:i2u,2)./a2l(i2l:i2u,1));
            a0_2u = mean(a2u(i2l:i2u,2)./a2l(i2l:i2u,1));

            a0_root = interp1(range, [a0_2l, a0_2u], Re);
            a0_tip = interp1(range, [a0_6l, a0_6u], Re);

        end
        function [C_lmax_root, C_lmax_tip] = clmax(Re)
            airfoil.loader();
            [range, a6l, a6u, a2l, a2u] = airfoil.datasetup(Re);
            C_lmax_root = interp1(range, [max(a2l(:,2)), max(a2u(:,2))], Re);
            C_lmax_tip = interp1(range, [max(a6l(:,2)), max(a6u(:,2))], Re);
        end
    end
end