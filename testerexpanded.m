% close all
clear
tic

% Another comment

%% Constants

AIM120C = 348/2.2; % AIM-120C weight, kg
AIM9X = 186/2.2; % AIM-9X weight, kg
MK83JDAM = 1000/2.2;% JDAM weight, kg

m_p = 2500/2.2; % Payload mass, kg
m_c = 250; % Crew mass, kg
t_ll = 20*60; % Loiter to landing time, s
K_rf = 1.06; % Fuel reserve factor, []
dogfight.C_lo = 34*const.g*1e-6;
dogfight.m_d = 6*AIM120C+2*AIM9X;
strike.C_lo = 53*const.g*1e-6;
strike.m_d = 4*MK83JDAM+2*AIM9X;
C_cr = 42*const.g*1e-6;
C_ll = 21*const.g*1e-6;
h_cr = 30e3*0.3048; % Cruise altitude, m
h_ll = 10e3*0.3048; % Loiter-to-landing altitude, m
h_lo = 20e3*0.3048; % Dogfight loiter altitude, m


% Design Constants

% WILD ASSUMPTIONS!
A_max = 4; % Maximum cross-sectional area, m^2 (assumption)
A_fus_top = 10; % Fuselage top area, m^2 (assumption)
A_fus_side = 10; % Fuselage side area, m^2 (assumption)
d_canopy = 0.5; % Canopy hydraulic diameter, m (assumption)
L_canopy = 2; % Canopy length, m (assumption)


%% Solver

% Design Variables
n = 1;
nvar = 0;
N = n^nvar;

% b = linspace(10,18.288,n);
b = 18.288;
% lambda = linspace(0, 1, n); % Taper ratio, []
lambda = 0.88889;
% Lambda_LE = linspace(0, 60, n); % Leading edge sweep, deg
Lambda_LE = 46.667;

% Customer Variables
% M_dash_lo = linspace(0.85, 0.9, n); % Strike dash mach, []
M_dash_lo = 0.9;
% M_dash_cr = linspace(1.6, 2, n); % Dogfight dash mach, []
M_dash_cr = 1.7;
% psidot_sturn = linspace(8, 10, n); % Dogfight sustained turn rate, deg/s
psidot_sturn = 8.5;
% R = linspace(700*1852, 1000*1852, n); % Range, m
R = 1419.9e3;
% t_lo = linspace(2*60, 5*60, n); % Dogfight loiter time, s
t_lo = 5*60;
% Nz = linspace(7, 8, n); % Vertical load factor, []
Nz = 7;
% m_d = linspace(max([strike.m_d, dogfight.m_d]), 10e3, n); % Drop weight, kg
% m_d = 3000;
m_d = strike.m_d;

% % Guess Initializations
% m_0g = 20e3; % MTOW, kg
% K_g = 0.07; % Performance coefficient guess, []
% V_crg = 300; % Cruise speed guess, m/s
% C_L_max_0g = 1.9; % Takeoff maximum lift coefficient, []
% C_D0_crg = 0.008; % Cruise parasitic drag coefficient, []
% C_D0_maxV_crg = 0.013; % Cruise dash parasitic drag coefficient, []
% C_D0_maxV_lg = 0.011; % Loiter dash parasitic drag coefficient, []
% W_0crg = 0.9; % 0->Cruise weight fraction, []
% W_0mmg = 0.8; % 0->Mid-mission weight fraction, []
% W_0lg = 0.7; % 0->Landing weight fraction, []

%V,b,Lambda_LE,lambda
m = NaN.*ones(n,n,n,n);
LD = m;
Ro = m;

for i = 1:N
    i
    % subs = cell(1,nvar);
    % [subs{:}] = ind2sub(repmat(n,1,nvar), i);

    % ib = subs{1}; % b
    ib = 1;
    % il = subs{2}; % lambda
    il = 1;
    % iL = subs{3}; % Lambda_LE
    iL = 1;
    % iMl = subs{1}; % Loiter dash
    iMl = 1;
    % iMc = subs{2}; % Cruise dash
    iMc = 1;
    % ip = subs{3}; % Psidot
    ip = 1;
    % iR = subs{4}; % Range
    iR = 1;
    % it = subs{1}; % Dogfight loiter time
    it = 1;
    % iN = subs{2}; % Load factor
    iN = 1;
    % imd = subs{3}; % Drop weight
    imd = 1;

    strike.t_lo = 50*1852/atmosphere.V(M_dash_lo(iMl), 0); % 50nm dash at M_dash_lo
    
    % Guess Initializations
    m_0g = 20e3; % MTOW, kg
    K_g = 0.07; % Performance coefficient guess, []
    V_crg = 300; % Cruise speed guess, m/s
    C_L_max_0g = 1.9; % Takeoff maximum lift coefficient, []
    C_D0_crg = 0.008; % Cruise parasitic drag coefficient, []
    C_D0_maxV_crg = 0.013; % Cruise dash parasitic drag coefficient, []
    C_D0_maxV_lg = 0.011; % Loiter dash parasitic drag coefficient, []
    W_0crg = 0.9; % 0->Cruise weight fraction, []
    W_0mmg = 0.8; % 0->Mid-mission weight fraction, []
    W_0lg = 0.7; % 0->Landing weight fraction, []

    % Set up loop
    j = 0;
    tol = 0.001;
    rmsre = 2*tol;
    ers = 100;
    Derr = tol + 1;
    while Derr > tol && j < 1e2
        j = j + 1;
        if mod(j,10) == 0
            j
        end
        % Performance
        % if m_0g*2.2 > 100e3 || m_0g*W_0lg*2.2 > 70e3
        %     warning('Aircaft too heavy')
        % end
        if m_0g*2.2 > 2*90e3
            break
        end

        warning('off','all')
        [S_ref, T_0, V_maxturn, V_cr, V_be, V_s, WS_0, TW_0] = performance.fighter(psidot_sturn(ip), Nz(iN)/1.5, M_dash_cr(iMc),...
            M_dash_lo(iMl), V_crg, 74.6/1.1, K_g, m_0g, C_L_max_0g, W_0lg, C_D0_crg, W_0crg, h_cr, W_0mmg,...
            C_D0_maxV_crg, C_D0_maxV_lg, false);
        warning('on','all')
        if isnan(S_ref)
            break
        end
        % end
        M_cr = atmosphere.M(V_cr,h_cr);

        % Component Sizing
        subsystems = subsystems();
        subsystems.m = m_0g;
        subsystems.MAC = b(ib)/S_ref;
        subsystems.b = b(ib);
        subsystems.S_ref = S_ref;

        % Aerodynamics
        aero = aa();
        aero.b = b(ib);
        aero.S_ref = S_ref;
        aero.Lambda_LE = Lambda_LE(iL);
        aero.lambda = lambda(il);
        aero.d = subsystems.d;
        aero.m = m_0g;
        aero.L_fus = subsystems.L_fus;
        aero.A_max = A_max;
        aero.A_fus_top = A_fus_top;
        aero.A_fus_side = A_fus_side;
        aero.d_canopy = d_canopy;
        aero.MAC_tail = 0.3*b(ib)/(subsystems.S_ref_VT+subsystems.S_ref_C); % Assuming tail span is 30% of wingspan. Also assuming not canard i guess (fix me)
        aero.L_canopy = L_canopy;
        aero.S_canard = subsystems.S_ref_C;
        aero.S_flapped_TE = 0.15*S_ref; % Assuming flapped area is 15% of ref
        aero.S_flapped_LE = 0.1*S_ref; % Assuming slat area is 10% of ref

        % Cruise
        if imag(V_cr) ~= 0
            break
        end
        aero.V = V_cr;
        aero.h = h_cr;
        aero.m = m_0g*W_0crg;
        LD_cr = aero.LD;
        K = aero.K;
        C_D0_cr = aero.C_D0;

        % Loiter to landing
        aero.V = V_be;
        aero.h = h_ll;
        aero.m = m_0g*W_0lg;
        LD_ll = aero.LD;

        % Dogfight Loiter
        aero.V = atmosphere.V(M_dash_cr(iMc),h_lo);
        aero.h = h_lo;
        aero.m = m_0g*W_0mmg;
        dogfight.LD_lo = aero.LD;
        C_D0_maxV_cr = aero.C_D0;

        % Strike Loiter
        aero.V = atmosphere.V(M_dash_lo(iMl),0);
        aero.h = 0;
        aero.m = m_0g*W_0mmg;
        strike.LD_lo = aero.LD;
        C_D0_maxV_l = aero.C_D0;

        % Takeoff
        aero.V = 1.1*V_s;
        aero.h = 0;
        C_L_max_0 = aero.C_L_max;

        if isnan(LD_cr) || isnan(LD_ll) || isnan(dogfight.LD_lo) || isnan(strike.LD_lo)
            break
        end

        % Weights
        [m_0df, m_fdf, m_fidf, W_0crdf, W_0mmdf, W_0ldf] = weight.weights(m_0g, m_p, m_d(imd), m_c, K_rf, dogfight.LD_lo,...
            dogfight.C_lo, M_dash_cr(iMc), t_lo(it), LD_cr, C_cr, M_cr, V_cr, R(iR), C_ll, LD_ll, t_ll,...
            aero.A, TW_0, WS_0);
        [m_0s, m_fs, m_fis, W_0crs, W_0mms, W_0ls] = weight.weights(m_0g, m_p, m_d(imd), m_c, K_rf, strike.LD_lo,...
            strike.C_lo, M_dash_lo(iMl), strike.t_lo, LD_cr, C_cr, M_cr, V_cr, R(iR), C_ll, LD_ll, t_ll,...
            aero.A, TW_0, WS_0);
        m_0 = max([m_0df, m_0s]);
        adfiojdsaof = size(m_0);
        if m_fdf > m_fs
            m_f = m_fdf;
            m_fi = m_fidf;
            W_0cr = W_0crdf;
            W_0mm = W_0mmdf;
            W_0l = W_0ldf;
        else
            m_f = m_fs;
            m_fi = m_fis;
            W_0cr = W_0crs;
            W_0mm = W_0mms;
            W_0l = W_0ls;
        end
        if isnan(m_0) || isnan(sum(m_fi)) || isnan(W_0cr) || isnan(W_0mm) || isnan(W_0l)
            break
        end

        % Subsystems
        subsystems.m_f = m_f;
        subsystems.D_to = aero.C_D*aero.q*aero.S_ref;
        subsystems.T0 = T_0;
        subsystems.M_ac = aero.C_m*aero.q*aero.MAC^2*b(ib);
        subsystems.a_to = (performance.V_end(m_0) + sqrt(2*TW_0*performance.l_to))^2/(2*performance.l_to);

        subsystems.L_to = aero.C_L*aero.q*aero.S_ref;
        subsystems.xbar_ac = aero.xbar_ac;
        subsystems.alpha_to = rad2deg(aero.alpha);
        subsystems.c_root = aero.C_root;
        subsystems.Lambda_05c = aero.Lambda_05c;
        subsystems.Lambda_LE = Lambda_LE(iL);
        subsystems.lambda = lambda(il);
        subsystems.S_wet_fus = aero.S_wet_fus;
        subsystems.W0mm = W_0mm;
        subsystems.S_exposed = aero.S_exposed;

        subsystems.L_canard = aero.C_l_root*aero.q*subsystems.S_ref_C; % Assuming canard airfoil is the same as wing root airfoil, at its maximum lift

        % Calculate Error
        rmsre = sqrt(mean([
            (m_0g - m_0)/m_0;
            (K_g - K)/K;
            (V_crg - V_cr)/V_cr;
            (C_L_max_0g - C_L_max_0)/C_L_max_0;
            (C_D0_crg - C_D0_cr)/C_D0_cr;
            (C_D0_maxV_crg - C_D0_maxV_cr)/C_D0_maxV_cr;
            (C_D0_maxV_lg - C_D0_maxV_l)/C_D0_maxV_l;
            (W_0crg - W_0cr)/W_0cr;
            (W_0mmg - W_0mm)/W_0mm;
            (W_0lg - W_0l)/W_0l
            ].^2));
        ers(j+1) = rmsre; %#ok<SAGROW>
        Derr = abs(rmsre-ers(j))/ers(j);

        % Update Guesses
        m_0g = m_0;
        K_g = K;
        C_L_max_0g = C_L_max_0;
        C_D0_crg = C_D0_cr;
        C_D0_maxV_crg = C_D0_maxV_cr;
        C_D0_maxV_lg = C_D0_maxV_l;
        W_0crg = W_0cr;
        W_0mmg = W_0mm;
        W_0lg = W_0l;
    end

    m(ib, il, iL, iMl, iMc, ip, iR, it, iN, imd) = m_0g;
    LD(ib, il, iL, iMl, iMc, ip, iR, it, iN, imd) = LD_cr;
    Ro(ib, il, iL, iMl, iMc, ip, iR, it, iN, imd) = R(iR);


    [~, ~, ~, ~, ~, ~, WS, TW] = performance.fighter(psidot_sturn(ip), Nz(iN)/1.5, M_dash_cr(iMc),...
        M_dash_lo(iMl), V_crg, 74.6/1.1, K_g, m_0g, C_L_max_0g, W_0lg, C_D0_crg, W_0crg, h_cr, W_0mmg,...
        C_D0_maxV_crg, C_D0_maxV_lg, true)
    figure
    plot(m_f-cumsum(m_fi))
    ylabel('Fuel Stored, kg')
    xlabel('Mission Segment')
    % close all
end
toc
%% Sensitivity Analysis

%% Trade Studies

%{
% Loiter dash
figure
plot(M_dash_lo,reshape(LD(1,1,1,:,1,1,1,1,1,1),1,10))
ylabel('L/D, []')
yyaxis right
plot(M_dash_lo,reshape(m(1,1,1,:,1,1,1,1,1,1),1,10))
xlabel('M_{dash,lo}, []')
ylabel('m_0, kg')
% Cruise dash
figure
plot(M_dash_cr,reshape(LD(1,1,1,1,:,1,1,1,1,1),1,10))
xlabel('M_{dash,cr}, []')
ylabel('L/D, []')
yyaxis right
plot(M_dash_cr,reshape(m(1,1,1,1,:,1,1,1,1,1),1,10))
ylabel('m_0, kg')
% Psidot
figure
plot(psidot_sturn,reshape(LD(1,1,1,1,1,:,1,1,1,1),1,10))
xlabel('(d\psi /dt)_{sturn}, deg/s')
ylabel('L/D, []')
yyaxis right
plot(psidot_sturn,reshape(m(1,1,1,1,1,:,1,1,1,1),1,10))
ylabel('m_0, kg')
% Range
figure
plot(R,reshape(LD(1,1,1,1,1,1,:,1,1,1),1,10))
xlabel('range, m')
ylabel('L/D, []')
yyaxis right
plot(R,reshape(m(1,1,1,1,1,1,:,1,1,1),1,10))
ylabel('m_0, kg')
% Dogfight loiter time
figure
plot(t_lo,reshape(LD(1,1,1,1,1,1,1,:,1,1),1,10))
xlabel('loiter time, s')
ylabel('L/D, []')
yyaxis right
plot(t_lo,reshape(m(1,1,1,1,1,1,1,:,1,1),1,10))
ylabel('m_0, kg')
yline(90e3/2.2,'k--','Max m0')
% Load factor
figure
plot(Nz,reshape(LD(1,1,1,1,1,1,1,1,:,1),1,10))
xlabel('Load factor, []')
ylabel('L/D, []')
yyaxis right
plot(Nz,reshape(m(1,1,1,1,1,1,1,1,:,1),1,10))
ylabel('m_0, kg')
yline(90e3/2.2,'k--','Max m0')
% Drop weight
figure
plot(m_d,reshape(LD(1,1,1,1,1,1,1,1,1,:),1,10))
xlabel('Drop weight, kg')
ylabel('L/D, []')
yyaxis right
plot(m_d,reshape(m(1,1,1,1,1,1,1,1,1,:),1,10))
ylabel('m_0, kg')
yline(90e3/2.2,'k--','Max m0')
%}

%% Score

% close all
% L_fus < 50ft
% H_fus < 18.5ft

% Data Setup/Validation
mr = reshape(m(:, :, :, 1, 1, 1, :, 1, 1, 1), n, n, n, n);
mr2 = mr;
mr2(mr2>90e3/2.2) = NaN;
mr2(mr2<20e3/2.2) = NaN;

LDr = reshape(m(:, :, :, 1, 1, 1, :, 1, 1, 1), n, n, n, n);
LDr2 = LDr;
LDr2(mr2>90e3/2.2) = NaN;
LDr2(mr2<20e3/2.2) = NaN;
LDr2(LDr2>1e6) = NaN;

Ror = reshape(m(:, :, :, 1, 1, 1, :, 1, 1, 1), n, n, n, n);
Ror2 = Ror;
Ror2(mr2>90e3/2.2) = NaN;
Ror2(mr2<20e3/2.2) = NaN;
Ror2(LDr2>1e6) = NaN;
mr2(LDr2>1e6) = NaN;


% Scoring
Z = 0.3.*(max(mr2,[],"all")-mr2)./sum(mr2, 'all','omitnan') + 1.*LDr2./sum(LDr2, 'all','omitnan') + 1.5.*Ror2./sum(Ror2, 'all','omitnan');
[i, j, k, l] = ind2sub(size(Z), find(Z==max(Z, [], "all")));
% LDmax = maxk(reshape(LDr2,1,[]),10)

%%
% close all
mout = mr2(i, j, k, l);
bout = b(i);
lambdaout = lambda(j);
Lambdaout = Lambda_LE(k);
Rangeout = R(l);
MAX_SCORE = max(Z, [], "all");
table([mout; bout; lambdaout; Lambdaout; Rangeout; MAX_SCORE],'VariableNames',{'Value'},'RowNames',{'m_0','b','lambda','Lambda_LE','R','Z'})

figure
hold on
surf(b, lambda, mr2(:,:,1,1), Z(:,:,1,1))
surf(b, lambda, mr2(:,:,k,l), Z(:,:,k,l),'EdgeColor','none')
plot3(b(i),lambda(j),mr2(i,j,k,l),'k*')
hold off
c = colorbar;
c.Label.String = 'Score';
xlabel("b, m")
ylabel("\lambda")
zlabel("m_0, kg")

figure
hold on
surf(Lambda_LE, R, reshape(mr2(1,1,:,:), [10, 10]),reshape(Z(1,1,:,:), [10, 10]))
surf(Lambda_LE, R, reshape(mr2(i,j,:,:), [10, 10]),reshape(Z(i,j,:,:), [10, 10]),'EdgeColor','none')
plot3(Lambda_LE(k),R(l),mr2(i,j,k,l),'k*')
hold off
c = colorbar;
c.Label.String = 'Score';
xlabel("\Lambda_{LE}, deg")
ylabel("R, m")
zlabel("m_0, kg")

%%

% figure
% contourf(b, lambda, Z(:,:,1,1))
% c = colorbar;
% c.Label.String = "Z";
% xlabel('b')
% ylabel('\lambda')

% figure
% hold on
% surf(b, lambda, Z(:,:,1,1),'EdgeColor','none')
% surf(b, lambda, Z(:,:,round(size(Z,3)/2),round(size(Z,4)/2)))
% surf(b, lambda, Z(:,:,end,end))
% hold off
% xlabel('b')
% ylabel('\\lamdba')
% zlabel('Z')
% legend('1','mid','end')



figure
hold on
surf(b, lambda, Z(:,:,k,l),'EdgeColor','none')
surf(b, lambda, Z(:,:,round(size(Z,3)/2),round(size(Z,4)/2)))
surf(b, lambda, Z(:,:,1,1))
plot3(b(i),lambda(l),Z(i, j, k, l),'k*')
hold off
xlabel('b')
ylabel("\lambda")
zlabel('Z')
legend('max','mid','maximum')

figure
hold on
surf(Lambda_LE, R, reshape(Z(i, j, :, :), [n, n]),'EdgeColor','none')
surf(Lambda_LE, R, reshape(Z(1,1,:,:),[n, n]))
plot3(Lambda_LE(j),R(k),Z(i, j, k, l),'k*')
ylim([-Inf, 1600e3])
hold off
legend('Max','1','Max Score')
xlabel("\Lambda_{LE}, deg")
ylabel('R, m')
zlabel('Z')

% figure
% hold on
% contour(Lambda_LE, R, reshape(Z(i, j, :, :), [n, n]))
% plot(Lambda_LE(j),R(k),Z(i, j, k, l),'k*')
% hold off