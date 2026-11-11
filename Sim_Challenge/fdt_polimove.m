%% fdt_polimove.m
% Funzioni di trasferimento linearizzate del veicolo PoliMOVE SimChallenge26.
% Legge i parametri da id_results.mat (prodotto da identify_polimove.m); se il file o un campo
% mancano usa i valori dell'ultimo log (riportati qui sotto).
% Richiede la Control System Toolbox (tf, ss, bode, step).
%
% Per ogni velocita' di linearizzazione v0 in V0 crea nel workspace:
%   LAT(k).Gr      r/delta            [(rad/s)/rad]  bicicletta 2 stati + ritardo sterzo
%   LAT(k).Gay     ay/delta           [(m/s^2)/rad]  ay = ACCELERAZIONE laterale (non l'errore)
%   LAT(k).Gepsi  e_psi/delta        [rad/rad]  errore di direzione = Gr/s   (riferimento dritto)
%   LAT(k).Gey     e_y/delta          [m/rad]    errore laterale = Gay/s^2    (riferimento dritto)
%   LAT(k).Gbeta   beta/delta         [rad/rad]
%   LAT(k).Gr_rel  r/delta con rilassamento pneumatici (4 stati) + ritardo
%   LAT(k).Gr_dp   r/delta come ritardo puro + polo singolo (modello del collega, da E6)
%   LON(k).Gped    dv/dpedale         [(m/s)/unita' di pedale]  + ritardo pedale
%   LON(k).Gbrk    dv/dfreno          [(m/s)/unita' di comando freni]  + attuatore freni
%   LON(k).GF      dv/dFx             [(m/s)/N]   (forza a terra, senza attuatori)
% e la struct BRK con saturazione e limiti di velocita' dell'attuatore freni, per i blocchi
% Saturation / Rate Limiter di Simulink (la parte NON lineare che le fdt non vedono).
%
% In Simulink: blocco "LTI System" con LAT(k).Gr, LON(k).Gbrk ecc. (il ritardo e' dentro come
% InputDelay), oppure "Transfer Fcn" con num/den + "Transport Delay".
% ATTENZIONE: sono fdt su VARIAZIONI attorno al punto di lavoro (v0, pedale0, freno 0):
% all'ingresso va sottratto il valore di regime, all'uscita va sommato v0.
clear; clc; close all;

%% 0) SCELTE ----------------------------------------------------------------
V0 = [20 40 60];          % [m/s] velocita' di linearizzazione (v max identificata ~72 m/s)
BIAS_F = 0.6;             % ripartizione del comando freni sull'anteriore (come nel nodo)
% Mappa del pedale: frazione della coppia (al netto del freno motore) vs pedale.
% La mappa stimata dall'ultimo log e' rumorosa e non monotona: uso quella misurata dal collega.
% Metti USE_PEDAL_FROM_RES = true per usare RES.pedal quando sara' pulita.
PED_P = [0 0.35 0.40 0.55 0.70 1.00];
PED_F = [0 0.00 0.14 0.48 0.78 1.00];
USE_PEDAL_FROM_RES = false;

%% 1) PARAMETRI: valori dell'ultimo log, sovrascritti da id_results.mat ------
P.m = 720;  P.a = 1.760;  P.b = 1.355;  P.L = P.a + P.b;  P.g = 9.81;  P.rho = 1.225;
P.ClfA = 1.832;  P.ClrA = 1.948;  P.Rw = 0.309;
P.Fz0f = P.m*P.g*P.b/P.L;  P.Fz0r = P.m*P.g*P.a/P.L;

% strada e freno motore (E1)
X.c0 = 39;  X.c2 = 0.7472;  X.e0 = 58.07;  X.e1 = 0.00020;
% powertrain (E2)
X.G  = [8.689 6.442 4.973 4.005 3.355 2.919];
X.knots = 1500:500:8000;
X.Te    = [80 89 95 104 121 141 171 221 319 384 380 350 354 359];
X.J  = 0.027;
X.upshift = [7922 7913 7915 7921 7932];
% freni (E3): pressione = gain * comando fino a psat; K per ruota [Nm/Pa]
X.brk_gain = 1e5;  X.psat = 1.2e7;  X.Kf = 5.771e-5;  X.Kr = 1.398e-4;
X.slew_up = 3.68e8;  X.slew_down = 9.26e7;     % [Pa/s]
% pneumatici (E5) [B C E mu0 mu1]
X.thf = [21.847 1.100 -1.673 1.357 -0.276];
X.thr = [21.400 1.100 -1.334 1.571 -0.154];
% imbardata (E6)
X.Iz = 668;  X.sigma = 0.455;  X.L_dp = 0.070;  X.tau_dp = 0.031;
% attuatori (check unita': ritardo comando -> feedback)
X.Td_steer = 0.06;  X.Td_ped = 0.06;  X.Td_brk = 0.07;

if isfile('id_results.mat')
    S = load('id_results.mat');  R = S.RES;
    fprintf('Parametri letti da id_results.mat\n');
    if isfield(R,'coast')
        X.c2 = R.coast.c2;  X.e0 = R.coast.e0;  X.e1 = R.coast.e1;
        if R.coast.c0 > 0, X.c0 = R.coast.c0;
        else, fprintf('  c0 = %.1f N non fisico (negativo): uso %.0f N\n', R.coast.c0, X.c0); end
    end
    if isfield(R,'ratios'), X.G(isfinite(R.ratios.G)) = R.ratios.G(isfinite(R.ratios.G)); end
    if isfield(R,'pt') && isfield(R.pt,'Te') && any(isfinite(R.pt.Te))
        ok = isfinite(R.pt.Te);  X.knots = R.pt.knots(ok);  X.Te = R.pt.Te(ok);
        X.J = max(R.pt.J, 0);
        if isfield(R.pt,'upshift_rpm'), u = R.pt.upshift_rpm(1:5); X.upshift(isfinite(u)) = u(isfinite(u)); end
    end
    if USE_PEDAL_FROM_RES && isfield(R,'pedal') && numel(R.pedal.p) >= 3
        PED_P = R.pedal.p;  PED_F = cummax(R.pedal.f);   % resa monotona
    end
    if isfield(R,'brkmap')
        if isfinite(R.brkmap.gain), X.brk_gain = R.brkmap.gain; end
        if isfinite(R.brkmap.psat), X.psat = R.brkmap.psat; end
        if isfield(R.brkmap,'slew_up') && isfinite(R.brkmap.slew_up)
            X.slew_up = R.brkmap.slew_up;  X.slew_down = R.brkmap.slew_down;
        end
    end
    if isfield(R,'brk')
        if isfinite(R.brk.Kf), X.Kf = R.brk.Kf; end
        if isfinite(R.brk.Kr), X.Kr = R.brk.Kr; end
    end
    if isfield(R,'tires'), X.thf = R.tires.front.th;  X.thr = R.tires.rear.th; end
    if isfield(R,'yaw'), X.Iz = R.yaw.Iz;  X.sigma = R.yaw.sigma; end
    if isfield(R,'delaypole') && ~isempty(R.delaypole)
        X.L_dp = median(R.delaypole(:,3));  X.tau_dp = median(R.delaypole(:,4));
    end
else
    fprintf('id_results.mat non trovato: uso i valori dell''ultimo log scritti nello script\n');
end
Te_of = @(n) interp1(X.knots, X.Te, min(max(n, X.knots(1)), X.knots(end)), 'linear');
Teb_of = @(n) X.e0 + X.e1*n;                         % freno motore lato motore [Nm]
f_of = @(p) interp1(PED_P, PED_F, min(max(p, 0), 1), 'linear');
[Fu, iu] = unique(PED_F, 'last');  Pu = PED_P(iu);  % inversa: per f = 0 prende il pedale piu' alto che da' 0
p_of = @(f) interp1(Fu, Pu, min(max(f, Fu(1)), Fu(end)), 'linear');

% attuatore freni (parte lineare): ritardo puro + primo ordine sulla salita
BRK.gain = X.brk_gain;  BRK.psat = X.psat;  BRK.cmd_sat = X.psat/X.brk_gain;
BRK.slew_up = X.slew_up;  BRK.slew_down = X.slew_down;
BRK.tau = X.psat/X.slew_up;          % costante di tempo equivalente in salita [s]
BRK.Td = X.Td_brk;
BRK.Kf = X.Kf;  BRK.Kr = X.Kr;  BRK.bias_f = BIAS_F;
% forza frenante a terra per unita' di comando totale u (anteriore BIAS_F*u, posteriore (1-BIAS_F)*u)
BRK.KF_per_cmd = 2*X.brk_gain/P.Rw * (X.Kf*BIAS_F + X.Kr*(1 - BIAS_F));    % [N per unita' di comando]

%% 2) COSTRUZIONE ------------------------------------------------------------
s = tf('s');
for k = 1:numel(V0)
    v = V0(k);  m = P.m;  a = P.a;  b = P.b;  Iz = X.Iz;
    % ---- LATERALE: rigidezze di deriva = pendenza della Pacejka (B*C*D) al carico a v
    Fzf = P.Fz0f + 0.5*P.rho*P.ClfA*v^2;
    Fzr = P.Fz0r + 0.5*P.rho*P.ClrA*v^2;
    Cf = X.thf(1)*X.thf(2)*(X.thf(4) + X.thf(5)*(Fzf - P.Fz0f)/P.Fz0f)*Fzf;
    Cr = X.thr(1)*X.thr(2)*(X.thr(4) + X.thr(5)*(Fzr - P.Fz0r)/P.Fz0r)*Fzr;
    A = [ -(Cf+Cr)/(m*v),        -v - (a*Cf - b*Cr)/(m*v);
          -(a*Cf - b*Cr)/(Iz*v), -(a^2*Cf + b^2*Cr)/(Iz*v) ];
    B = [ Cf/m; a*Cf/Iz ];
    C = [ 0 1;  A(1,1) A(1,2)+v;  1/v 0 ];
    D = [ 0; Cf/m; 0 ];
    sysl = ss(A, B, C, D, 'InputDelay', X.Td_steer, 'InputName', 'delta', 'OutputName', {'r','ay','beta'});
    Gl = tf(sysl);
    LAT(k).v = v;  LAT(k).Cf = Cf;  LAT(k).Cr = Cr;  LAT(k).K_us = m/P.L*(b/Cf - a/Cr); %#ok<SAGROW>
    LAT(k).Gr = Gl(1,1);  LAT(k).Gay = Gl(2,1);  LAT(k).Gbeta = Gl(3,1);
    [LAT(k).num, LAT(k).den] = tfdata(LAT(k).Gr, 'v');
    % errori rispetto a un riferimento DRITTO (angoli piccoli): e_psi' = r, e_y'' = a_y
    LAT(k).Gepsi = LAT(k).Gr / s;          % errore di direzione / delta   [rad/rad]
    LAT(k).Gey   = LAT(k).Gay / s^2;       % errore laterale / delta       [m/rad]
    tau = X.sigma/v;                                   % rilassamento: stati [vy r alpha_f alpha_r]
    A4 = [ 0, -v, Cf/m, Cr/m;
           0,  0, a*Cf/Iz, -b*Cr/Iz;
          -1/(v*tau), -a/(v*tau), -1/tau, 0;
          -1/(v*tau),  b/(v*tau), 0, -1/tau ];
    LAT(k).Gr_rel = tf(ss(A4, [0;0;1/tau;0], [0 1 0 0], 0, 'InputDelay', X.Td_steer));
    Kdc = dcgain(LAT(k).Gr);                           % stesso guadagno statico della bicicletta
    LAT(k).Gr_dp = Kdc/(X.tau_dp*s + 1);  LAT(k).Gr_dp.InputDelay = X.L_dp;

    % ---- LONGITUDINALE: marcia secondo le soglie di salita identificate
    rpm_g = v/P.Rw*X.G*60/(2*pi);
    g = find(rpm_g(1:5) <= X.upshift, 1);  if isempty(g), g = 6; end
    G = X.G(g);  n0 = rpm_g(g);
    meq = m + X.J*(G/P.Rw)^2;
    % smorzamento: derivata della resistenza aerodinamica + del freno motore (giri ~ v)
    kd_road = 2*X.c2*v;
    kd_eb   = (G/P.Rw) * X.e1 * (G/P.Rw)*60/(2*pi);
    Froad = X.c0 + X.c2*v^2;
    % pedale a regime per tenere v0: T_net = F_road*Rw/G = f*(T_e + T_eb) - T_eb
    Te = Te_of(n0);  Teb = Teb_of(n0);
    f0 = (Froad*P.Rw/G + Teb)/(Te + Teb);
    p0 = p_of(f0);
    dp = 0.02;  dfdp = (f_of(min(p0+dp,1)) - f_of(max(p0-dp,0)))/(min(p0+dp,1) - max(p0-dp,0));
    KFp = (Te + Teb)*dfdp*G/P.Rw;                      % [N per unita' di pedale] attorno a p0
    LON(k).v = v;  LON(k).gear = g;  LON(k).rpm = n0;  LON(k).meq = meq; %#ok<SAGROW>
    LON(k).pedale0 = p0;  LON(k).dFdp = KFp;  LON(k).tau = meq/(kd_road + kd_eb);
    LON(k).GF   = 1/(meq*s + kd_road);                 % dv/dFx (forza a terra), frizione indipendente
    LON(k).Gped = KFp/(meq*s + kd_road);  LON(k).Gped.InputDelay = X.Td_ped;
    % freni (pedale 0, marcia inserita: c'e' anche il freno motore come smorzamento)
    LON(k).Gbrk = -BRK.KF_per_cmd/(meq*s + kd_road + kd_eb) * 1/(BRK.tau*s + 1);
    LON(k).Gbrk.InputDelay = BRK.Td;
    LON(k).brk_decel_per_cmd = BRK.KF_per_cmd/meq;    % [m/s^2 per unita' di comando] a regime di forza
end

%% 3) STAMPA -----------------------------------------------------------------
fprintf('\n===== FRENI (attuatore, per Simulink) =====\n');
fprintf('pressione = %.3g * comando, satura a %.3g Pa (comando %.0f)\n', BRK.gain, BRK.psat, BRK.cmd_sat);
fprintf('velocita'' max pressione: salita %.3g Pa/s, discesa %.3g Pa/s  (Rate Limiter)\n', BRK.slew_up, BRK.slew_down);
fprintf('K_f = %.3e, K_r = %.3e Nm/Pa per ruota, bias ant. %.2f -> %.1f N di forza per unita'' di comando\n', ...
    BRK.Kf, BRK.Kr, BRK.bias_f, BRK.KF_per_cmd);
fprintf('(comando massimo utile %.0f -> %.0f N, cioe'' %.1f m/s^2)\n', BRK.cmd_sat, BRK.KF_per_cmd*BRK.cmd_sat, BRK.KF_per_cmd*BRK.cmd_sat/P.m);
for k = 1:numel(V0)
    fprintf('\n===== v0 = %d m/s =====\n', V0(k));
    fprintf('LATERALE: C_f = %.3g, C_r = %.3g N/rad, K_us = %.2e; poli %s\n', LAT(k).Cf, LAT(k).Cr, LAT(k).K_us, ...
        mat2str(round(pole(LAT(k).Gr).', 2)));
    fprintf('r/delta:');  LAT(k).Gr   %#ok<NOPTS>
    fprintf('ay/delta:');  LAT(k).Gay  %#ok<NOPTS>
    fprintf('e_y/delta (errore laterale):');  LAT(k).Gey  %#ok<NOPTS>
    fprintf('LONGITUDINALE: %da marcia, %.0f rpm, m_eq %.0f kg, pedale di regime %.2f, dF/dpedale %.0f N, tau %.1f s\n', ...
        LON(k).gear, LON(k).rpm, LON(k).meq, LON(k).pedale0, LON(k).dFdp, LON(k).tau);
    fprintf('v/pedale:');  LON(k).Gped  %#ok<NOPTS>
    fprintf('v/freno:');   LON(k).Gbrk  %#ok<NOPTS>
end

%% 4) GRAFICI ----------------------------------------------------------------
lg = arrayfun(@(x) sprintf('%d m/s', x), V0, 'UniformOutput', false);
figure('Name', 'Sterzo -> imbardata');
subplot(1,2,1); hold on; for k = 1:numel(V0), bode(LAT(k).Gr); end; grid on; legend(lg); title('r/\delta');
subplot(1,2,2); hold on; for k = 1:numel(V0), step(0.01*LAT(k).Gr, 0.6); end; grid on; title('r per \delta = 0.01 rad');
figure('Name', 'Sterzo -> accelerazione ed errore laterale');
subplot(1,2,1); hold on; for k = 1:numel(V0), bode(LAT(k).Gay); end; grid on; legend(lg); title('a_y/\delta');
subplot(1,2,2); hold on; for k = 1:numel(V0), step(0.001*LAT(k).Gey, 1.5); end; grid on; title('e_y per \delta = 0.001 rad');
figure('Name', 'Pedale e freni -> velocita''');
subplot(1,2,1); hold on; for k = 1:numel(V0), step(0.05*LON(k).Gped, 5); end; grid on; legend(lg);
title('\Deltav per pedale +0.05');
subplot(1,2,2); hold on; for k = 1:numel(V0), step(10*LON(k).Gbrk, 5); end; grid on;
title('\Deltav per comando freni +10');

save('fdt_polimove.mat', 'LAT', 'LON', 'BRK', 'X', 'P', 'V0');
