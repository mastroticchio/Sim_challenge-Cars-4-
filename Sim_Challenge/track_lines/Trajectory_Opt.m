%% ========================================================================
% OTTIMIZZAZIONE TRAIETTORIA CON SPLINE ANALITICHE E QP (CIRCUITO MONZA)
% ========================================================================
clear; clc; close all;

%% 1. Importazione e Preparazione Dati
fprintf('1. Importazione dati...\n');
D_innertrackline = readmatrix("inner_track_line.csv");
D_outertrackline = readmatrix("outer_track_line.csv");

% Estraiamo le righe 2 (X) e 3 (Y) e FORZIAMO subito a vettori colonna (N x 1)
X_inner = D_innertrackline(2, :)';
Y_inner = D_innertrackline(3, :)';
X_outer = D_outertrackline(2, :)';
Y_outer = D_outertrackline(3, :)';

% Calcolo Centerline grezza
Xc_raw = (X_inner + X_outer) / 2;
Yc_raw = (Y_inner + Y_outer) / 2;

% Numero di punti
N = length(Xc_raw);

% Passo ds costante dalla colonna S (Riga 1)
ds = mean(diff(D_innertrackline(1, :)));

%% 2. Pulizia e Derivate Analitiche tramite Spline Cubica
fprintf('2. Calcolo derivate analitiche tramite Spline...\n');
% Applichiamo un leggerissimo filtro per smussare gli scalini di campionamento del CSV
% in modo che la spline non "serpeggi" per rincorrere il rumore
Xc_smooth = smoothdata(Xc_raw, 'sgolay', 15);
Yc_smooth = smoothdata(Yc_raw, 'sgolay', 15);

% Creiamo un parametro arbitrario t per parametrizzare la curva
t = (1:N)';

% Interpoliamo con spline cubica
ppX = spline(t, Xc_smooth);
ppY = spline(t, Yc_smooth);

% Estrazione dei coefficienti [a, b, c, d] del polinomio: a*t^3 + b*t^2 + c*t + d
cX = ppX.coefs;
cY = ppY.coefs;

% Inizializziamo i vettori delle derivate
dX_c = zeros(N, 1); ddX_c = zeros(N, 1);
dY_c = zeros(N, 1); ddY_c = zeros(N, 1);

% La derivata prima in t=0 è il coeff 'c'. La derivata seconda è '2*b'
dX_c(1:N-1)  = cX(:, 3);
ddX_c(1:N-1) = 2 * cX(:, 2);
dY_c(1:N-1)  = cY(:, 3);
ddY_c(1:N-1) = 2 * cY(:, 2);

% L'ultimo punto (N) va valutato alla fine dell'ultimo intervallo (t=1)
dX_c(N)  = 3 * cX(end, 1) + 2 * cX(end, 2) + cX(end, 3);
ddX_c(N) = 6 * cX(end, 1) + 2 * cX(end, 2);
dY_c(N)  = 3 * cY(end, 1) + 2 * cY(end, 2) + cY(end, 3);
ddY_c(N) = 6 * cY(end, 1) + 2 * cY(end, 2);

% Curvatura della Centerline (ora liscia e perfetta)
num_c = (dX_c .* ddY_c) - (dY_c .* ddX_c);
den_c = (dX_c.^2 + dY_c.^2).^(3/2);
K_c = num_c ./ den_c;

% Vettori tangenti e normali (Rotazione +90° antioraria, verso l'ESTERNO)
norma = sqrt(dX_c.^2 + dY_c.^2);
tx = dX_c ./ norma;
ty = dY_c ./ norma;
nx = -ty; 
ny =  tx;

%% 3. Setup del Quadratic Program (QP)
fprintf('3. Configurazione Problema di Ottimizzazione...\n');
c_val = 1 / (ds^2);

% Costruzione matrice Tridiagonale
diag_minus1 =  c_val * ones(N, 1);
diag_main   = -2*c_val * ones(N, 1);
diag_plus1  =  c_val * ones(N, 1);
D = spdiags([diag_minus1, diag_main, diag_plus1], [-1, 0, 1], N, N);

% Chiusura ad anello per Monza
D(1, N) = c_val;
D(N, 1) = c_val;

% Matrice Hessiana (con Regolarizzazione di Tikhonov)
w_reg = 1e-5; 
H = 2 * (D' * D) + w_reg * speye(N);

% Gradiente
f = 2 * (D' * K_c);

% Limiti (Bounds)
w_car = 1.6059;            
safety_margin = 0.2;    
clearance = (w_car / 2) + safety_margin;

d_inner = sqrt((Xc_raw - X_inner).^2 + (Yc_raw - Y_inner).^2);
d_outer = sqrt((Xc_raw - X_outer).^2 + (Yc_raw - Y_outer).^2);

% I bound rispetto al vettore normale (che punta all'esterno)
lb = -d_inner + clearance;
ub =  d_outer - clearance;

%% 4. Esecuzione Ottimizzazione
options = optimoptions('quadprog', ...
    'Algorithm', 'interior-point-convex', ...
    'StepTolerance', 1e-12, ...        
    'OptimalityTolerance', 1e-12, ...  
    'MaxIterations', 3000, ...
    'Display', 'iter');

fprintf('4. Avvio Solutore QP...\n');
alpha_opt = quadprog(H, f, [], [], [], [], lb, ub, [], options);

if isempty(alpha_opt)
    error('Ottimizzazione fallita!');
end

% Generazione Coordinate X, Y finali
X_opt = Xc_raw + alpha_opt .* nx;
Y_opt = Yc_raw + alpha_opt .* ny;

%% 5. Calcolo del Profilo di Velocità Filtrato G(s) = 1/(tau*s + 1)
fprintf('5. Calcolo e filtraggio del profilo di velocità...\n');

% Calcoliamo la curvatura della nuova linea direttamente dalla matrice (ZERO rumore)
K_opt = K_c + (D * alpha_opt);

% Eliminazione infinitesimi sui rettilinei per non far esplodere la formula
K_opt(abs(K_opt) < 1e-4) = 1e-4;

% Parametri veicolo
g = 9.81;               
mu = 1.2;               
V_max = 80.0; % [m/s] -> 288 km/h

% Velocità geometrica ideale
V_opt_raw = sqrt((g * mu) ./ abs(K_opt));
V_opt_raw(V_opt_raw > V_max) = V_max;

% Applicazione filtro discreto 1/(tau*s + 1)
tau = 2.0; % Costante spaziale [m]. Aumentala (es. 10.0 o 20.0) se vuoi curve più dolci
alpha_filt = ds / (tau + ds);

V_filtered = zeros(N, 1);
V_filtered(1) = V_opt_raw(1);

for k = 2:N
    V_filtered(k) = alpha_filt * V_opt_raw(k) + (1 - alpha_filt) * V_filtered(k - 1);
end
% Chiusura anello: facciamo un secondo giro a vuoto per annullare lo scalino al traguardo
for k = 2:N
    V_filtered(k) = alpha_filt * V_opt_raw(k) + (1 - alpha_filt) * V_filtered(k - 1);
end
%% 6. Esportazione File per ROS 2 (3 x N)
fprintf('6. Esportazione CSV...\n');
Export_Matrix = [X_opt, Y_opt, V_filtered]'; % L'apice traspone da [N x 3] a [3 x N]
writematrix(Export_Matrix, 'optimized_racing_line_with_velocity.csv');
fprintf('Successo! File salvato.\n');

%% 7. Grafici di Verifica
figure('Name', 'Ottimizzazione Racing Line', 'NumberTitle', 'off', 'Position', [100 100 1200 500]);

% Plot 1: Mappa del Circuito
subplot(1, 2, 1); hold on; axis equal; grid on;
plot(X_inner, Y_inner, 'k-', 'LineWidth', 1);
plot(X_outer, Y_outer, 'k-', 'LineWidth', 1);
plot(Xc_raw, Yc_raw, 'b--', 'DisplayName', 'Centerline');
plot(X_opt, Y_opt, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Racing Line');
title('Traiettoria Ottimizzata'); xlabel('Est [m]'); ylabel('Nord [m]'); legend;

% Plot 2: Profilo di Velocità
subplot(1, 2, 2); hold on; grid on;
plot(V_opt_raw .* 3.6, 'Color', [0.8 0.8 0.8], 'DisplayName', 'V limite teorico');
plot(V_filtered .* 3.6, 'r-', 'LineWidth', 2, 'DisplayName', 'V filtrata G(s)');
title('Profilo di Velocità (km/h)'); xlabel('Punti del tracciato'); ylabel('Velocità [km/h]');
ylim([0, 310]); legend('Location', 'best');