%% PARAMETER FITTING DA DATI REALI (ROS 2 SIMULATOR)
clear; clc; close all;

%% 1. Importazione Dati dal Simulatore
% Assicurati che il file abbia il tempo nella colonna 1 e la velocità (m/s) nella colonna 2
fprintf('1. Caricamento dati di telemetria...\n');
data = readmatrix('telemetry_coastdown.csv');

t_meas = data(:, 1);
v_meas = data(:, 2);

% Normalizziamo il tempo per farlo partire esattamente da t = 0
t_meas = t_meas - t_meas(1);
v0 = v_meas(1); % Velocità iniziale dell'esperimento (es. 60 m/s)

%% 2. Parametri Fisici Noti (da impostare in base alla tua vettura)
m = 720;                 % Massa veicolo [kg]
g = 9.81;                % Gravità [m/s^2]
rho = 1.225;             % Densità aria [kg/m^3]
CdA = 1.194;              % Prodotto nominale C_d * Area [m^2]
ClA = 0.8;               % Prodotto C_l * Area (Downforce) [m^2]

%% 3. Definizione del Modello Matematico (ODE)
% dv/dt = -(1/m) * [ F_aero + crr*(mg + F_down) + c0 ]
coastdown_ode = @(t, v, crr, c0) -(1/m) * ...
    (0.5 * rho * CdA * v.^2 + crr * (m*g + 0.5 * rho * ClA * v.^2) + c0);

%% 4. Setup Algoritmo di Fitting
fprintf('2. Avvio ottimizzazione per trovare crr e c0...\n');

% Funzione di Costo
% Inserendo 't_meas' in ode45, MATLAB calcola la velocità simulata esattamente 
% negli stessi istanti di tempo dei tuoi sensori ROS 2, permettendo un confronto perfetto.
cost_func = @(params) sum((v_meas - run_sim(t_meas, v0, coastdown_ode, params(1), params(2))).^2);

% Ipotesi di partenza per aiutare l'algoritmo (valori tipici)
initial_guess = [0.015, 50]; % [crr_iniziale, c0_iniziale]

% Opzioni solutore
options_opt = optimset('Display', 'iter', 'TolFun', 1e-6, 'TolX', 1e-6);

% Esecuzione ricerca
best_params = fminsearch(cost_func, initial_guess, options_opt);

crr_fit = best_params(1);
c0_fit  = best_params(2);

fprintf('\n--- RISULTATI DEL SIMULATORE ROS 2 ---\n');
fprintf('Coefficiente di rotolamento (crr) stimato : %.5f\n', crr_fit);
fprintf('Attrito meccanico (c0) stimato          : %.2f N\n', c0_fit);

%% 5. Plot di Verifica (Cruciale per confermare il risultato)
% Generiamo la curva fittata usando i parametri appena trovati
v_fit = run_sim(t_meas, v0, coastdown_ode, crr_fit, c0_fit);

figure('Name', 'Verifica Fitting', 'NumberTitle', 'off');
hold on; grid on;
plot(t_meas, v_meas .* 3.6, 'k.', 'MarkerSize', 10, 'DisplayName', 'Dati Simulatore (ROS 2)');
plot(t_meas, v_fit .* 3.6, 'r-', 'LineWidth', 2, 'DisplayName', 'Modello Fittato (Ottimizzato)');
title('Verifica del Modello vs Dati Reali');
xlabel('Tempo [s]');
ylabel('Velocità [km/h]');
legend('Location', 'best');

%% Funzione Helper per eseguire l'ODE in modo pulito
function v_out = run_sim(t_array, v_init, ode_func, crr, c0)
    % Risolve l'equazione differenziale e restituisce solo il vettore velocità
    [~, v_sim] = ode45(@(t,v) ode_func(t, v, crr, c0), t_array, v_init);
    v_out = v_sim;
end