%% Generazione Reference Line per Esperimento E1: Coast-Down
clear; clc; close all;

fprintf('Generazione Reference Line per Coast-Down...\n');

% 1. Parametri Spaziali
% Creiamo un rettilineo molto lungo (es. 4 km) per dare all'auto tutto lo 
% spazio necessario sia per l'accelerazione che per il lungo coasting.
L = 2000;         % Lunghezza totale del rettilineo [m]
ds = 0.5;         % Passo di campionamento spaziale [m]

% 2. Generazione Coordinate Geometria (Rettilineo perfetto)
X_ref = 0 : ds : L;   
N = length(X_ref);
Y_ref = zeros(1, N);  % Y sempre a 0 per mantenere l'auto in linea retta

% 3. Profilo di Velocità Target
% Chiediamo al controllore longitudinale di puntare a 65 m/s. 
% In questo modo il motore spingerà forte fino a superare la soglia dei 60 m/s 
% in cui il nodo ROS 2 disabiliterà i comandi (throttle=0, brake=0).
V_target = 65; % [m/s] -> 234 km/h
V_ref = V_target * ones(1, N);
V_ref(1,4*N/5:N)=0;

% 4. Costruzione Matrice da Esportare
% Il nodo C++ si aspetta [3 righe x N colonne]
% Riga 1: X
% Riga 2: Y
% Riga 3: V
Export_Data = [X_ref; Y_ref; V_ref];

% 5. Salvataggio File
filename = 'coastdown_reference_line.csv';
writematrix(Export_Data, filename);

fprintf('File "%s" generato con successo!\n', filename);
fprintf('Dimensioni CSV: 3 Righe x %d Colonne\n', N);

%% 6. Plot di Verifica
figure('Name', 'Reference Line - E1 Coast Down', 'NumberTitle', 'off');

% Plot Traiettoria
subplot(2,1,1);
plot(X_ref, Y_ref, 'b', 'LineWidth', 2);
title('Traiettoria (Input per Controllore Laterale)');
xlabel('X [m]'); ylabel('Y [m]');
ylim([-5 5]); grid on;

% Plot Velocità
subplot(2,1,2);
plot(X_ref, V_ref, 'r', 'LineWidth', 2);
title('Velocità Target (Input per Controllore Longitudinale)');
xlabel('X [m]'); ylabel('Velocità Target [m/s]');
ylim([0 80]); grid on;