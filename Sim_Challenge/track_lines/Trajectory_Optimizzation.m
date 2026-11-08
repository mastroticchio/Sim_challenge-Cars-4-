D_innertrackline = readmatrix("inner_track_line.csv");
D_outertrackline = readmatrix("outer_track_line.csv");
X_inner = D_innertrackline(2,:);
Y_inner = D_innertrackline(3,:);
X_outer = D_outertrackline(2,:);
Y_outer = D_outertrackline(3,:);
Xc = (X_inner + X_outer)/2;
Yc = (Y_inner + Y_outer)/2;
ds = mean(diff(D_innertrackline(1,:)));
c = 1 / (ds^2);       % La nostra costante
N = size(Y_inner,2);
Reference_OLD = readmatrix("reference_line.csv");
V_Mean_OLD = mean(Reference_OLD(3,:));

% Forziamo le coordinate originali ad essere rigorosamente vettori colonna (N x 1)
Xc = Xc(:);
Yc = Yc(:);

X_inner = X_inner(:);
Y_inner = Y_inner(:);

X_outer = X_outer(:);
Y_outer = Y_outer(:);

% Supponendo di avere i vettori colonna Xc e Yc già calcolati

% Derivate prime e componenti dei vettori tangenti
dX_c = gradient(Xc);
dY_c = gradient(Yc);
norma = sqrt(dX_c.^2 + dY_c.^2);
tx = dX_c./norma;
ty = dY_c./norma;
%Calcolo vettori normali (positivo esterno pista)
nx=-ty;
ny=tx;
% Assicuriamoci che anche i vettori normali siano lisci prima di usarli
nx = smoothdata(nx, 'sgolay', 21);
ny = smoothdata(ny, 'sgolay', 21);
% Derivate seconde
ddX_c = gradient(dX_c);
ddY_c = gradient(dY_c);

% Numeratore e Denominatore della formula
numeratore = (dX_c .* ddY_c) - (dY_c .* ddX_c);
denominatore = (dX_c.^2 + dY_c.^2).^(3/2);

% Curvatura della mezzeria (Vettore colonna Nx1)
K_c = numeratore ./ denominatore;
K_c = smoothdata(K_c, 'sgolay', 21);
% 1. Costruiamo le 3 colonne per le diagonali di D
diag_minus1 =  c * ones(N, 1);
diag_main   = -2*c * ones(N, 1);
diag_plus1  =  c * ones(N, 1);

% 2. Assembliamo la matrice sparsa D (Tridiagonale)
D = spdiags([diag_minus1, diag_main, diag_plus1], [-1, 0, 1], N, N);

% AGGIUNTA FONDAMENTALE: Chiusura ad anello per circuito di Monza
D(1, N) = c;
D(N, 1) = c;

% 3. Calcoliamo H con Regolarizzazione di Tikhonov
w_reg = 1e-5; % Elastico per evitare il blocco del solutore
H = 2 * (D' * D) + w_reg * speye(N);

% 4. Calcolo del vettore f (Gradiente)
% f rappresenta la parte lineare legata alla curvatura originaria della centerline (K_c)
% Il solutore cerca di far sì che alpha'' compensi K_c
f = 2 * (D' * K_c);

% Parametri fisici dell'auto
w_car = 1.6059;            % Larghezza auto in metri
safety_margin = 0.2;    % Margine di sicurezza dai cordoli in metri
clearance = (w_car / 2) + safety_margin;

% 1. Calcolo delle distanze fisiche dalla centerline ai bordi
d_inner = sqrt((Xc - X_inner).^2 + (Yc - Y_inner).^2);
d_outer = sqrt((Xc - X_outer).^2 + (Yc - Y_outer).^2);

% 2. Costruzione dei vettori dei limiti (N x 1)
% Assumendo che il vettore normale punti verso l'ESTERNO della curva
lb = -d_inner + clearance;
ub =  d_outer - clearance;

% Controllo di sicurezza: se la pista è troppo stretta lb potrebbe superare ub
% mandando in crash il solutore. Verifichiamo che ci sia spazio fisico:
if any(lb >= ub)
    warning('Attenzione: in alcuni punti la pista è più stretta dell''auto!');
end

%% 1. Configurazione del solutore
% Diciamo a MATLAB quale algoritmo usare. 
% 'interior-point-convex' è un algoritmo potentissimo e velocissimo 
% studiato appositamente per le matrici sparse di grandi dimensioni come la nostra H.
options = optimoptions('quadprog', ...
    'Algorithm', 'interior-point-convex', ...
    'StepTolerance', 1e-12, ...        % Cruciale per evitare stop prematuri
    'OptimalityTolerance', 1e-12, ...  % Cruciale per tollerare gradienti piccoli
    'MaxIterations', 3000, ...
    'Display', 'iter');

%% 2. Lancio dell'Ottimizzazione
% Sintassi completa: quadprog(H, f, A, b, Aeq, beq, lb, ub, x0, options)
% I parametri A, b (disequazioni) e Aeq, beq (equazioni) non ci servono e li lasciamo vuoti [].
fprintf('Avvio ottimizzazione traiettoria (solutore QP)...\n');

% quadprog restituisce il vettore ottimizzato degli alpha
alpha_opt = quadprog(H, f, [], [], [], [], lb, ub, [], options);

%% 3. Controllo del fallimento
% Se la pista è troppo stretta (es. lb > ub in un punto) il solutore si arrende e restituisce vuoto
if isempty(alpha_opt)
    error('Ottimizzazione fallita! Il solutore non ha trovato una via d''uscita. Controlla lb e ub.');
else
    fprintf('Ottimizzazione conclusa con successo!\n');
end

%% 4. Riconversione in Coordinate Globali (X, Y)
% Usiamo la formula vettoriale element-wise vista prima per tornare al mondo reale
% Assicurati di avere in memoria Xc, Yc, nx, ny
X_opt = Xc + alpha_opt .* nx;
Y_opt = Yc + alpha_opt .* ny;

%% 5. Visualizzazione Grafica
figure('Name', 'Ottimizzazione Racing Line', 'NumberTitle', 'off');
hold on; 
axis equal; % FONDAMENTALE: rende le proporzioni X e Y uguali, altrimenti la curva sembra schiacciata
grid on;

% Disegniamo la pista
plot(X_inner, Y_inner, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Inner Bound');
plot(X_outer, Y_outer, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Outer Bound');

% Disegniamo la linea centrale originale
plot(Xc, Yc, 'b--', 'DisplayName', 'Centerline Originale');

% Disegniamo il capolavoro: la traiettoria ottimizzata
plot(X_opt, Y_opt, 'r-', 'LineWidth', 2, 'DisplayName', 'Racing Line (Ottimizzata)');

legend('Location', 'best');
title('Risultato Ottimizzazione Curvatura Minima');
xlabel('Est [m]');
ylabel('Nord [m]');

%% 6. Calcolo del Profilo di Velocità Ottimo

% Ricalcoliamo la curvatura ESATTA della nuova traiettoria ottimizzata
dX_opt = gradient(X_opt);
dY_opt = gradient(Y_opt);
ddX_opt = gradient(dX_opt);
ddY_opt = gradient(dY_opt);

num_opt = (dX_opt .* ddY_opt) - (dY_opt .* ddX_opt);
den_opt = (dX_opt.^2 + dY_opt.^2).^(3/2);
K_opt = num_opt ./ den_opt; % Curvatura della racing line
K_opt_clean = smoothdata(K_opt, 'sgolay', 50);
% Parametri fisici (adattali al veicolo PoliMove)
g = 9.81;               % Accelerazione di gravità [m/s^2]
mu = 1;               % Coefficiente di aderenza (stima per gomme racing)
V_max = 80.0;           % Velocità massima fisica dell'auto in rettilineo [m/s] (80 m/s = 288 km/h)

% Calcolo della velocità limite (Steady-state cornering)
% Usiamo abs() per evitare numeri complessi nelle curve a destra
V_opt = sqrt((g * mu) ./ abs(K_opt_clean))./2.5;

% Saturazione: tagliamo tutti i picchi infiniti nei rettilinei alla V_max dell'auto
V_opt(V_opt > V_max) = V_max;
V_Mean_NEW = mean(V_opt);
% Filtro passa basso velocità


%% 7. Esportazione nel formato richiesto (3 righe x N colonne)

% Affianchiamo i vettori colonna [N x 3]
Data_Matrix = [X_opt, Y_opt, V_opt];

% TRASPOSIZIONE: 
% Trasponiamo la matrice usando l'apice (') per farla diventare [3 x N]
% Riga 1: X
% Riga 2: Y
% Riga 3: V
Export_Data = Data_Matrix'; 

% Salvataggio finale
writematrix(Export_Data, 'optimized_racing_line_with_velocity.csv');
fprintf('Esportazione completata: "optimized_racing_line_with_velocity.csv" (3 Righe, %d Colonne)\n', N);

%% (Opzionale) Plot del Profilo di Velocità
figure('Name', 'Profilo di Velocità', 'NumberTitle', 'off');
plot(V_opt .* 3.6, 'LineWidth', 1.5); % Moltiplico per 3.6 per plottare in km/h
grid on;
title('Profilo di Velocità lungo il giro');
xlabel('Indice del punto (da 1 a N)');
ylabel('Velocità [km/h]');

