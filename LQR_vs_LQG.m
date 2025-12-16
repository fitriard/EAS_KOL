%% ========================================================================
%  PERBANDINGAN LQR vs LQG - SISTEM POULTRY HOUSE
%  Fitri Ardianti - EAS Kontrol Optimal
%% ========================================================================

clear; clc; close all;

fprintf('=== PERBANDINGAN LQR vs LQG - SISTEM POULTRY HOUSE ===\n\n');

%% ------------------------------------------------------------------------
%  BAGIAN 1 : MODEL STATE-SPACE
%% ------------------------------------------------------------------------

Rcb   = 0.015;     % thermal resistance (K/W)
Kg    = 50;        % gain heating
a     = 1000;      % konstanta a
vb    = 500;       % volume udara 1
vc    = 600;       % volume udara 2
Rcv   = 0.01;      % coupling suhu -> kelembaban
x1e   = 25;        % suhu ekuilibrium (°C)
u3_op = 100;       % ventilasi nominal

A = [-(Rcb + Kg/a + u3_op/vb),      0;
      Rcv*(0.52*x1e - 6.46),   -u3_op/vc];

B  = [1/a,   -1/a,  -x1e/vb;
      0,     1/vb, -x1e/vc];

Bd = [Kg/a + u3_op/vb,      0;
      0,               u3_op/vc];

C = eye(2);
D = zeros(2,3);

fprintf('Matriks A:\n'); disp(A);
fprintf('Matriks B:\n'); disp(B);

% cek controllability
Co = ctrb(A,B);
if rank(Co)==size(A,1)
    fprintf('✓ Sistem controllable.\n\n');
else
    error('Sistem tidak controllable.');
end

%% ------------------------------------------------------------------------
%  BAGIAN 2 : DESAIN LQR
%% ------------------------------------------------------------------------

fprintf('=== DESAIN LQR ===\n');

Q = diag([100 50]);          % penalti state
R = diag([1 1 0.1]);         % penalti sinyal kontrol

[K_lqr, P_lqr, eig_lqr] = lqr(A,B,Q,R);
Acl_lqr = A - B*K_lqr;

fprintf('Gain LQR K:\n'); disp(K_lqr);
fprintf('Pole closed-loop LQR:\n'); disp(eig_lqr);

%% ------------------------------------------------------------------------
%  BAGIAN 3 : DESAIN LQG
%% ------------------------------------------------------------------------

fprintf('\n=== DESAIN LQG ===\n');

% Kalman filter
W = diag([0.05 0.05]);       % kovariansi noise proses
V = diag([0.5 0.5]);         % kovariansi noise pengukuran

[L, P_est, eig_est] = lqe(A, eye(2), C, W, V);

fprintf('Gain Kalman L:\n'); disp(L);
fprintf('Pole estimator:\n'); disp(eig(A-L*C));

%% ------------------------------------------------------------------------
%  BAGIAN 4 : SIMULASI LQR
%% ------------------------------------------------------------------------

fprintf('\n=== SIMULASI LQR ===\n');

dt = 0.1;
t  = 0:dt:500;
N  = numel(t);

% Inisialisasi LQR
x_lqr = zeros(2,N);
u_lqr = zeros(3,N);
x_lqr(:,1) = [5; 10];   % deviasi awal

% gangguan
d1 = 3*ones(1,N);
d2 = 5*ones(1,N);
d2(t>250) = -2;

for k = 1:N-1
    d = [d1(k); d2(k)];
    
    % kontrol LQR (full state feedback)
    u_lqr(:,k) = -K_lqr*x_lqr(:,k);
    
    % dinamika plant
    x_dot = A*x_lqr(:,k) + B*u_lqr(:,k) + Bd*d;
    x_lqr(:,k+1) = x_lqr(:,k) + dt*x_dot;
end
u_lqr(:,N) = u_lqr(:,N-1);

fprintf('Simulasi LQR selesai.\n');

%% ------------------------------------------------------------------------
%  BAGIAN 5 : SIMULASI LQG
%% ------------------------------------------------------------------------

fprintf('=== SIMULASI LQG ===\n');

x_lqg      = zeros(2,N);   % state nyata
xhat_lqg   = zeros(2,N);   % state estimasi
u_lqg      = zeros(3,N);   % sinyal kontrol
y_meas     = zeros(2,N);   % pengukuran dengan noise

x_lqg(:,1)    = [5; 10];   % deviasi awal plant
xhat_lqg(:,1) = [0; 0];    % tebakan awal estimator

% noise pengukuran
Rv = chol(V,'lower');
v  = Rv*randn(2,N);

for k = 1:N-1
    d = [d1(k); d2(k)];
    
    % pengukuran dengan noise
    y_meas(:,k) = C*x_lqg(:,k) + v(:,k);

    % kontrol menggunakan state estimasi
    u_lqg(:,k) = -K_lqr*xhat_lqg(:,k);
    
    % dinamika plant
    x_dot = A*x_lqg(:,k) + B*u_lqg(:,k) + Bd*d;
    x_lqg(:,k+1) = x_lqg(:,k) + dt*x_dot;
    
    % dinamika estimator (Kalman)
    xhat_dot = A*xhat_lqg(:,k) + B*u_lqg(:,k) + L*(y_meas(:,k) - C*xhat_lqg(:,k));
    xhat_lqg(:,k+1) = xhat_lqg(:,k) + dt*xhat_dot;
end
y_meas(:,N) = C*x_lqg(:,N) + v(:,N);
u_lqg(:,N)  = u_lqg(:,N-1);

fprintf('Simulasi LQG selesai.\n\n');

%% ------------------------------------------------------------------------
%  BAGIAN 6 : ANALISIS PERFORMANSI
%% ------------------------------------------------------------------------

fprintf('=== ANALISIS PERFORMANSI ===\n\n');

% --- LQR ---
idx_ts1_lqr = find(abs(x_lqr(1,:)) < 0.02*abs(x_lqr(1,1)),1);
if isempty(idx_ts1_lqr), Ts_lqr = t(end); else, Ts_lqr = t(idx_ts1_lqr); end

Ess1_lqr = x_lqr(1,end);
Ess2_lqr = x_lqr(2,end);

state_cost_lqr   = trapz(t, sum((x_lqr.'*Q).*x_lqr.',2));
control_cost_lqr = trapz(t, sum((u_lqr.'*R).*u_lqr.',2));
J_lqr = state_cost_lqr + control_cost_lqr;

fprintf('PERFORMANSI LQR:\n');
fprintf('  Settling time (x1) : %.2f s\n', Ts_lqr);
fprintf('  Ess (x1)           : %.4f °C\n', Ess1_lqr);
fprintf('  Ess (x2)           : %.4f %%\n', Ess2_lqr);
fprintf('  Fungsi biaya J     : %.2f\n\n', J_lqr);

% --- LQG ---
idx_ts1_lqg = find(abs(x_lqg(1,:)) < 0.02*abs(x_lqg(1,1)),1);
if isempty(idx_ts1_lqg), Ts_lqg = t(end); else, Ts_lqg = t(idx_ts1_lqg); end

Ess1_lqg = x_lqg(1,end);
Ess2_lqg = x_lqg(2,end);

state_cost_lqg   = trapz(t, sum((x_lqg.'*Q).*x_lqg.',2));
control_cost_lqg = trapz(t, sum((u_lqg.'*R).*u_lqg.',2));
J_lqg = state_cost_lqg + control_cost_lqg;

fprintf('PERFORMANSI LQG:\n');
fprintf('  Settling time (x1) : %.2f s\n', Ts_lqg);
fprintf('  Ess (x1)           : %.4f °C\n', Ess1_lqg);
fprintf('  Ess (x2)           : %.4f %%\n', Ess2_lqg);
fprintf('  Fungsi biaya J     : %.2f\n\n', J_lqg);

%% ------------------------------------------------------------------------
%  BAGIAN 7 : KONVERSI KE NILAI AKTUAL
%% ------------------------------------------------------------------------

% Setpoint aktual (fisik)
T_sp  = 25;   % °C
RH_sp = 60;   % %

% Konversi deviasi -> nilai aktual
T_lqr  = x_lqr(1,:) + T_sp;     % suhu aktual LQR
RH_lqr = x_lqr(2,:) + RH_sp;    % kelembaban aktual LQR

T_lqg  = x_lqg(1,:) + T_sp;     % suhu aktual LQG
RH_lqg = x_lqg(2,:) + RH_sp;    % kelembaban aktual LQG

%% ------------------------------------------------------------------------
%  BAGIAN 8 : PLOT PERBANDINGAN LQR vs LQG
%% ------------------------------------------------------------------------

fig = figure('Name','Perbandingan LQR vs LQG','Position',[50 50 1400 900]);

% Warna
color_lqr = [0 0.4470 0.7410];      % biru
color_lqg = [0.8500 0.3250 0.0980]; % merah

%% --- Subplot 1: Suhu Aktual (x1) ---
subplot(3,3,1);
plot(t, T_lqr, '-', 'Color', color_lqr, 'LineWidth', 2); hold on;
plot(t, T_lqg, '--', 'Color', color_lqg, 'LineWidth', 2);
yline(T_sp, 'k:', 'LineWidth', 1.5);
grid on;
xlabel('Waktu (s)');
ylabel('Suhu (°C)');
title('x_1 - Suhu Internal : LQR vs LQG', 'FontWeight', 'bold');
legend('LQR', 'LQG', 'Setpoint', 'Location', 'best');
xlim([0 500]);

%% --- Subplot 2: Kelembaban Aktual (x2) ---
subplot(3,3,2);
plot(t, RH_lqg, '--', 'Color', color_lqg, 'LineWidth', 2); hold on;
plot(t, RH_lqr, '-', 'Color', color_lqr, 'LineWidth', 2);
yline(RH_sp, 'k:', 'LineWidth', 1.5);
grid on;
xlabel('Waktu (s)');
ylabel('Kelembaban (%)');
title('x_2 - Kelembaban Internal : LQR vs LQG', 'FontWeight', 'bold');
legend('LQG', 'LQR', 'Setpoint', 'Location', 'best');
xlim([0 500]);

%% --- Subplot 3: Ringkasan Performansi (TEXT) ---
subplot(3,3,3);
axis off;
text_str = {
    '\bf{Ringkasan Performansi}', '', ...
    sprintf('Ts LQR (x1) = %.2f s', Ts_lqr), ...
    sprintf('Ts LQG (x1) = %.2f s', Ts_lqg), '', ...
    sprintf('J_{LQR}     = %.2f', J_lqr), ...
    sprintf('J_{LQG}     = %.2f', J_lqg)
};
text(0.1, 0.5, text_str, 'FontSize', 11, 'VerticalAlignment', 'middle');

%% --- Subplot 4: Sinyal Kontrol u1 (Heating) ---
subplot(3,3,4);
plot(t, u_lqr(1,:), '-', 'Color', color_lqr, 'LineWidth', 2); hold on;
plot(t, u_lqg(1,:), '--', 'Color', color_lqg, 'LineWidth', 2);
yline(0, 'k--', 'LineWidth', 1);
grid on;
xlabel('Waktu (s)');
ylabel('u_1 (Heating)');
title('Sinyal Kontrol u_1 : LQR vs LQG', 'FontWeight', 'bold');
legend('LQR', 'LQG', 'Location', 'best');
xlim([0 500]);

%% --- Subplot 5: Sinyal Kontrol u2 (PAD) ---
subplot(3,3,5);
plot(t, u_lqr(2,:), '-', 'Color', color_lqr, 'LineWidth', 2); hold on;
plot(t, u_lqg(2,:), '--', 'Color', color_lqg, 'LineWidth', 2);
yline(0, 'k--', 'LineWidth', 1);
grid on;
xlabel('Waktu (s)');
ylabel('u_2 (PAD)');
title('Sinyal Kontrol u_2 : LQR vs LQG', 'FontWeight', 'bold');
legend('LQR', 'LQG', 'Location', 'best');
xlim([0 500]);

%% --- Subplot 6: Kosong atau informasi tambahan ---
subplot(3,3,6);
axis off;

%% --- Subplot 7: Sinyal Kontrol u3 (Ventilasi) ---
subplot(3,3,7);
plot(t, u_lqr(3,:), '-', 'Color', color_lqr, 'LineWidth', 2); hold on;
plot(t, u_lqg(3,:), '--', 'Color', color_lqg, 'LineWidth', 2);
yline(0, 'k--', 'LineWidth', 1);
grid on;
xlabel('Waktu (s)');
ylabel('u_3 (Ventilasi)');
title('Sinyal Kontrol u_3 : LQR vs LQG', 'FontWeight', 'bold');
legend('LQR', 'LQG', 'Location', 'best');
xlim([0 500]);
ylim([0 20]);          % MAKS Y-AXIS = 20
yticks(0:5:20);        % interval 5

%% --- Subplot 8 & 9: Kosong ---
subplot(3,3,8);
axis off;

subplot(3,3,9);
axis off;

%% Judul keseluruhan
sgtitle('Perbandingan LQR vs LQG - Sistem Poultry House', ...
        'FontSize', 16, 'FontWeight', 'bold');

fprintf('=== PLOT PERBANDINGAN SELESAI ===\n');