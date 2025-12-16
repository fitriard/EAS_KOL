%% ========================================================================
%  DESAIN LQG, ANALISIS KESTABILAN (LYAPUNOV) & PERFORMANSI
%  SISTEM POULTRY HOUSE
%  Fitri Ardianti - EAS Kontrol Optimal
%% ========================================================================

clear; clc; close all;

fprintf('=== DESAIN LQG SISTEM POULTRY HOUSE ===\n\n');

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
    fprintf('✓ Sistem controllable → LQR (bagian dari LQG) bisa diterapkan.\n\n');
else
    error('Sistem tidak controllable.');
end

%% ------------------------------------------------------------------------
%  BAGIAN 2 : DESAIN LQG  (LQR + KALMAN FILTER)
%% ------------------------------------------------------------------------

fprintf('=== DESAIN LQG (LQR + KALMAN FILTER) ===\n');

% ----- 2.1 LQR (controller) -----
Q = diag([100 50]);          % penalti state
R = diag([1 1 0.1]);         % penalti sinyal kontrol

[K,P_ctrl,eig_cl] = lqr(A,B,Q,R);
Acl = A - B*K;

fprintf('Gain LQR K:\n'); disp(K);
fprintf('Pole closed-loop LQR (eig(A-BK)):\n'); disp(eig_cl);

% ----- 2.2 Kalman filter (state estimator) -----
% W : kovariansi noise proses, V : kovariansi noise pengukuran
W = diag([0.05 0.05]);       % misal gangguan model relatif kecil
V = diag([0.5 0.5]);         % noise sensor agak lebih besar

% LQE (continuous-time Kalman filter)
[L,P_est,eig_est] = lqe(A, eye(2), C, W, V);
% (MATLAB: L matriks gain estimator, P_est kovariansi error steady-state)

fprintf('Gain Kalman L:\n'); disp(L);
fprintf('Pole estimator (eig(A-LC)):\n'); disp(eig(A-L*C));

% ----- 2.3 Sistem LQG gabungan (controller + estimator) -----
% augmented system [x; xhat]
A_lqg = [A-B*K,       B*K;
         zeros(2),  A-L*C];    % 4x4
Bd_lqg = [Bd;
          Bd];                 % gangguan masuk ke plant dan estimator

fprintf('Pole sistem LQG (augmented):\n'); disp(eig(A_lqg));
fprintf('\n');

%% ------------------------------------------------------------------------
%  BAGIAN 3 : ANALISIS KESTABILAN (LYAPUNOV)
%% ------------------------------------------------------------------------

fprintf('=== ANALISIS KESTABILAN LQG (LYAPUNOV) ===\n');

% --- 3.1 Controller LQR ---
eig_P = eig(P_ctrl);
fprintf('\nEigen P (controller Riccati):\n'); disp(eig_P);
if all(eig_P>0)
    fprintf('✓ P_ctrl positive definite → V_c(x)=x''P_cx > 0.\n');
end

M_c = Q + K'*R*K;           % untuk LQR: Vdot = -x''M_c x
eig_Mc = eig(M_c);
fprintf('Eigen M_c = Q + K''RK:\n'); disp(eig_Mc);
if all(eig_Mc>0)
    fprintf('✓ M_c positive definite → Vdot_c < 0 (controller asimtotik stabil).\n');
end

% --- 3.2 Estimator (Kalman filter) ---
eig_Pe = eig(P_est);
fprintf('\nEigen P_est (error covariance estimator):\n'); disp(eig_Pe);
if all(eig_Pe>0)
    fprintf('✓ P_est positive definite → kandidat fungsi Lyapunov untuk error estimasi.\n');
end

% cek persamaan Riccati estimator
ARE_est = A*P_est + P_est*A' - P_est*C'/V*C*P_est + W;
fprintf('Norm residual ARE estimator = %.3e\n',norm(ARE_est,'fro'));

% --- 3.3 Sistem LQG gabungan (plant + estimator) ---
% solusi Lyapunov untuk A_lqg
Q_lyap = eye(4);                      % sembarang positive definite
P_lqg  = lyap(A_lqg', Q_lyap);        % A'P + P A = -Q
eig_Plqg = eig(P_lqg);

fprintf('\nEigen P_lqg (solusi Lyapunov augmented):\n'); disp(eig_Plqg);
if all(eig_Plqg>0) && all(real(eig(A_lqg))<0)
    fprintf('✓ P_lqg>0 & semua pole(A_lqg) di LHP → sistem LQG ASIMTOTIK STABIL.\n\n');
else
    fprintf('⚠ Periksa kembali, stabilitas belum terjamin.\n\n');
end

%% ------------------------------------------------------------------------
%  BAGIAN 4 : SIMULASI CLOSED-LOOP LQG
%% ------------------------------------------------------------------------

fprintf('=== SIMULASI CLOSED-LOOP LQG ===\n');

dt = 0.1;
t  = 0:dt:500;
N  = numel(t);

x      = zeros(2,N);   % state nyata
xhat   = zeros(2,N);   % state estimasi
u      = zeros(3,N);   % sinyal kontrol
y_meas = zeros(2,N);   % pengukuran (dengan noise)

x(:,1)    = [5; 10];   % deviasi awal plant
xhat(:,1) = [0; 0];    % tebakan awal estimator

% gangguan seperti di kasus LQR
d1 = 3*ones(1,N);      % suhu luar +3°C
d2 = 5*ones(1,N);      % RH luar +5%%
d2(t>250) = -2;        % step ke -2%%

% noise pengukuran v(k) ~ N(0,V)
Rv = chol(V,'lower');
v  = Rv*randn(2,N);    % satu realisasi noise

for k = 1:N-1
    d = [d1(k); d2(k)];
    
    % pengukuran dengan noise
    y_meas(:,k) = C*x(:,k) + v(:,k);

    % kontrol menggunakan state estimasi
    u(:,k) = -K*xhat(:,k);
    
    % dinamika plant
    x_dot = A*x(:,k) + B*u(:,k) + Bd*d;
    x(:,k+1) = x(:,k) + dt*x_dot;
    
    % dinamika estimator (Kalman)
    xhat_dot = A*xhat(:,k) + B*u(:,k) + L*(y_meas(:,k) - C*xhat(:,k));
    xhat(:,k+1) = xhat(:,k) + dt*xhat_dot;
end
y_meas(:,N) = C*x(:,N) + v(:,N);
u(:,N)      = u(:,N-1);

fprintf('Simulasi selesai.\n\n');

%% ------------------------------------------------------------------------
%  BAGIAN 5 : ANALISIS PERFORMANSI LQG
%% ------------------------------------------------------------------------

fprintf('=== ANALISIS PERFORMANSI LQG ===\n');

% 1) Settling time (2%% band dari nilai awal)
idx_ts1 = find(abs(x(1,:)) < 0.02*abs(x(1,1)),1);
idx_ts2 = find(abs(x(2,:)) < 0.02*abs(x(2,1)),1);

if isempty(idx_ts1), Ts1 = t(end); else, Ts1 = t(idx_ts1); end
if isempty(idx_ts2), Ts2 = t(end); else, Ts2 = t(idx_ts2); end

fprintf('Settling time (2%%%% band):\n');
fprintf('  x1 (suhu)      : %.2f s\n',Ts1);
fprintf('  x2 (kelembaban): %.2f s\n\n',Ts2);

% 2) Overshoot (terhadap nilai awal)
OS1 = (max(x(1,:)) - x(1,1))/abs(x(1,1))*100;
OS2 = (max(x(2,:)) - x(2,1))/abs(x(2,1))*100;

fprintf('Overshoot:\n');
fprintf('  x1 (suhu)      : %.2f %%%%\n',OS1);
fprintf('  x2 (kelembaban): %.2f %%%%\n\n',OS2);

% 3) Steady-state error
Ess1 = x(1,end);
Ess2 = x(2,end);
fprintf('Steady-state error:\n');
fprintf('  x1 (suhu)      : %.4f °C\n',Ess1);
fprintf('  x2 (kelembaban): %.4f %%%%\n\n',Ess2);

% 4) Energi kontrol
E_u1 = trapz(t,u(1,:).^2);
E_u2 = trapz(t,u(2,:).^2);
E_u3 = trapz(t,u(3,:).^2);
E_tot = E_u1 + E_u2 + E_u3;

fprintf('Energi kontrol (∫u^2 dt):\n');
fprintf('  Heating (u1)   : %.2f\n',E_u1);
fprintf('  PAD (u2)       : %.2f\n',E_u2);
fprintf('  Vent (u3)      : %.2f\n',E_u3);
fprintf('  Total          : %.2f\n\n',E_tot);

% 5) Fungsi biaya LQR (pakai state nyata & kontrol LQG)
state_cost   = trapz(t, sum((x.'*Q).*x.',2));
state_cost   = state_cost(end);   % trik karena barusan double integrasi
control_cost = trapz(t, sum((u.'*R).*u.',2));
control_cost = control_cost(end);
J = state_cost + control_cost;

fprintf('Fungsi biaya (LQG, pakai Q,R yang sama):\n');
fprintf('  State cost     : %.2f\n',state_cost);
fprintf('  Control cost   : %.2f\n',control_cost);
fprintf('  Total J        : %.2f\n\n',J);

%% ------------------------------------------------------------------------
%  BAGIAN 6 : PLOT RESPON & SINYAL KONTROL (NILAI AKTUAL)
%% ------------------------------------------------------------------------

% Setpoint aktual (fisik)
T_sp  = 25;   % °C
RH_sp = 60;   % %

% Konversi deviasi -> nilai aktual (PAKAI STATE NYATA x)
T_actual  = x(1,:) + T_sp;
RH_actual = x(2,:) + RH_sp;

%% --- Plot respon suhu & kelembaban aktual ---
figure('Name','Respon Suhu dan Kelembaban Aktual - LQG');

yyaxis left
plot(t, T_actual, 'b', 'LineWidth', 2); hold on;
yline(T_sp, 'b--', 'LineWidth', 1.5);
ylabel('Suhu Internal (°C)');

yyaxis right
plot(t, RH_actual, 'r', 'LineWidth', 2);
yline(RH_sp, 'r--', 'LineWidth', 1.5);
ylabel('Kelembaban Internal (%)');

xlabel('Waktu (s)');
grid on;
title('Respon Suhu dan Kelembaban Internal Aktual (LQG)');
legend('Suhu Aktual','Setpoint Suhu',...
       'Kelembaban Aktual','Setpoint Kelembaban',...
       'Location','best');

%% --- Plot sinyal kontrol aktual ---
figure('Name','Sinyal Kontrol Aktual - LQG');

plot(t, u(1,:), 'b', 'LineWidth', 1.8); hold on;
plot(t, u(2,:), 'r', 'LineWidth', 1.8);
plot(t, u(3,:), 'g', 'LineWidth', 1.8);

yline(0,'k--');
xlabel('Waktu (s)');
ylabel('Sinyal Kontrol');
grid on;
title('Sinyal Kontrol LQG (Heating, PAD Cooling, Ventilasi)');
legend('Heating','PAD Cooling','Ventilasi','Location','best');
