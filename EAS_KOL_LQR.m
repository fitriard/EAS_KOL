%% ========================================================================
%  DESAIN LQR, ANALISIS KESTABILAN (LYAPUNOV) & PERFORMANSI
%  SISTEM POULTRY HOUSE
%  Fitri Ardianti - EAS Kontrol Optimal
%% ========================================================================

clear; clc; close all;

fprintf('=== DESAIN LQR SISTEM POULTRY HOUSE ===\n\n');

%% ------------------------------------------------------------------------
%  BAGIAN 1 : MODEL STATE-SPACE
%% ------------------------------------------------------------------------

% Parameter contoh (bisa kamu sesuaikan dengan jurnal)
Rcb   = 0.015;     % thermal resistance (K/W)
Kg    = 50;        % gain heating
a     = 1000;      % konstanta a
vb    = 500;       % volume udara 1
vc    = 600;       % volume udara 2
Rcv   = 0.01;      % coupling suhu -> kelembaban
x1e   = 25;        % suhu ekuilibrium (°C)
u3_op = 100;       % ventilasi nominal

% Matriks dinamik (sesuai slide kamu)
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

% Kontrolabilitas
Co = ctrb(A,B);
if rank(Co)==size(A,1)
    fprintf('✓ Sistem controllable → LQR bisa diterapkan.\n\n');
else
    error('Sistem tidak controllable.');
end

%% ------------------------------------------------------------------------
%  BAGIAN 2 : DESAIN LQR (SATU SKENARIO UTAMA: BALANCED)
%% ------------------------------------------------------------------------

% Bobot LQR (boleh kamu tulis di slide sebagai "algoritma LQR")
Q = diag([100 50]);        % prioritas suhu > kelembaban
R = diag([1 1 0.1]);       % penalti energi kontrol

[K,P,eig_cl] = lqr(A,B,Q,R);      % K gain LQR, P solusi Riccati

Acl = A - B*K;                    % matriks closed-loop

fprintf('Gain LQR K:\n'); disp(K);
fprintf('Pole closed-loop (eig(A-BK)):\n'); disp(eig_cl);

% Estimasi waktu tunak dari pole dominan
[~,idx_dom] = max(real(eig_cl));
lambda_dom  = eig_cl(idx_dom);
ts_theory   = 4/abs(real(lambda_dom));
fprintf('Estimasi Ts (teori, ~4/|Re(pole dominan)|) = %.2f s\n\n',ts_theory);

%% ------------------------------------------------------------------------
%  BAGIAN 3 : ANALISIS KESTABILAN DENGAN LYAPUNOV
%% ------------------------------------------------------------------------

fprintf('=== ANALISIS KESTABILAN (METODE LYAPUNOV) ===\n');

% 1. Cek P > 0  → V(x)=x''Px > 0
eig_P = eig(P);
fprintf('Eigen P:\n'); disp(eig_P);
if all(eig_P>0)
    fprintf('✓ P positive definite → kandidat fungsi Lyapunov valid.\n');
else
    warning('P tidak positive definite.');
end

% 2. Hitung M = Q + K''RK, harus >0 → dotV = -x''Mx < 0
M = Q + K'*R*K;
eig_M = eig(M);
fprintf('Eigen M = Q + K''RK:\n'); disp(eig_M);
if all(eig_M>0)
    fprintf('✓ M positive definite → Vdot = -x''Mx < 0 (x≠0).\n');
else
    warning('M tidak positive definite.');
end

% 3. Cek residual Riccati (harus ≈ 0)
ARE_res = A'*P + P*A - P*B*(R\B')*P + Q;
fprintf('Norm residual ARE = %.3e\n\n',norm(ARE_res,'fro'));

fprintf('→ Karena P>0 dan M>0, sistem closed-loop ASIMTOTIK STABIL (Lyapunov).\n\n');

%% ------------------------------------------------------------------------
%  BAGIAN 4 : SIMULASI CLOSED-LOOP (PERFORMANSI)
%% ------------------------------------------------------------------------

fprintf('=== SIMULASI PERFORMANSI CLOSED-LOOP LQR ===\n');

t  = 0:0.1:500;                 % waktu simulasi
x0 = [ 5; 10];                  % deviasi awal: suhu +5°C, RH +10%%

% gangguan temperatur & kelembaban (contoh)
d1 = 3*ones(size(t));           % suhu luar +3°C (konstan)
d2 = 5*ones(size(t));           % RH luar +5%
d2(t>250) = -2;                 % step ke -2% di t=250 s
dist = [d1; d2];

% Sistem tertutup terhadap gangguan d
sys_cl = ss(Acl,Bd,C,zeros(2,2));
[y,t,x] = lsim(sys_cl,dist',t,x0);  % x : state, y: keluaran = x

% Hitung sinyal kontrol u = -Kx
u = zeros(length(t),3);
for k = 1:length(t)
    u(k,:) = -K*x(k,:)';
end

%% ------------------------------------------------------------------------
%  BAGIAN 5 : METRIK PERFORMANSI
%% ------------------------------------------------------------------------

% 1) Settling time (2% dari nilai awal)
idx_ts1 = find(abs(x(:,1)) < 0.02*abs(x0(1)),1);
idx_ts2 = find(abs(x(:,2)) < 0.02*abs(x0(2)),1);

if isempty(idx_ts1), Ts1 = t(end); else, Ts1 = t(idx_ts1); end
if isempty(idx_ts2), Ts2 = t(end); else, Ts2 = t(idx_ts2); end

fprintf('Settling time (2%% band):\n');
fprintf('  x1 (suhu)      : %.2f s\n',Ts1);
fprintf('  x2 (kelembaban): %.2f s\n',Ts2);
fprintf('  Estimasi teori : %.2f s\n\n',ts_theory);

% 2) Overshoot
OS1 = (max(x(:,1)) - x0(1))/abs(x0(1))*100;
OS2 = (max(x(:,2)) - x0(2))/abs(x0(2))*100;

fprintf('Overshoot:\n');
fprintf('  x1 (suhu)      : %.2f %%\n',OS1);
fprintf('  x2 (kelembaban): %.2f %%\n\n',OS2);

% 3) Steady-state error (akhir simulasi)
Ess1 = x(end,1);
Ess2 = x(end,2);
fprintf('Steady-state error:\n');
fprintf('  x1 (suhu)      : %.4f °C\n',Ess1);
fprintf('  x2 (kelembaban): %.4f %%\n\n',Ess2);

% 4) Energi kontrol ∫ u^2 dt
E_u1 = trapz(t,u(:,1).^2);
E_u2 = trapz(t,u(:,2).^2);
E_u3 = trapz(t,u(:,3).^2);
E_tot = E_u1 + E_u2 + E_u3;

fprintf('Energi kontrol (∫u^2 dt):\n');
fprintf('  Heating (u1)   : %.2f\n',E_u1);
fprintf('  PAD (u2)       : %.2f\n',E_u2);
fprintf('  Vent (u3)      : %.2f\n',E_u3);
fprintf('  Total          : %.2f\n\n',E_tot);

% 5) Fungsi biaya LQR J = ∫ (x''Qx + u''Ru) dt
state_cost   = trapz(t, sum((x*Q).*x,2));
control_cost = trapz(t, sum((u*R).*u,2));
J = state_cost + control_cost;

fprintf('Fungsi biaya LQR:\n');
fprintf('  State cost     : %.2f\n',state_cost);
fprintf('  Control cost   : %.2f\n',control_cost);
fprintf('  Total J        : %.2f\n\n',J);

%% ------------------------------------------------------------------------
%  BAGIAN 6 : PLOT RESPON
%% ------------------------------------------------------------------------

figure('Name','Respon State','Position',[100 100 1100 400]);

subplot(1,2,1);
plot(t,x(:,1),'b','LineWidth',2); hold on;
yline(0,'r--','Setpoint');
yline( 0.02*abs(x0(1)),'k:');
yline(-0.02*abs(x0(1)),'k:');
grid on;
xlabel('Waktu (s)');
ylabel('Deviasi Suhu (°C)');
title('Respon Suhu Internal x_1');

subplot(1,2,2);
plot(t,x(:,2),'r','LineWidth',2); hold on;
yline(0,'b--','Setpoint');
yline( 0.02*abs(x0(2)),'k:');
yline(-0.02*abs(x0(2)),'k:');
grid on;
xlabel('Waktu (s)');
ylabel('Deviasi Kelembaban (%)');
title('Respon Kelembaban Internal x_2');

figure('Name','Sinyal Kontrol','Position',[100 100 900 350]);
plot(t,u,'LineWidth',2);
grid on;
xlabel('Waktu (s)');
ylabel('Control Input');
legend('u_1 (Heating)','u_2 (PAD)','u_3 (Ventilasi)');
title('Sinyal Kontrol LQR');
