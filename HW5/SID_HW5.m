% #1
fmin = 0.1; % hz
fmax = 10.0; % hz
dt = .01; % s
T = 30; % s
amp = 10; % N-m
Fs = 1/dt;
fu = [fmin:0.5: fmax]

% frequency sweep input/output %
[u,t,pf,f] = mksswp(amp,fmin,fmax,dt,T);

% multisine input/output %
[u2,t2,pf,f,M,ph] = mkmsswp(amp,fmin,fmax,dt,T, 1, fu);

% plotting %
figure(1)
plot(t, u)
title('Frequency-Sweep Command')
xlabel('Time (s)')
ylabel('Torque (N-m)')
figure(2)
plot(t2, u2)
title('multisine Command')
xlabel('Time (s)')
ylabel('Torque (N-m)')
% Frequency Sweep FFT %
NFFT = length(u);
Y = fft(u,NFFT);
F = ((0:1/NFFT:1-1/NFFT)*Fs).';
magnitudeY= abs(Y); % Magnitude of the FFT
phaseY = unwrap(angle(Y)); % Phase of the FFT
dB_mag=mag2db(magnitudeY);
figure(3)
subplot(2,1,1);plot(F(1:NFFT/2),dB_mag(1:NFFT/2));title('Magnitude response of Frequency-Sweep');
ylabel('Magnitude(dB)');
subplot(2,1,2);plot(F(1:NFFT/2),phaseY(1:NFFT/2));title('Phase response of Frequency-Sweep');
xlabel('Frequency in Hz')
ylabel('radians');
% Frequency Sweep FFT %
NFFT2 = length(u2);
Y2 = fft(u2,NFFT2);
F2 = ((0:1/NFFT2:1-1/NFFT2)*Fs).';
magnitudeY2= abs(Y2); % Magnitude of the FFT
phaseY2 = unwrap(angle(Y2)); % Phase of the FFT
dB_mag2=mag2db(magnitudeY2);
figure(4)
% hold on
subplot(2,1,1);plot(F2(1:NFFT2/2),dB_mag2(1:NFFT2/2));title('Magnitude response of multisine');
ylabel('Magnitude(dB)');
subplot(2,1,2);plot(F2(1:NFFT2/2),phaseY2(1:NFFT2/2));title('Phase response of multisine');
xlabel('Frequency in Hz')
ylabel('radians');
% combined Mag plot %
figure(5)
plot(F(1:NFFT/2),dB_mag(1:NFFT/2),F2(1:NFFT2/2),dB_mag2(1:NFFT2/2))
title('Magnitude response of frequency sweep and multisine');
ylabel('Magnitude(dB)');
xlabel('Frequency in Hz')
legend('frequency sweep', 'multisine')
% #2-4
fmin = 0.1; % hz
fmax = 10.0; % hz
dt = .01; % s
T = 30; % s
amp = 10; % N-m
Fs = 1/dt;
control_on = 0; % no conrol
% frequency sweep input/output %
[u,t,pf,f] = mksswp(amp,fmin,fmax,dt,T);
N = length(t);
[y,yd, ydd, u_cmd, t] = pend(u, dt, T, control_on);
% model from frequency sweep %
x = [ones(N,1), sin(y), yd, u_cmd];
T_hat = (x'*x)\x'*ydd
Y_hat = x*T_hat;
% frequency-sweep coefficient of determination %
disp('The model using frequency sweep data R squared value is:')
R_sqfw = (T_hat'*x'*ydd - N*mean(ydd)^2) / (ydd'*ydd - N*mean(ydd)^2)
% frequency sweep confidence interval %
v = (ydd - Y_hat);
s_sq = sum(v.^2)/(length(v) -length(T_hat));
T_hat_confid = 2 * sqrt(s_sq) * diag((x'*x)^-1);
disp('and the confidence intervals are (frequency sweep):')
T_hat_confid
% frequency sweep residual %
R = ydd - Y_hat;
% Normalized Regressor %
% for k=1:length(T_hat)
% init = strcat('v',num2str(k));
% init_ar = strcat('v',num2str(k));
% sjj.(init) = 0;
% qj.(init_ar) = [];
% end
% for k=1:length(T_hat)
% for i = 1:length(t)
% sjj.(strcat('v',num2str(k))) = sjj.(strcat('v',num2str(k))) + (x(i,k) - mean(x(:,k)))^2;
%
% end
% end
% for k = 1:length(T_hat)
% qj.(strcat('v',num2str(k))) = (x(:,k) - mean(x(:,k)))/sqrt(sjj.(strcat('v',num2str(k))));
% end
% xr2 = [];
% np = length(T_hat);
% np_temp = np -1;
% for k = 1:np_temp
% xr2 = [xr2, qj.(strcat('v',num2str(k)))];
% end
% T_hatnorm2 = (xr'*xr)\xr'*qj.(strcat('v',num2str(length(T_hat))))
sjj_1 = 0;
sjj_2 = 0;
sjj_3 = 0;
sjj_4 = 0;
sjj_5 = 0;
qj_1 = [];
qj_2 = [];
qj_3 = [];
qj_4 = [];
qj_5 = [];
for i = 1:length(t)
 sjj_1 = sjj_1 + (x(i,1) - mean(x(:,1)))^2;
 sjj_2 = sjj_2 + (x(i,2) - mean(x(:,2)))^2;
 sjj_3 = sjj_3 + (x(i,3) - mean(x(:,3)))^2;
 sjj_4 = sjj_4 + (x(i,4) - mean(x(:,4)))^2;
 sjj_5 = sjj_5 + (ydd(i) - mean(ydd))^2;
end
qj_1 = (x(:,1) - mean(x(:,1)))/sqrt(sjj_1);
qj_2 = (x(:,2) - mean(x(:,2)))/sqrt(sjj_2);
qj_3 = (x(:,3) - mean(x(:,3)))/sqrt(sjj_3);
qj_4 = (x(:,4) - mean(x(:,4)))/sqrt(sjj_4);
qj_5 = (ydd - mean(ydd))/sqrt(sjj_5);
xr = [qj_2, qj_3, qj_4];
T_hatnorm = (xr'*xr)\xr'*qj_5 % #3
% Frequency Sweep Plot %
figure(1)
plot(t,ydd,t,Y_hat, t, Y_hat + 2 * sqrt(s_sq),'--r', t, Y_hat - 2 * sqrt(s_sq), '--r')
title('actual and modeled output (frequency sweep) vs. time')
xlabel('time (s)')
ylabel('angular acceleration (rad/s^2)')
legend('actual','modeled')
figure(2)
plot(t,R,'.')
title('(frequency sweep)residual vs. time')
xlabel('time(s)')
ylabel('residual (rad/s^2)')
% multisine input/output %
fu = [fmin:0.5: fmax];
[u2,t2,pf,f,M,ph] = mkmsswp(amp,fmin,fmax,dt,T,1,fu);
N2 = length(t2);
[y2,yd2, ydd2, u_cmd2, t2] = pend(u2, dt, T, control_on);
% model from multisine %
x2 = [ones(N2,1), sin(y2), yd2, u_cmd2];
T_hat2 = (x2'*x2)\x2'*ydd2
Y_hat2 = x2*T_hat2;
% multisine coefficient of determination %
disp('The model using multisine data R squared value is:')
R_sqms = (T_hat2'*x2'*ydd2 - N2*mean(ydd2)^2) / (ydd2'*ydd2 - N2*mean(ydd2)^2)
% frequency sweep confidence interval %
v2 = (ydd2 - Y_hat2);
s_sq2 = sum(v2.^2)/(length(v2) -length(T_hat2));
T_hat_confid2 = 2 * sqrt(s_sq2) * diag((x2'*x2)^-1);
disp('and the confidence intervals are (multisine):')
T_hat_confid2
% multisine residual %
R2 = ydd2 - Y_hat2;
% multisine Normalized Regressor %
sjj2_1 = 0;
sjj2_2 = 0;
sjj2_3 = 0;
sjj2_4 = 0;
sjj2_5 = 0;
qj2_1 = [];
qj2_2 = [];
qj2_3 = [];
qj2_4 = [];
qj2_5 = [];
for i = 1:length(t2)
 sjj2_1 = sjj2_1 + (x2(i,1) - mean(x2(:,1)))^2;
 sjj2_2 = sjj2_2 + (x2(i,2) - mean(x2(:,2)))^2;
 sjj2_3 = sjj2_3 + (x2(i,3) - mean(x2(:,3)))^2;
 sjj2_4 = sjj2_4 + (x2(i,4) - mean(x2(:,4)))^2;
 sjj2_5 = sjj2_5 + (ydd2(i) - mean(ydd2))^2;
end
qj2_1 = (x2(:,1) - mean(x2(:,1)))/sqrt(sjj2_1);
qj2_2 = (x2(:,2) - mean(x2(:,2)))/sqrt(sjj2_2);
qj2_3 = (x2(:,3) - mean(x2(:,3)))/sqrt(sjj2_3);
qj2_4 = (x2(:,4) - mean(x2(:,4)))/sqrt(sjj2_4);
qj2_5 = (ydd2 - mean(ydd2))/sqrt(sjj2_5);
xr2 = [qj2_2, qj2_3, qj2_4];
T_hatnorm2 = (xr2'*xr2)\xr2'*qj2_5 % #3
% Frequency Sweep Plot %
figure(3)
plot(t2, ydd2, t2, Y_hat2, t2, Y_hat2 + 2 * sqrt(s_sq2),'--r', t2, Y_hat2 - 2 * sqrt(s_sq2), '--r')
title('actual and modeled output (multisine) vs. time')
xlabel('time (s)')
ylabel('angular acceleration (rad/s^2)')
legend('actual','modeled')
figure(4)
plot(t2,R2,'.')
title('(multisine)residual vs. time')
xlabel('time(s)')
ylabel('residual (rad/s^2)')