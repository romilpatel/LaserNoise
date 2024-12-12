% function [OUT] = estimate_noise_variance(phase_signal,H,Fs)

clear all
close all
clc

%% Input Parser
data_all = load('\\FS1\Docs5\romil.patel\My Documents\GitHub\Tyndall\prog\data\laserChar\AOM_200MHz\C1330mA PP00000.dat');

addpath(genpath('..\..\LaserNoise\'))

%% 
amplitude_signal(1,:) = data_all(:,2)-mean(data_all(:,2));
nUpSample             = 2;
amp_sig1              = resample(amplitude_signal,nUpSample,1);

%% Apply bandpass filter
Fs = nUpSample*(1/(abs(data_all(2,1)-data_all(1,1))));

nTaps = 700;
% amp_sig2 = bandpass_filter_signal(amp_sig1, Fs, 200e6, 60e6, nTaps);
% amp_sig3 = bandpass_filter_signal(amp_sig2, Fs, 200e6, 80e6, nTaps);
% amp_sig = amp_sig2(1*nTaps+1:end-1*nTaps);

amp_sig = amp_sig1;
dt = 1/Fs;
N_meas = length(amp_sig);
f_meas = (-N_meas/2 : N_meas/2-1)/(dt*N_meas);
t_meas = (0:N_meas-1)*dt;

%%
IQ = hilbert(amp_sig);
I = real(IQ);
Q = imag(IQ);

phase_signal = unwrap(angle(IQ) - (2*pi*t_meas*200e6));

% phase_signal = correct_phase_jumps(phase_signal);

[p11,f11] = pwelch(phase_signal, hann(1e6), [], [], Fs,'centered','psd');
figure
loglog(f11,p11)
xlim([min(f11) 1000e6])
xline(200e6,'r:')

%%
S_zz   = (2*pi*f_meas).^2 .* PSD(phase_signal, 1/Fs, 1); % phase_signal =  dphi(t) - dphi(t-tau)   + noise
figure
loglog(f_meas,S_zz)
xlim([min(f11) 1e9])
xline(200e6,'r:')

%%
Fn_signal = (1/(2*pi)) * gradient(phase_signal,1/Fs);
[p11,f11] = pwelch(Fn_signal, hann(1e6), [], [], Fs,'centered','psd');
% close all
figure
loglog(f11,p11)
xlim([min(f11) 1e9])
xline(200e6,'r:')
title('Frequency Noise')
ylabel('FN-PSD [Hz^{2}/Hz]')
xlabel('Frequency [Hz]')

%%
tau_d               = fiber_delay(1e3);
l_meas              = round(tau_d/dt);
tau_d               = l_meas*dt;
     

h                   = zeros(N_meas,1);
h(1)                = 1;  % instantaneous signal
h(N_meas+1-l_meas)  = -1; % retarded signal
H                   = transpose(fftshift(fft(h)));

% direct
tau_d1 = 2.5*1e-6;
Hsq = 2*(1 - cos(2*pi*f_meas.*tau_d));

figure
loglog(f_meas,(abs(Hsq)))
% xlim(floor(length(Hsq))/2+[1 1e4])

%% DSH measurement raw data
S_zz   = (2*pi*f_meas).^2 .* PSD(phase_signal, 1/Fs,1); % phase_signal =  dphi(t) - dphi(t-tau)   + noise
close all     

loglog(f_meas,S_zz)

% find pole indices
idx = find(abs(Hsq).^2 < 1e-8 & f_meas > 0);
% plot values of S_zz at pole indices
figure
plot(f_meas(idx), S_zz(idx), 'bo', 'DisplayName', 'S_{z,z} (f_{pole})')
set(gca,'XScale','log')
set(gca,'YScale','log')
hold on
% xlim([1e6 200e6])

%%
% quadratic least squares fit
X = f_meas(idx);
Y = S_zz(idx);
sigma_fit = sqrt( sum(X.^2 .* Y)./sum(X.^4) );

plot(X, sigma_fit^2 .* X.^2, 'm-','LineWidth',2, 'DisplayName', 'quadratic fit: \sigma^2 f^2')

legend('Location','northwest')
title(['fit \sigma = ',num2str(sigma_fit)])
set(gca,'XScale','log')
set(gca,'YScale','log')

%% Debug
plot(f_meas(f_meas>0), sigma_fit^2 * f_meas(f_meas>0).^2, 'm-','LineWidth',2, 'DisplayName','noise fit')
% 
% end