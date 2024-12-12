function estimate_linewidth(scopeCapturedSignal,Fs,taud_d)


    %% Input parser


    %% Initial vatiable calculations
    I = scopeCapturedSignal;
    Q = imag(hilbert(scopeCapturedSignal));

    % time and frequency for measurement data (is shorter since leading delay period is erased): N_meas =   - l
    dt        = 1/Fs;
    N_meas    = length(I);
    t_meas    = [0:N_meas-1]*dt;
    f_meas    = [-N_meas/2 : N_meas/2-1]/(dt*N_meas);

    dt_meas   = dt;
    l_meas    = round(tau_d/dt);
    tau_d     = l_meas*dt;

    %% convolution function
    % convolution function h(t) and Fourier transform H(\omega) for phase signal
    h                   = zeros(N_meas,1);
    h(1)                = 1;  % instantaneous signal
    h(N_meas+1-l_meas)  = -1; % retarded signal

    H                   = transpose(fftshift(fft(h)));

    % direct
    Hsq = 2*(1 - cos(2*pi*f_meas.*par.tau_d));


    %% Get phase signal
    phase_signal = atan(Q./I); % + par.deltaOmega * t_meas - Omega_CW * par.tau_d;
    
    %% Correct phase jump
    phase_signal = correct_phase_jumps(phase_signal);




end