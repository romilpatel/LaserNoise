function x_filtered = bandpass_filter_signal(x, Fs, F_center, BW, nTaps)
    % Inputs:
    % x         - Input signal (real-valued)
    % Fs        - Sampling rate in Hz (e.g., 2.5e9 for 2.5 GSa/s)
    % F_center  - Center frequency of the passband in Hz (e.g., 200e6 for 200 MHz)
    % BW        - Bandwidth around the center frequency in Hz (e.g., 50e6 for 50 MHz)
    % nTaps     - Number of filter taps

    % Outputs:
    % x_filtered - Bandpass filtered signal

    % Calculate the lower and upper cutoff frequencies
    F_low = F_center - BW/2; % Lower cutoff frequency
    F_high = F_center + BW/2; % Upper cutoff frequency
    
    % Normalize the frequencies based on the Nyquist frequency
    nyquist = Fs / 2;
    f_low_norm = F_low / nyquist;   % Normalized lower cutoff frequency
    f_high_norm = F_high / nyquist; % Normalized upper cutoff frequency

    % Filter Design Parameters
    N = nTaps; % Filter order (adjustable parameter for sharper roll-off)

    % Use firpm (Parks-McClellan) to design a linear phase FIR bandpass filter
    % Specify passband and stopband frequencies (normalized by Nyquist frequency)
    frequencies = [0, f_low_norm*0.9, f_low_norm, f_high_norm, f_high_norm*1.1, 1];
    amplitudes = [0, 0, 1, 1, 0, 0]; % Desired amplitude response in each frequency band
    weights = [1, 1, 1]; % Weights for passband and stopband regions

    % Design the FIR filter using Parks-McClellan algorithm
    b = firpm(N, frequencies, amplitudes, weights);

    % Apply the bandpass filter to the input signal
    x_filtered = filter(b, 1, x);

    % Plot the frequency response of the filter (optional)
    fvtool(b, 1); % Visualize the frequency response of the filter
    
    % Plot the original and filtered signals in the frequency domain (optional)
    plot_filtered_signal_spectrum(x, x_filtered, Fs);
end

function plot_filtered_signal_spectrum(x, x_filtered, Fs)
    % This function plots the original and filtered signal spectra
    % Inputs:
    % x          - Original signal
    % x_filtered - Filtered signal
    % Fs         - Sampling rate

    % Frequency vectors for plotting
    N = length(x);
    f = linspace(0, Fs, N);

    % FFT of the original signal
    X = fft(x);

    % FFT of the filtered signal
    X_filtered = fft(x_filtered);

    % Optionally, plot a comparison of the spectra
    figure;
    hold on;
    plot(f(1:N/2)/1e6, abs(X(1:N/2))/max(abs(X)), 'r'); % Original signal (normalized)
    plot(f(1:N/2)/1e6, abs(X_filtered(1:N/2))/max(abs(X_filtered)), 'b'); % Filtered signal (normalized)
    title('Comparison of Original and Filtered Signal Spectra');
    xlabel('Frequency (MHz)');
    ylabel('Normalized Magnitude');
    legend('Original Signal', 'Filtered Signal');
    set(gca,'YScale','log')
    grid on;
end
