function [Pxx_custom, f] = custom_PSD(signal, Fs, NFFT, window, noverlap)
    % custom_PSD computes the power spectral density using Welch's method manually.
    %
    % Inputs:
    %   - signal: The input signal (vector)
    %   - Fs: The sampling frequency of the signal (scalar)
    %   - NFFT: The number of FFT points
    %   - window: The window function (vector)
    %   - noverlap: Number of overlapping samples between segments (scalar)
    %
    % Outputs:
    %   - Pxx_custom: The custom computed Power Spectral Density
    %   - f: Frequency vector corresponding to Pxx_custom

    % Ensure the signal is a column vector
    signal = signal(:);
    
    % Length of the window
    win_length = length(window);
    
    % Number of segments
    step = win_length - noverlap;
    num_segments = floor((length(signal) - noverlap) / step);
    
    % Initialize variables
    Pxx_custom = zeros(NFFT, 1);
    
    % Loop through segments
    for i = 1:num_segments
        % Extract segment
        segment_start = (i-1) * step + 1;
        segment_end = segment_start + win_length - 1;
        segment = signal(segment_start:segment_end);
        
        % Apply window function
        windowed_segment = segment .* window;
        
        % Compute the FFT of the windowed segment
        segment_fft = fft(windowed_segment, NFFT);
        
        % Compute the periodogram (normalized power of FFT)
        segment_periodogram = (abs(segment_fft).^2) / (win_length * Fs);
        
        % Accumulate the periodograms
        Pxx_custom = Pxx_custom + segment_periodogram;
    end
    
    % Average the periodograms
    Pxx_custom = Pxx_custom / num_segments;
    
    % Extract the one-sided PSD (only positive frequencies)
    Pxx_custom = Pxx_custom(1:NFFT/2+1);
    
    % Frequency vector
    f = Fs * (0:(NFFT/2)) / NFFT;
    
    % Normalize PSD to account for window power
    Pxx_custom = Pxx_custom / (sum(window.^2) / win_length);
end

function compare_PSD(signal, Fs)
    % compare_PSD compares the custom PSD function with MATLAB's pwelch function.
    %
    % Inputs:
    %   - signal: The input signal (vector)
    %   - Fs: The sampling frequency of the signal (scalar)
    
    % Parameters
    NFFT = 2^nextpow2(length(signal));  % Number of FFT points
    window = hamming(256);  % Hamming window of length 256
    noverlap = 128;  % 50% overlap
    
    % Custom Welch PSD
    [Pxx_custom, f_custom] = custom_PSD(signal, Fs, NFFT, window, noverlap);
    
    % MATLAB's pwelch function
    [Pxx_pwelch, f_pwelch] = pwelch(signal, window, noverlap, NFFT, Fs);
    
    % Plot comparison
    figure;
    subplot(2,1,1);
    plot(f_custom, 10*log10(Pxx_custom));
    title('Custom PSD Estimate');
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');
    grid on;
    
    subplot(2,1,2);
    plot(f_pwelch, 10*log10(Pxx_pwelch));
    title('MATLAB pwelch PSD Estimate');
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');
    grid on;
end
