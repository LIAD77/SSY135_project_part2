% ADC sample frequency
Fs = 1e6;
Ts = 1/Fs;

% Doppler/Hz
fd = 100;
fdTs = fd*Ts;

% Simulation Eb/No range
% EbN0_dB = linspace(0,25,50);
% EbN0 = 10.^(EbN0_dB/10);

% Channel tap delays
tap_delay = [0 4]';
max_delay = max(tap_delay);

% Set the number of channels
channels = 64;

% Set the cyclic prefix length
prefix_length = 8;

% Power delay profile
power_delay_profile = [0.5 0.5]';

% The number of symbols
modulation_order = 16;
num_symbols = 20 * channels;

% Generate random input
bits_tx = randi([0 modulation_order - 1], 1, num_symbols);

% Modulate
symbols_tx = qammod(bits_tx, modulation_order, 'UnitAveragePower',true)';

% Separate the input symbols into channels
parallel_tx = parallelize(symbols_tx, channels);

% Create the OFDM symbols
ofdm_symbols_tx = zeros(size(parallel_tx));
for col = 1:size(parallel_tx, 2)
    ofdm_symbols_tx(:, col) = sqrt(channels / Ts) * ifft(parallel_tx(:, col));
end

% Add the cyclic prefix
ofdm_symbols_prefix_tx = zeros(size(ofdm_symbols_tx, 1) + prefix_length, ...
                               size(ofdm_symbols_tx, 2));
for col = 1:size(parallel_tx, 2)
    ofdm_symbols_prefix_tx(:, col) = add_cp(ofdm_symbols_tx(:, col), prefix_length);
end

% Serialize the ofdm symbols
channel_input = serialize(ofdm_symbols_prefix_tx);

% Pass the input through the channel
[channel_output, tap_process] = Fading_Channel(channel_input, tap_delay, fdTs, ...
                                               power_delay_profile);

channel_output = channel_output(1:end - max_delay, :);

% Verify that the rows of the tap process are equal
for row = 2:size(tap_process, 1)
    if sum(abs(tap_process(1, :) - tap_process(row, :)) > eps) ~= 0
        %error('Not all rows of the tap process are equal!')
    end
end

% Parallelize the channel output
ofdm_symbols_prefix_rx = parallelize(channel_output, channels + prefix_length);

% Remove the prefix
ofdm_symbols_rx = zeros(size(ofdm_symbols_prefix_rx, 1) - prefix_length, ...
                        size(ofdm_symbols_prefix_rx, 2));
for col = 1:size(ofdm_symbols_rx, 2)
    ofdm_symbols_rx(:, col) = remove_cp(ofdm_symbols_prefix_rx(:, col), prefix_length);
end

% Compute the symbol-wise fft
parallel_rx = zeros(size(ofdm_symbols_rx));
for col = 1:size(parallel_rx, 2)
    parallel_rx(:, col) = sqrt(Ts / channels) * fft(ofdm_symbols_rx(:, col));

    % Compute the channel impulse response for this OFDM symbols from the tap process
    channel_response = zeros(1, max_delay);
    for index = 1:length(tap_delay)
        delay = tap_delay(index);
        channel_response(delay + 1) = ...
            tap_process(1 + (channels + prefix_length) * (col - 1), index);
    end

    % Compute the phase correction matrix
    phase_correction_matrix = diag(fft(channel_response, channels));

    % Compute the amplitude correction matrix
    amplitude_correction_matrix = diag(phase_correction_matrix' * phase_correction_matrix);
    amplitude_correction_matrix = diag(1 ./ amplitude_correction_matrix);

    parallel_rx(:, col) = phase_correction_matrix' * parallel_rx(:, col);
    parallel_rx(:, col) = amplitude_correction_matrix * parallel_rx(:, col);
end

% Serialize the symbols
symbols_rx = serialize(parallel_rx);

figure; hold on;
% plot(abs(channel_input))
% plot(abs(channel_output))

subplot(1,2,1);
plot(real(symbols_tx), imag(symbols_tx), 'o')
grid on;

subplot(1,2,2);
plot(real(symbols_rx), imag(symbols_rx), 'o')
grid on;