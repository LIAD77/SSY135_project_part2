clear; close;

% ADC sample frequency
Fs = 1e6;
Ts = 1/Fs;

% Doppler/Hz
fd = 100;
fdTs = fd*Ts;

% Channel tap delays
tap_delay = [0 4]';
max_delay = max(tap_delay);

% Set the cyclic prefix length
prefix_length = max_delay;

% Power delay profile
power_delay_profile = [0.5 0.5]';

% Set the number of channels
channels = floor(1 / 50 / fdTs - prefix_length);

% The number of symbols
modulation_order = 16;
num_symbols = 10 * channels;

% Simulation Eb/No range
EbN0_sequence = [0:1:25];
EbN0_sequence = [0:1:2];
BER = zeros(size(EbN0_sequence));
MAX_RUNS = 1e5;
MAX_ERRORS = 1e5;

for EbN0_index = 1:length(EbN0_sequence)
    disp(['Running simulation for ', num2str(EbN0_sequence(EbN0_index)), 'dB.'])
    iteration_count = 0;
    iteration_errors = 0;
    while iteration_count < MAX_RUNS && iteration_errors < MAX_ERRORS
        iteration_count = iteration_count + 1;

        % Generate random input
        bits_tx = randi([0 1], num_symbols * log2(modulation_order), 1);

        % Modulate
        symbols_tx = qammod(bits_tx, modulation_order, 'InputType', 'bit', 'UnitAveragePower', true);

        % Encode with a repetition coder
        coded_tx = repenc(symbols_tx, 1);

        % Separate the input symbols into channels
        parallel_tx = parallelize(coded_tx, channels);

        % Create the OFDM symbols
        ofdm_symbols_tx = zeros(size(parallel_tx));
        for col = 1:size(parallel_tx, 2)
            ofdm_symbols_tx(:, col) = sqrt(channels) * ifft(parallel_tx(:, col));
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
        % Add noise
        EbN0 = EbN0_sequence(EbN0_index);
        channel_output = add_awgn(channel_output, EbN0, channels / (prefix_length + channels), log2(modulation_order));

        % Remove the max_delay last symbols
        channel_output = channel_output(1:end - max_delay, :);

        % Parallelize the channel output
        ofdm_symbols_prefix_rx = parallelize(channel_output, channels + prefix_length);

        % Remove the prefix
        ofdm_symbols_rx = zeros(size(ofdm_symbols_prefix_rx, 1) - prefix_length, ...
                                size(ofdm_symbols_prefix_rx, 2));
        for col = 1:size(ofdm_symbols_rx, 2)
            ofdm_symbols_rx(:, col) = remove_cp(ofdm_symbols_prefix_rx(:, col), prefix_length);
        end

        % Store the signal amplitude for use in the decoder
        channel_gain = zeros(size(ofdm_symbols_rx));

        % Compute the symbol-wise fft
        parallel_rx = zeros(size(ofdm_symbols_rx));
        for col = 1:size(parallel_rx, 2)
            parallel_rx(:, col) = sqrt(1 / channels) * fft(ofdm_symbols_rx(:, col));

            % Verify that the rows of the tap process corresponding to this OFDM symbol are equal
            first_row = 1 + (channels + prefix_length) * (col - 1);
            last_row = 1 + (channels + prefix_length) * col;
            for row = first_row:last_row
                if sum(abs(tap_process(row, :) - tap_process(first_row, :)) > eps) > 10
                    error('Not all rows of the tap process are equal!')
                end
            end

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
            signal_amplitude = diag(phase_correction_matrix' * phase_correction_matrix);
            amplitude_correction_matrix = diag(1 ./ signal_amplitude);

            channel_gain(:, col) = signal_amplitude;

            parallel_rx(:, col) = phase_correction_matrix' * parallel_rx(:, col);
            parallel_rx(:, col) = amplitude_correction_matrix * parallel_rx(:, col);
        end

        % Serialize the symbols
        coded_rx = serialize(parallel_rx);
        channel_gain = serialize(channel_gain);

        % Decode the repetition code
        symbols_rx = repdec(coded_rx, 1, channel_gain);

        % Demodulate into a bit sequence
        bits_rx = qamdemod(symbols_rx, modulation_order, 'OutputType', 'bit');
        bits_rx = bits_rx(1:length(bits_tx));

        % Compute the error rate
        iteration_errors = iteration_errors + sum(abs(bits_tx - bits_rx) ~= 0);
    end

    BER(EbN0_index) = iteration_errors / (iteration_count * length(bits_tx));
end

figure;
plot(BER)
grid on;
return

% Plot the results
% TODO: Remove
figure; hold on;

subplot(2,2,1); hold on;
plot(real(symbols_tx), imag(symbols_tx), 'b.')
plot(real(symbols_rx), imag(symbols_rx), 'ro')
legend('Transmitted Symbols', 'Received Symbols')
grid on; hold off;

subplot(2,2,2); hold on;
plot(abs(channel_input(1:channels + prefix_length)))
plot(abs(channel_output(1:channels + prefix_length)))
legend('Transmitted Signal', 'Received Signal')
grid on; hold off;

subplot(2,2,3); hold on;
plot(abs(parallel_tx(:, 1)))
plot(abs(parallel_rx(:, 1)))
grid on; hold off;