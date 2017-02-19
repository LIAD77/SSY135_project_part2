function output = add_awgn(input, EbN0, rate, bits_per_symbol)
%ADD_AWGN - Add AWGN noise to a signal.
% input - The input signal. Assumed to have unit power.
% EbN0 - The target SNR in power per source bit. Measured in dB.
% rate - The overall code rate.
EbN0_W = 10 ^ (EbN0 / 10);
EsN0_W = EbN0_W * bits_per_symbol  * rate;
N0 = sqrt(1 / (EsN0_W * 2));
noise = N0 * (randn(length(input), 1) + 1i * randn(length(input), 1));
output = input + noise;
end