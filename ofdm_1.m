%% define constants
% channel
Fs = 1E6;   %ADC sample freq.
Ts = 1/Fs;
N0 = 2*2.07E-20; %Noise PSD(in W/Hz)
fd = 100;   %doppler/Hz
EbN0_dB = linspace(0,25,50);    %simulation Eb/No range
EbN0 = 10.^EbN0_dB;
L = 10^(-101/10);    %path loss
tau = [0 4]';
fdTs = fd*Ts;
P = [0.5 0.5]'; %power delay profile
% ofdm tx setup
Pt = 0.1;   %tx power
mod_type = 4;   %modulation type (4,16,64 - QAM)
fc = 2E9;   %carrier frequency
Ncp = 8;    %twice the delay spread
Bc = 0;
Tc = 0;
N = 64;
Nsym = 8;
Tsym = 0.01;    %symbol interval?
% check setup
if (N+Ncp)*fdTs > 0.01
    warning('time-varing frequency response assumption invalid')
end
%% ofdm BER simulation
%transmitting energy
E = Pt * (N+Ncp) * Ts;
%bit error counter
err = zeros(length(EbN0),1);
%BER simulation symbol by symbol
for j = 1:length(EbN0)
    %calculate noise power
    sigma = 0;
    for i = 1:Nsym
        %bit
        b = randi([0 1],1,N*log2(mod_type));
        %ofdm generator
        ofdm_symbol = ofdm_sym_gen(b,N,Ncp,mod_type,E,Ts );
        %path loss
        ofdm_symbol = ofdm_symbol * L;
        %fading
        [y, h] = Fading_Channel(ofdm_symbol,tau,fdTs,P);
        %add noise
        y = y + (randn(length(y),1)*sigma + 1j * randn(length(y),1)*sigma);
        %remove cp
        y = y(Ncp+1:end);
        %fft
        r = sqrt(Ts/N)*fft(y);
        %correct phase and gain
        C = diag(fft(h(:,1),N));
        A = 1;
        %demodulate
        s_hat = qamdemod(A*C'*r,mod_type,'PlotConstellation',true);
        b_hat = de2bi(s_hat,log2(mod_type));
        b_hat = reshape(b_hat,1,[]);
        %counting error
        err(j) = err(j) + sum(b_hat ~= b);
    end
end