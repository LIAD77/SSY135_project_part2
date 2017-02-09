%% define constants
% channel
Fs = 1E6;   %ADC sample freq.
Ts = 1/Fs;
N0 = 2*2.07E-20; %Noise PSD(in W/Hz)
fd = 100;   %doppler/Hz
EbNo_dB = linspace(0,25,50);    %simulation Eb/No range
EbNo = 10.^EbNo_dB;
L = 10^(101/10);    %path loss
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
%% ofdm BER simulation
%transmitting energy
E = Pt * (N+Ncp) * Ts;
%bit error counter
err = zeros(length(EbN0),1);
%BER simulation symbol by symbol
for j = 1:length(EbNo)
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
        [r, h] = Fading_Channel(ofdm_symbol,tau,fdTs,P);
        %add noise
        r = r + (randn(length(r),1)*sigma + 1j * randn(length(r),1)*sigma);
        %remove cp
        %correct phase and gain
        %demodulate
        %counting error
    end
end