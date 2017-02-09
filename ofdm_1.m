%% define constants
% channel
fc = 2E9;   %carrier
N0 = 2*2.07E-20; %Noise PSD(in W/Hz)
fd = 100;   %doppler/Hz
EbNo_dB = linspace(0,25,50);    %simulation Eb/No range
EbNo = 10.^EbNo_dB;
mod_type = 4;   %modulation type (4,16,64 - QAM)
L = 10^(101/10);    %path loss
Pt = 0.1;   %tx power
Fs = 1E6;   %ADC sample freq.
Ts = 1/Fs;
% ofdm
Ncp = 4;
Bc = 0;
Tc = 0;
N = 64;
Nsym = 8;
%% ofdm symbol generator
%generate random bits
ofdm_symbols = zeros(Nsym,N+Ncp);%each column is an OFDM symbol
for i = 1:Nsym
    b = randi([0 1],1,N*log2(mod_type));
    ofdm_symbols(:,i) = ofdm_sym_gen(b,N,Ncp,mod_type,E,Ts );
end
%% time varing channel
