%% define constants
% channel
Fs = 1E6;   %ADC sample freq.
Ts = 1/Fs;
%N0 = 2*2.07E-20; %Noise PSD(in W/Hz)
fd = 100;   %doppler/Hz
EbN0_dB = linspace(0,25,50);    %simulation Eb/No range
EbN0 = 10.^(EbN0_dB/10);
%L = 10^(-101/10);    %path loss
tau = [0 4]';
fdTs = fd*Ts;
P = [0.5 0.5]'; %power delay profile
% ofdm tx setup
%Pt = 0.1;   %tx power
mod_type = 64;   %modulation type (4,16,64 - QAM)
%fc = 2E9;   %carrier frequency
Ncp = 5;    %twice the delay spread
%Bc = 0;
%Tc = 0;
N = 128;
r1 = N/(N+Ncp); % the effective power
Nsym = 1;%number of ofdm symbols in one packet
%Tsym = 0.01;    %symbol interval?
% check setup
if ((N+Ncp)*fdTs > 0.01 ||(N+Ncp)*Nsym>9000 )
    warning('time-varing frequency response assumption invalid')
end
%simulation round
max_err_num = 100; %stop simulation when 500 bit errors
maxNum = 100000;%or stop when reach this number of bits
k = 1;
%% ofdm BER simulation
%transmitting energy
%E = Pt * (N+Ncp) * Ts;
Es = 1; %fix the Es, change N0
%bit error counter
BER = zeros(length(EbN0),1);
%BER simulation symbol by symbol
for j = 1:length(EbN0)
    %calculate noise power
    EsN0 = EbN0(j) * log2(mod_type) * r1;%energy per tx bit
    sigma = sqrt(1/(EsN0*2*1));
 
    %counters
    totErr = 0; % number of errors observed
    num = 0; %number of bits processed
    while (totErr <= max_err_num || num<=maxNum)
        %information bit packet
        b = randi([0 1],N*log2(mod_type),Nsym);
        %ofdm packet
        ofdm_packet = ofdm_pkt_gen(b,N,Ncp,Nsym,mod_type);%the matrix contains Nsym ofdm symbols
        %fading chanel
        tx_signal = reshape(ofdm_packet,[],1);%vector contains Nsym ofdm symbols
        [rx_signal, h] = Fading_Channel(tx_signal,tau,fdTs,P);
        %add noise
        rx_signal = rx_signal + sigma*(randn(length(rx_signal),1) + 1j * randn(length(rx_signal),1));
        %remove transient
        rx_signal = rx_signal(1:end-max(tau));
        rx_signal = reshape(rx_signal,N+Ncp,Nsym);
        Err_per_b_mat = 0;
        %demodulate each symbol
        for i = 1:Nsym
            %remove cp,fft
            y_temp = rx_signal(:,i);
            r = sqrt(1/N)*fft(y_temp(Ncp+1:end));
            %correct phase and gain
            c = fft([h(1 + (i-1)*(Ncp+N),1) 0 0 0 h(1 + (i-1)*(Ncp+N),2)],N);
            C = diag(c);
            A = diag(1./(c.*conj(c)));
            if isnan(A)
                warning('channel gain is 0')
            end
            %demodulate
            r_ = A*C'*r;
            %scatter(real(r_),imag(r_))
            s_hat = qamdemod(r_,mod_type,'UnitAveragePower',true);
            b_hat = de2bi(s_hat,log2(mod_type));
            b_hat = reshape(b_hat,[],1);
            %counting error per ofdm symbol
            Err_per_b_mat = Err_per_b_mat + sum(b_hat ~= b(:,i));
        end
        totErr = totErr + Err_per_b_mat;
        num = num + numel(b);
    end
    BER(j) = totErr/num;
end