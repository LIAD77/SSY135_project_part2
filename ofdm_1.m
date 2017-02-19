%simulation of coded OFDM transmission:
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
mod_type = 4;   %modulation type (4,16,64 - QAM)
%fc = 2E9;   %carrier frequency
Ncp = 5;    %twice the delay spread
%Bc = 0;
%Tc = 0;
N = 128;
r1 = N/(N+Ncp); % the effective power
diversity = 3;
r2 = 1/diversity; %code rate
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
    EsN0 = EbN0(j) * log2(mod_type);%energy per info bit
    sigma = sqrt(1/(EsN0*2*1 * r1 *r2));
    N0 = sigma^2*2;
    %counters
    totErr = 0; % number of errors observed
    num = 0; %number of bits processed
    while (totErr <= max_err_num || num<=maxNum)
        %information bit packet
        b = randi([0 1],N*log2(mod_type),1);
        %ofdm packet
        ofdm_symbol = r2 * ofdm_sym_gen(b,N,Ncp,mod_type);%equally divide the tx power to 3 ofdm symbols
        %fading chanel
        rx_signal = zeros(length(ofdm_symbol)+max(tau),diversity); % the tx signal (each column is a branch)
        h = zeros(2,diversity); %the inpulse response of each branch
        for branch = 1:diversity
            [rx_signal(:,branch), h_temp] = Fading_Channel(ofdm_symbol,tau,fdTs,P);
            h(:,branch) = transpose(h_temp(1,:));
        end
        %add noise
        rx_signal = rx_signal + sigma*(randn(size(rx_signal)) + 1j * randn(size(rx_signal)));
        %remove transient
        rx_signal = rx_signal(1:end-max(tau),:);
        %remove cp
        rx_signal = rx_signal(Ncp+1:end,:);
        %fft
        y = sqrt(1/N)*fft(rx_signal);
        Err_per_b_mat = 0;
        %initialize some variables
        c = zeros(size(y)); %the frequency response of each branch
        for i = 1:diversity
            c(:,i) = fft([h(1,i) 0 0 0 h(2,i)],N);
        end
        mrc_r = zeros(N,1);%MRC cofficients
        A = zeros(N,1);%scaling factor
        r_ = zeros(N,1); %after MRC
        %MRC and decoding
        for i = 1:diversity
            mrc_r = conj(c(:,i)) ;%./ (c(:,i).*conj(c(:,i)));
            A = A + r2 * c(:,i).*conj(c(:,i));
            if isnan(mrc_r)
                warning('channel gain is 0')
            end
            %MRC
            r_ = r_ + mrc_r .* y(:,i);
            %scatter(real(r_),imag(r_))
        end
        r_ = r_ ./ A;
        %scaling factor
        s_hat = qamdemod(r_,mod_type,'UnitAveragePower',true);
        b_hat = de2bi(s_hat,log2(mod_type));
        b_hat = reshape(b_hat,[],1);
        %counting error per ofdm symbol
        Err_per_b_mat = Err_per_b_mat + sum(b_hat ~= b);
        
        totErr = totErr + Err_per_b_mat;
        num = num + numel(b);
    end
    BER(j) = totErr/num;
end