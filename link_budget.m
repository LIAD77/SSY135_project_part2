function [EbN0_dB,R] = link_budget(N,Ncp,mod_type)
%calculate the average Eb/N0 PER INFORMATION BIT at the receiver
Ts = 1E-6;
M = log2(mod_type);
r1 = N/(N+Ncp); %loss due to cyclic prefix
r2 = 1; %code rate
%path loss
L = 10^(101/10);
%rx power
Pt = 0.1 * r1 * r2 / L;
Es = Pt * Ts;
Eb = Es/M;
N0 = 2.07E-20 * 2;
EbN0_dB = 10*log10(Eb/N0);
%data rate
R = (N*M)/((N+Ncp) * Ts);
end