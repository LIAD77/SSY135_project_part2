%plot all BER curves in one
load('th_BER_4.mat');
load('th_BER_16.mat')
load('th_BER_64.mat')
load('BER_measured.mat');
EbN0_1 = linspace(0,25,50);
EbN0_2 = th_BER_16.data{1,1};
figure
semilogy(EbN0_1,[BER_4_5_128';BER_16_5_128';BER_64_5_128'])
hold on
semilogy(EbN0_2,[th_BER_4.data{1,2};th_BER_16.data{1,2};th_BER_64.data{1,2}])
grid on
legend('QPSK(s)','16QAM(s)','64QAM(s)','QPSK(t)','16QAM(t)','64QAM(t)')