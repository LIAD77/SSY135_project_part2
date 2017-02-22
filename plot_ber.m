clear
modulation = 64;
%repetition = '1';
EbN0_sequence = 0:1:20;
BER_mat = zeros(3,21);
for i = 1:3
    ber_file = sprintf('./data/BER%d_%d.mat',modulation,i);
    data = load(ber_file);
    BER_mat(i,:) = data.BER;   
end
figure
f = semilogy(EbN0_sequence,BER_mat);
title(sprintf('%dQAM',modulation))
xlabel('Eb/N0(dB)')
ylabel('BER')
legend('1','1/2','1/3')
grid on
%saveas(f,sprintf('%dQAM.png',modulation));