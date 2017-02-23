clear
modulation = 4;
%repetition = '1';
EbN0_sequence = 0:1:20;
BER_mat = zeros(6,21);
for i = 1:2:6
    ber_file = sprintf('./data/BER%d_%d.mat',modulation,(i+1)/2);
    ber_file_bertool = sprintf('./data/BER%d_%d_bertool.mat',modulation,(i+1)/2);
    data = load(ber_file);
    data_bertool = load(ber_file_bertool);
    BER_mat(i,:) = data.BER;
    BER_mat(i+1,:) = data_bertool.bers0.data{1,2};
end
figure
f = semilogy(EbN0_sequence,BER_mat);
title(sprintf('%dQAM',modulation))
xlabel('Eb/N0(dB)')
ylabel('BER')
legend('1(s)','1(t)','1/2(s)','1/2(t)','1/3(s)','1/3(t)')
grid on
%saveas(f,sprintf('%dQAM.png',modulation));