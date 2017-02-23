%calculate the theoratical BER for task B
M = [4 16 64]; %modulation
%average EbN0 of the sum of all branches
avg_EbN0 = 10.^([19.65 16.64 14.88]/10);
N = [1 2 3];    %diversity order
A = 0;
B = pi/2;
%define constants and functions
%c1  = @(M) 4/(log2(M)*pi);
%M_rb = @(s,M,avg_rb) (1-log2(M).^2 * avg_rb).^(-1);
%c2_phi = @(phi,M) 3*log2(M)./(2*(M-1).*sin(phi).^2);
to_be_calc = @(phi,M,N,avg_rb) 4/(pi*log2(M)) * (1+ avg_rb * 3 ./ (2*(M-1)*sin(phi).^2)).^-N ;
%cauclate the integral
for i = 1:3
    for j = 1:3
        % note that the EbN0 per branch here
        q = integral(@(phi) to_be_calc(phi,M(i),N(j),avg_EbN0(i)/N(j)),A,B);
        q_ = berfading(avg_EbN0(i)/N(j),'qam',M(i),N(j));
        disp(sprintf('Mod %d div %d BER %f vs berfading %f\n',M(i),N(j),q,q_));        
    end
end