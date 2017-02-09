function [ z ] = ofdm_sym_gen( b,N,Ncp,mod_type,E,Ts )
%OFDM_SYM Summary of this function goes here
%   Detailed explanation goes here
%   generate one ofdm symbol
%   b: bit string
%   N: number of sub-carriers
%   Ncp: length of CP
%   mod_type: modulation type
%   k_norm: power normalization factor
%   E: energe per symbol?
%   Ts: sample interval
%map bit to decimal
bs = reshape(b(:),log2(mod_type),N);
sb = bi2de(bs');
%map decimal to MQAM symbol
%don't need to calculate k_norm here
sp = sqrt(E) * qammod(sb,mod_type,'UnitAveragePower',true);
%ifft
z = sqrt(N/Ts) * ifft(sp);
%add cp
z = [z(end-Ncp+1:end); z];
end

