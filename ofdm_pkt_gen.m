function [ ofdm_symbol_mat ] = ofdm_pkt_gen( b,N,Ncp,Nsym,mod_type )
%OFDM_PKT_GEN Summary of this function goes here
%   Detailed explanation goes here
%   Generate ofdm packet
%   b: bit matrix, each column is one ofdm symbol
%   N: # of sub carriers
%   Ncp: # of cyclic prefix
%   Nsym = # of OFDM symbols per packet
%   mod_type: 4,16,64 QAM
%   output: ofdm packet, each column is an ofdm symbol
ofdm_symbol_mat = zeros(N+Ncp,Nsym);
        for i = 1:Nsym
            %ofdm generator
            ofdm_symbol_mat(:,i) = ofdm_sym_gen(b(:,i),N,Ncp,mod_type);
        end
end

