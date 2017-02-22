function decoded = repdec(coded, repetitions, gain)
%REPDEC - Repetition decoder of a sequence of symbols. The output is interleaved.
% coded - The repetition-coded symbols.
% repetitions - The number of repetitions.
% gain - The channel gain coefficients corresponding to each coded symbol.

if size(coded, 1) ~= numel(coded)
    error('Coded must be a solumn vector.')
end

if mod(repetitions, 1) ~= 0
    error('Repetitions must be integer.')
end

if mod(numel(coded), repetitions) ~= 0
    error('The length of coded must be a multiple of repetitions.')
end

if size(coded) ~= size(gain)
    error('coded and gain must be of the same size.')
end

deinterleaved = reshape(coded, numel(coded) / repetitions, repetitions);
gain = reshape(gain, numel(gain) / repetitions, repetitions);
decoded = zeros(numel(coded) / repetitions, 1);
for index = 1:size(deinterleaved, 1)
    % gain(index, :)
    % deinterleaved(index, :)
    % deinterleaved(index, :) * gain(index, :)' / sum(gain(index, :))
    % pause(1)
    decoded(index) = deinterleaved(index, :) * gain(index, :)' / sum(gain(index, :));
end

end