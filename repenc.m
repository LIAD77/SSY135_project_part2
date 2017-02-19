function coded = repenc(source, repetitions)
%REPENC - Repetition encode a sequence of symbols. The output is interleaved.
% source - The source symbols.
% repetitions - The number of repetitions.

if size(source, 1) ~= numel(source)
    error('Input must be a solumn vector.')
end

if mod(repetitions, 1) ~= 0
    error('Repetitions must be integer.')
end

coded = repmat(source, 1, repetitions);
coded = reshape(coded, numel(coded), 1);

end