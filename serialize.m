function serial = serialize(parallel)
%SERIALIZE Serialize a vector by reshaping it into a column vector
%by taking elements column-wise from the matrix.
% parallel - Input matrix.

serial = reshape(parallel, numel(parallel), 1);
% serial = reshape(parallel', numel(parallel), 1);

end