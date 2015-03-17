function s = input_size(varargin)
% Helper function for ellipsoid method. Calculates the input size (encoding
% length) of multiple arguments (matrices, vectors, integers) according to 
% formula by Grötschel et al. 1988, p. 30.
%
% Example: input_size([1 55;1 1], 55, 1, [1 1]) = 26

s = 0;
%iterate over all input arguments
for i=1:nargin
    %iterate over all entries of vector or matrix
    for j=1:numel(varargin{i})
        % formula for input size
        s = s + 1 + ceil( log2(abs(varargin{i}(j)) + 1) );
    end
end