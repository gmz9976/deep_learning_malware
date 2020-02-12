function g = gistGabor(img, param)
% 
% Input:
%   img = input image (it can be a block: [nrows, ncols, c, Nimages])
%   param.w = number of windows (w*w)
%   param.G = precomputed transfer functions
%
% Output:
%   g: are the global features = [Nfeatures Nimages], 
%                    Nfeatures = w*w*Nfilters*c

img = single(img);

w = param.numberBlocks;
G = param.G;
be = param.boundaryExtension;

if ndims(img)==2
    c = 1; 
    N = 1;
    [nrows ncols c] = size(img);
end
if ndims(img)==3
    [nrows ncols c] = size(img);
    N = c;
end
if ndims(img)==4
    [nrows ncols c N] = size(img);
    img = reshape(img, [nrows ncols c*N]);
    N = c*N;
end

[ny nx Nfilters] = size(G);
W = w*w;
g = zeros([W*Nfilters N]);

% pad image
img = padarray(img, [be be], 'symmetric');

img = single(fft2(img)); 
k=0;
for n = 1:Nfilters
    ig = abs(ifft2(img.*repmat(G(:,:,n), [1 1 N]))); 
    ig = ig(be+1:ny-be, be+1:nx-be, :);
    
    v = downN(ig, w);
    g(k+1:k+W,:) = reshape(v, [W N]);
    k = k + W;
    drawnow
end

if c == 3
    % If the input was a color image, then reshape 'g' so that one column
    % is one images output:
    g = reshape(g, [size(g,1)*3 size(g,2)/3]);
end