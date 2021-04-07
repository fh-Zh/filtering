f = imread('432.tif');

%不使用填充的频率低通滤波处理
[M,N] = size(f);
F = fft2(f);
H = lpfilter('gaussian',M,N,10);%高斯低通滤波器
G = H.*F;
g = real(ifft2(G));%取实部
subplot(1,3,1);imshow(f);title('原图');
%使用填充的频率低通滤波处理
PQ = paddedsize(size(f));          
Fp = fft2(f,PQ(1),PQ(2));                   
Hp = lpfilter('gaussian',PQ(1),PQ(2),20);%高斯低通滤波器
Gp = Hp.*Fp;                                
gp = real(ifft2(Gp));
gp1 = gp(1:size(f,1),1:size(f,2));%将左上部的矩形修剪为原始大小

subplot(1,3,1);imshow(f);title('原图');
subplot(1,3,2);imshow(g,[ ]);title('不使用填充的低通滤波');
subplot(1,3,3);imshow(gp1,[ ]);title('使用零填充的低通滤波');

function [ H, D ] = lpfilter( type,M,N,D0,n )%三种低通滤波器
%LPFILTER creates the transfer function of a lowpass filter.
%   Detailed explanation goes here

%use function dftuv to set up the meshgrid arrays needed for computing 
%the required distances.
[U, V] = dftuv(M,N);
 
%compute the distances D(U,V)
D = sqrt(U.^2 + V.^2);

%begin filter computations
switch type
    case 'ideal'%理想低通滤波器
        H = double(D <= D0);
    case 'btw'%巴特沃斯低通滤波器
        if nargin == 4
            n = 1;
        end
        H = 1./(1+(D./D0).^(2*n));
    case 'gaussian'%高斯低通滤波器
        H = exp(-(D.^2)./(2*(D0^2)));
    otherwise 
        error('Unkown filter type');

end
end

function [U, V] = dftuv(M, N)
%DFTUV Computes meshgrid frequency matrices.
%   [U, V] = DFTUV(M, N) computes meshgrid frequency matrices U and
%   V.  U and V are useful for computing frequency-domain filter
%   functions that can be used with DFTFILT.  U and V are both
%   M-by-N.
 
%   Copyright 2002-2004 R. C. Gonzalez, R. E. Woods, & S. L. Eddins
%   Digital Image Processing Using MATLAB, Prentice-Hall, 2004
%   $Revision: 1.3 $  $Date: 2003/04/16 22:30:34 $
 
% Set up range of variables.
u = 0:(M - 1);
v = 0:(N - 1);
 
% Compute the indices for use in meshgrid.
idx = find(u > M/2);
u(idx) = u(idx) - M;
idy = find(v > N/2);
v(idy) = v(idy) - N;
 
% Compute the meshgrid arrays.
[V, U] = meshgrid(v, u);
end

function PQ = paddedsize(AB,CD,~ )
%填充尺寸函数
if nargin == 1%输入变量的个数
    PQ = 2*AB;
elseif nargin ==2 && ~ischar(CD)
    PQ = QB +CD -1;
    PQ = 2*ceil(PQ/2);
elseif nargin == 2
    m = max(AB);
    P = 2^nextpow(2*m);
    PQ = [P,P];
elseif nargin == 3
    m = max([AB CD]);
    P = 2^nextpow(2*m);
    PQ = [P,P];
else 
    error('Wrong number of inputs');
end
end
