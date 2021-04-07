f = imread('432.tif');

%��ʹ������Ƶ�ʵ�ͨ�˲�����
[M,N] = size(f);
F = fft2(f);
H = lpfilter('gaussian',M,N,10);%��˹��ͨ�˲���
G = H.*F;
g = real(ifft2(G));%ȡʵ��
subplot(1,3,1);imshow(f);title('ԭͼ');
%ʹ������Ƶ�ʵ�ͨ�˲�����
PQ = paddedsize(size(f));          
Fp = fft2(f,PQ(1),PQ(2));                   
Hp = lpfilter('gaussian',PQ(1),PQ(2),20);%��˹��ͨ�˲���
Gp = Hp.*Fp;                                
gp = real(ifft2(Gp));
gp1 = gp(1:size(f,1),1:size(f,2));%�����ϲ��ľ����޼�Ϊԭʼ��С

subplot(1,3,1);imshow(f);title('ԭͼ');
subplot(1,3,2);imshow(g,[ ]);title('��ʹ�����ĵ�ͨ�˲�');
subplot(1,3,3);imshow(gp1,[ ]);title('ʹ�������ĵ�ͨ�˲�');

function [ H, D ] = lpfilter( type,M,N,D0,n )%���ֵ�ͨ�˲���
%LPFILTER creates the transfer function of a lowpass filter.
%   Detailed explanation goes here

%use function dftuv to set up the meshgrid arrays needed for computing 
%the required distances.
[U, V] = dftuv(M,N);
 
%compute the distances D(U,V)
D = sqrt(U.^2 + V.^2);

%begin filter computations
switch type
    case 'ideal'%�����ͨ�˲���
        H = double(D <= D0);
    case 'btw'%������˹��ͨ�˲���
        if nargin == 4
            n = 1;
        end
        H = 1./(1+(D./D0).^(2*n));
    case 'gaussian'%��˹��ͨ�˲���
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
%���ߴ纯��
if nargin == 1%��������ĸ���
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