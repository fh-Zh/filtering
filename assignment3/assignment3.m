I=imread('3_3.jpg');
f=double(I);
F=fft2(f);%对原图做二维傅立叶变换
S=fftshift(log(1+abs(F)));
figure,imshow(S,[]);title('原图频域');
 
h=fspecial ('sobel')  %3*3垂直方向Sobel滤波器
figure,freqz2(h);title('sobel算子的频域响应');    %显示sobel算子的频域响应

PQ= paddedsize(size(f));%得到零填充后的尺寸
H=freqz2(h, PQ(1), PQ(2));%生成与h相应的频域滤波器
H1=ifftshift(H);          
gs=imfilter(double(f), h);%空间域滤波得到图像
gf=dftfilt(f, H1);%频域滤波得到图像
figure, subplot(1,2,1);imshow(gf, [ ]);title('图2：频域滤波');%频域滤波
subplot(1,2,2);imshow(gs, [ ]);title('图3：空间域滤波');%直接在空间域滤波

d=abs(gs-gf);     %对空间滤波结果和频域滤波结果取差值
disp(['空间滤波结果和频域滤波结果的差值不大于:'] );  
max(d(:))
disp('差值最小为：')
min(d(:))
disp('综上，可验证空间域和频域滤波等价')

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

function g = dftfilt(f,H)
% 频域滤波函数
F=fft2(f,size(H,1),size(H,2));
g=real(ifft2(H.*F));
g=g(1:size(f,1),1:size(f,2));
end
