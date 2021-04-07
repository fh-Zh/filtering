I=imread('3_3.jpg');
f=double(I);
F=fft2(f);%��ԭͼ����ά����Ҷ�任
S=fftshift(log(1+abs(F)));
figure,imshow(S,[]);title('ԭͼƵ��');
 
h=fspecial ('sobel')  %3*3��ֱ����Sobel�˲���
figure,freqz2(h);title('sobel���ӵ�Ƶ����Ӧ');    %��ʾsobel���ӵ�Ƶ����Ӧ

PQ= paddedsize(size(f));%�õ�������ĳߴ�
H=freqz2(h, PQ(1), PQ(2));%������h��Ӧ��Ƶ���˲���
H1=ifftshift(H);          
gs=imfilter(double(f), h);%�ռ����˲��õ�ͼ��
gf=dftfilt(f, H1);%Ƶ���˲��õ�ͼ��
figure, subplot(1,2,1);imshow(gf, [ ]);title('ͼ2��Ƶ���˲�');%Ƶ���˲�
subplot(1,2,2);imshow(gs, [ ]);title('ͼ3���ռ����˲�');%ֱ���ڿռ����˲�

d=abs(gs-gf);     %�Կռ��˲������Ƶ���˲����ȡ��ֵ
disp(['�ռ��˲������Ƶ���˲�����Ĳ�ֵ������:'] );  
max(d(:))
disp('��ֵ��СΪ��')
min(d(:))
disp('���ϣ�����֤�ռ����Ƶ���˲��ȼ�')

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

function g = dftfilt(f,H)
% Ƶ���˲�����
F=fft2(f,size(H,1),size(H,2));
g=real(ifft2(H.*F));
g=g(1:size(f,1),1:size(f,2));
end
