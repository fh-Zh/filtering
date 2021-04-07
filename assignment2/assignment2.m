f = imread('436.tif');
D0=20;

f = mat2gray(f,[0 255]);
subplot(2,4,1);imshow(f);title('a');
% 1.����һ����СΪM*N������ͼ��f(x,y),�õ�������P = 2M,Q = 2N
[m,n] = size(f);
P = 2 * m;
Q = 2 * n; 
% 2.��f(x,y)��ӱ�Ҫ������0���γɴ�СΪP * Q�������ͼ��fp(x,y)
fp = zeros(P,Q);
fp(1:m,1:n) = f(1:m,1:n);
subplot(2,4,2);imshow(fp);title('b');
% 3.��(-1)^(x+y)����fp(x,y)�Ƶ��佻��������
for i0 = 1 : m
    for j0 = 1 : n 
        fp(i0,j0) = double(fp(i0,j0)*(-1)^(i0+j0));
    end
end
subplot(2,4,3);imshow(fp);title('c');
% 4.�������Բ���3��ͼ���DFT���õ�F(u,v)
F = fft2(fp,P,Q);
subplot(2,4,4);imshow(F);title('d');
% 5.����һ��ʵ�ġ��ԳƵ��˲�����H(u,v),���СΪP*Q�������ڣ�P/2,Q/2����;����������γɳ˻�G(u,v) = H(u,v)F(u,v);��G(i,k)=H(i,k)F(i,k)
H = zeros(P,Q);
a = 2 * D0^2;
for u = 1 :P
    for v = 1:Q
        D = (u-(m+1.0))^2+(v-(n+1.0))^2;
        H(u,v) = exp((-D)/a);
    end
end
G = H.*F;
subplot(2,4,5);imshow(H);title('e');
subplot(2,4,6);imshow(G);title('f');
% 6.�õ�������ͼ��
gp = ifft2(G); 
gp = real(gp);
subplot(2,4,7);imshow(gp);title('g');
for i0 = 1 : m
    for j0 = 1 : n 
        gp(i0,j0) = double(gp(i0,j0)*(-1)^(i0+j0));
    end
end
% 7.ͨ����gp(x,y)������������ȡM*N���򣬵õ����մ�����g(x,y)
g=zeros(m,n);
g(1:m,1:n) = gp(1:m,1:n);
subplot(2,4,8);imshow(g);title('h');