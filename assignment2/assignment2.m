f = imread('436.tif');
D0=20;

f = mat2gray(f,[0 255]);
subplot(2,4,1);imshow(f);title('a');
% 1.给定一幅大小为M*N的输入图像f(x,y),得到填充参数P = 2M,Q = 2N
[m,n] = size(f);
P = 2 * m;
Q = 2 * n; 
% 2.对f(x,y)添加必要数量的0，形成大小为P * Q的填充后的图像fp(x,y)
fp = zeros(P,Q);
fp(1:m,1:n) = f(1:m,1:n);
subplot(2,4,2);imshow(fp);title('b');
% 3.用(-1)^(x+y)乘以fp(x,y)移到其交换的中心
for i0 = 1 : m
    for j0 = 1 : n 
        fp(i0,j0) = double(fp(i0,j0)*(-1)^(i0+j0));
    end
end
subplot(2,4,3);imshow(fp);title('c');
% 4.计算来自步骤3的图像的DFT，得到F(u,v)
F = fft2(fp,P,Q);
subplot(2,4,4);imshow(F);title('d');
% 5.生成一个实的、对称的滤波函数H(u,v),其大小为P*Q，中心在（P/2,Q/2）处;用阵列相乘形成乘积G(u,v) = H(u,v)F(u,v);即G(i,k)=H(i,k)F(i,k)
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
% 6.得到处理后的图像
gp = ifft2(G); 
gp = real(gp);
subplot(2,4,7);imshow(gp);title('g');
for i0 = 1 : m
    for j0 = 1 : n 
        gp(i0,j0) = double(gp(i0,j0)*(-1)^(i0+j0));
    end
end
% 7.通过从gp(x,y)的左上象限提取M*N区域，得到最终处理结果g(x,y)
g=zeros(m,n);
g(1:m,1:n) = gp(1:m,1:n);
subplot(2,4,8);imshow(g);title('h');