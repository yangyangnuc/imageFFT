%% FFT of digital image
% ---IDFT is inverse discrete fourier transform
% ---DFT is discrete fourier transform
% ---f(x,y) denotes for spatial pixels of image
% ---F(u,v) denotes for frquency coefficients, having same size of f
% ---F(0,0) is direct current component of F(u,v), so F(0,0) =
% M*N*mean(f(x,y)). 
% ---M is source image row number
% ---N is source image column number
% ---Often using R(u,v), I(u,v) denote for real and image componet of F(u,v)
% So fourier spectra of image, |F(u,v)| = sqrt( R(u,v)^2 + I(u,v)^2 )
% ---angle of transform is Fai(u,v) = arctan( I(u,v)/R(u,v) )
% ---F(u,v) = |F(u,v)|*e^(-j*Fai(u,v))
% ---Power spectra is defined as power(|F(u,v)|,2), distribution of power 
% values as a function of frequency
% ---conjugate symmetric means F(u,v) = F(-u,-v)
% ---it can be proven that F(u,v) = F(u+k1*M,v+k2*N), but DFT only calculate
% one cycle, so we only deal with M*N image. k1, k2 are integers.
% center the spectra of FFT











clc
close all
clearvars


src = imread('originalFluorescenceImage.png');
size(src)
src = rgb2gray(src);
size(src)
figure, subplot 331, imshow(src,[]),title('Original image');

src = double(src);

F = fft2(src);
for i=1:size(src,1)
    for j=1:size(src,2)
        src(i,j) = power(-1,i+j-2)*src(i,j);
    end
end

F1 = fft2(src);
% calculate amplitude of each element, including real and image part
S = abs(F);
Sc = fftshift(S);
Sc_log = log(1+abs(Sc));
S = 255*( S-min(min(S)) )/(  max(max(S)) - min(min(S))  );
S_log = log(1 + abs(S));

S1 = abs(F1);
% S1 = 255*( S1-min(min(S1)) )/(  max(max(S1)) - min(min(S1))  );
S1_log = log(1 + abs(S1));

% 
subplot 332, imshow(S),title('Log mode Fourier specture');
subplot 333, imshow(S1_log,[]),title('Self-design Log Mode Fourier specture centered');
subplot 334, imshow(Sc_log,[]),title('Log mode Fourier specture build-in centered');
% figure, plot(-100:100,S(101,1:201))

%% demonstrate f(x)*(-1)^(x), then FFT, translate frquency 
% domain's original to u = M/2
x = 1:100;
z = (-1)*ones(1,size(x,2));
y1 = x.*x; 
y = x.*x.*power(z,x);
figure, plot(x,y1),title('orginal x,y'), hold on;
Y = fft(y);
Yc = fftshift(Y)
Y1 = fft(y1);
figure, plot(x,Yc),title('centered frequency spectra');
figure, plot(x,Y1),title('raw in frequency'), 
hold on;
plot(x,Y),title('translate to M/2, in frquency domain'), 
legend('raw in frequency','translate to M/2, in frquency domain');
hold off;



%% 

