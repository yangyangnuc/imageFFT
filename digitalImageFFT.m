%% FFT of digital image
% IDFT is inverse discrete fourier transform
% DFT is discrete fourier transform
% f(x,y) denotes for spatial pixels of image
% F(u,v) denotes for frquency coefficients, having same size of f
% F(0,0) is direct current component of F(u,v), so F(0,0) =
% M*N*mean(f(x,y)). 
% M is source image row number
% N is source image column number
% Often using R(u,v), I(u,v) denote for real and image componet of F(u,v)
% So fourier spectra of image, |F(u,v)| = sqrt( R(u,v)^2 + I(u,v)^2 )
% angle of transform is Fai(u,v) = arctan( I(u,v)/R(u,v) )
% F(u,v) = |F(u,v)|*e^(-j*Fai(u,v))
% Power spectra is defined as power(|F(u,v)|,2), distribution of power 
% values as a function of frequency
% conjugate symmetric means F(u,v) = F(-u,-v)
% it can be proven that











clc
close all


src = imread('originalFluorescenceImage.png');
size(src)
src = rgb2gray(src);
size(src)
figure, subplot 331, imshow(src,[0,255]),title('Original image');

F = fft2(src);

% calculate amplitude of each element, including real and image part
S = abs(F);
S = 255*( S-min(min(S)) )/(  max(max(S)) - min(min(S))  );
% 
subplot 332, imshow(S),title('Fourier specture');

