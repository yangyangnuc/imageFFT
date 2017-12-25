%% FFT of digital image
clc
close all


src = imread('originalFluorescenceImage.png');
size(src)
src = rgb2gray(src);
size(src)
figure, subplot 331, imshow(src,[0,255]),title('Original image');

F = fft2(src);
subplot 332, imshow(F),title('Fourier specture');

