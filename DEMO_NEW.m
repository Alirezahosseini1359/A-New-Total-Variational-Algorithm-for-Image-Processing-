function [out] = DEMO_NEW

addpath(genpath('./UTILS'));

%% Setting parameters
global  alpha1 gamma eta kappa Nbiter b R R_adj W m n;
Nbiter = 5;
kappa=8e-15;
eta= 1;
gamma = .75;
alpha1 = 1;
n = 2*128;
m = n;
%% Loadibg data
x0 = double(imread('knee.png')); % Load another image if needed 
x0 = imresize(x0, [n m]);
figure(1);imshow(uint8(x0));title('Ground Truth');
load('random16_256.mat');             % Load another mask if needed
mask = Q1;
figure(2);imshow(mask);title('Sampling Mask');
imwrite(mask,'SM.png', 'png');
b = mask.*fft2c(x0);
zf = ifft2c(b);
figure(3);imshow(uint8(abs(zf)));title('ZF');
imwrite(uint8(zf),'ZF.png', 'png');

%% Defining operators
R =@(x) mask.*fft2c(x);
W = Wavelet('Daubechies',4,4);
R_adj = @(r) ifft2c(mask.*r);
%% solution by TV+W 
out = TGVmain3(x0,b,alpha1,gamma,kappa,eta,R,R_adj,W,Nbiter,zf);