% The code can only be used for research purpose.

% Please cite the following paper when you use it:
   % Ning Li, Yongqiang Zhao, Quan Pan, and Seong G. Kong, 
   % "Demosaicking DoFP images using Newton¡¯s polynomial interpolation and polarization difference model" 
   % Optics Express 27, 1376-1391 (2019)

%Note:
%   The code is not optimized and may have bugs. There are ways to improve the efficiency of the algorithms. 
%   Your revision and improvement are welcome!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Newton Polynomial Interpolation                          %
%                                                          %
% Copyright (C) 2019 Ning Li. All rights reserved.         %
%                    ln_neo@mail.nwpu.edu.cn               %
%                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
close all;
clear all;
%%%%%%%%%%1. Read real DoFP image %%%%%%% 
I = double(imread('test.bmp'));

%%%%%%%%%%2. Demosaicking DoFP image with Newton's polynomial intepolation method %%%%%% 
[PI0,PI45,PI90,PI135] = Newton_Polynomial_Interpolation(I);

%%%%%%%%%%3. Calculating the Stokes parameters %%%% 
ii = 0.5*(PI0 + PI45 + PI90 + PI135);
S0P = im2uint8(mat2gray(ii));
q = PI0 - PI90;
u = PI45 - PI135;
dolp = sqrt(q.*q + u.*u);
Dolp = dolp./ii;
DOLPP = im2uint8(mat2gray(Dolp));
aop = (1/2) * atan(u./q);
AOPP = im2uint8(mat2gray(aop));
figure
imshow(S0P)
title('S0')
figure
imshow(DOLPP)
title('DoLP')
figure
imshow(AOPP)
title('AoP')
