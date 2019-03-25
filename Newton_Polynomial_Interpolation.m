% The code can only be used for research purpose.

% Please cite the following paper when you use it:
   % Ning Li, Yongqiang Zhao, Quan Pan, and Seong G. Kong, 
   % "Demosaicking DoFP images using Newton¡¯s polynomial interpolation and polarization difference model" 
   % Optics Express 27, 1376-1391 (2019)

%Note:
%   The code is not optimized and may have bugs. There are ways to improve the efficiency of the algorithms. 
%   Your revision and improvement are welcome!
%   All the notes in this code correspond to the cases explained in the
%   original paper.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Newton Polynomial Interpolation                          %
%                                                          %
% Copyright (C) 2019 Ning Li. All rights reserved.         %
%                    ln_neo@mail.nwpu.edu.cn               %
%                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [I0,I45,I90,I135] = Newton_Polynomial_Interpolation(I)

[m,n] = size(I);

R = zeros(m,n,4);

% label the different polarization channels
O = zeros(m,n);
O(1:2:m,1:2:n) = 1;
O(1:2:m,2:2:n) = 2;
O(2:2:m,2:2:n) = 3;
O(2:2:m,1:2:n) = 4;

% Store intermediate results
Y1 = I;
Y2 = I;
R(:,:,1) = I;
R(:,:,2) = I;
R(:,:,3) = I;
R(:,:,4) = I;

% Stage one interpolation: interpolate vertically for case Fig.6(b),
% interpolate horizontally for case Fig.6(c), interpolate in diagonal
% directions for case Fig.6(a). The Eqs.(14)-(17) are simplified in this
% code.
for i = 4:m-3
    for j = 4:n-3
        R(i,j,O(i,j)) = I(i,j);
        R(i,j,O(i,j+1)) = 0.5*I(i,j) + 0.0625*I(i,j-3) - 0.25*I(i,j-2) + 0.4375*I(i,j-1) + 0.4375*I(i,j+1) - 0.25*I(i,j+2) + 0.0625*I(i,j+3);
        R(i,j,O(i+1,j)) = 0.5*I(i,j) + 0.0625*I(i-3,j) - 0.25*I(i-2,j) + 0.4375*I(i-1,j) + 0.4375*I(i+1,j) - 0.25*I(i+2,j) + 0.0625*I(i+3,j);
        Y1(i,j) = 0.5*I(i,j) + 0.0625*I(i-3,j-3) - 0.25*I(i-2,j-2) + 0.4375*I(i-1,j-1) + 0.4375*I(i+1,j+1) - 0.25*I(i+2,j+2) + 0.0625*I(i+3,j+3);
        Y2(i,j) = 0.5*I(i,j) + 0.0625*I(i-3,j+3) - 0.25*I(i-2,j+2) + 0.4375*I(i-1,j+1) + 0.4375*I(i+1,j-1) - 0.25*I(i+2,j-2) + 0.0625*I(i+3,j-3);
    end
end

% One can adjust for better result
thao = 5.8;

% Fusion of the estimations with edge classifier for case Fig.6(a).
for i = 4:m-3
    for j = 4:n-3
        pha1 = 0;
        for k = -2:2:2
            for l = -2:2:2
                pha1 = pha1 + abs(Y1(i+k,j+l) - I(i+k,j+l));
            end
        end
        pha2 = 0;
        for k = -2:2:2
            for l = -2:2:2
                pha2 = pha2 + abs(Y2(i+k,j+l) - I(i+k,j+l));
            end
        end
        if (pha1/pha2)>thao
            R(i,j,O(i+1,j+1)) = Y2(i,j);
        end
        if (pha2/pha1)>thao
            R(i,j,O(i+1,j+1)) = Y1(i,j);
        end
        if ((pha1/pha2)<thao)&&((pha2/pha1)<thao)
            d1 = abs(I(i-1,j-1) - I(i+1,j+1)) + abs(2*I(i,j) -I(i-2,j-2) - I(i+2,j+2));
            d2 = abs(I(i+1,j-1) - I(i-1,j+1)) + abs(2*I(i,j) -I(i+2,j-2) - I(i-2,j+2));
            epsl = 0.000000000000001;
            w1 = 1/(d1 + epsl);
            w2 = 1/(d2 + epsl);
            R(i,j,O(i+1,j+1)) = (w1*Y1(i,j) + w2*Y2(i,j))/(w1 + w2);
        end
    end
end

RR = R;

XX1 = I;
XX2 = I;
YY1 = I;
YY2 = I;

% Stage two interpolation: interpolate horizontally for case Fig.6(b),
% interpolate vertically for case Fig.6(c).
for i = 4:m-3
    for j = 4:n-3
        XX1(i,j) = R(i,j,O(i,j+1));
        XX2(i,j) = 0.5*I(i,j) + 0.0625*R(i-3,j,O(i,j+1)) - 0.25*I(i-2,j) + 0.4375*R(i-1,j,O(i,j+1)) + 0.4375*R(i+1,j,O(i,j+1)) - 0.25*I(i+2,j) + 0.0625*R(i+3,j,O(i,j+1));
        YY1(i,j) = R(i,j,O(i+1,j));
        YY2(i,j) = 0.5*I(i,j) + 0.0625*R(i,j-3,O(i+1,j)) - 0.25*I(i,j-2) + 0.4375*R(i,j-1,O(i+1,j)) + 0.4375*R(i,j+1,O(i+1,j)) - 0.25*I(i,j+2) + 0.0625*R(i,j+3,O(i+1,j));
    end
end

% Fusion of the estimations with edge classifier for case Fig.6(b) and Fig.6(c).
for i = 4:m-3
    for j = 4:n-3
        pha1 = 0;
        for k = -2:2:2
            for l = -2:2:2
                pha1 = pha1 + abs(XX1(i+k,j+l) - I(i+k,j+l));
            end
        end
        pha2 = 0;
        for k = -2:2:2
            for l = -2:2:2
                pha2 = pha2 + abs(XX2(i+k,j+l) - I(i+k,j+l));
            end
        end
        if (pha1/pha2)>thao
            RR(i,j,O(i,j+1)) = XX2(i,j);
        end
        if (pha2/pha1)>thao
            RR(i,j,O(i,j+1)) = XX1(i,j);
        end
        if ((pha1/pha2)<thao)&&((pha2/pha1)<thao)
            d1 = abs(I(i,j-1) - I(i,j+1)) + abs(2*I(i,j) -I(i,j-2) - I(i,j+2));
            d2 = abs(I(i+1,j) - I(i-1,j)) + abs(2*I(i,j) -I(i+2,j) - I(i-2,j));
            epsl = 0.000000000000001;
            w1 = 1/(d1 + epsl);
            w2 = 1/(d2 + epsl);
            RR(i,j,O(i,j+1)) = (w1*XX1(i,j) + w2*XX2(i,j))/(w1 + w2);
        end
        pha1 = 0;
        for k = -2:2:2
            for l = -2:2:2
                pha1 = pha1 + abs(YY1(i+k,j+l) - I(i+k,j+l));
            end
        end
        pha2 = 0;
        for k = -2:2:2
            for l = -2:2:2
                pha2 = pha2 + abs(YY2(i+k,j+l) - I(i+k,j+l));
            end
        end
        if (pha1/pha2)>thao
            RR(i,j,O(i+1,j)) = YY2(i,j);
        end
        if (pha2/pha1)>thao
            RR(i,j,O(i+1,j)) = YY1(i,j);
        end
        if ((pha1/pha2)<thao)&&((pha2/pha1)<thao)
            d1 = abs(I(i,j-1) - I(i,j+1)) + abs(2*I(i,j) -I(i,j-2) - I(i,j+2));
            d2 = abs(I(i+1,j) - I(i-1,j)) + abs(2*I(i,j) -I(i+2,j) - I(i-2,j));
            epsl = 0.000000000000001;
            w1 = 1/(d1 + epsl);
            w2 = 1/(d2 + epsl);
            RR(i,j,O(i+1,j)) = (w1*YY2(i,j) + w2*YY1(i,j))/(w1 + w2);
        end
    end
end
R = RR;
I0 = R(:,:,1);
I45 = R(:,:,2);
I90 = R(:,:,3);
I135 = R(:,:,4);
end