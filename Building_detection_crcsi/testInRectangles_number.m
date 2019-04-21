function ret = testInRectangles_number(P1, P2, P3, P4, Img);

%                    L4
%        P1-----------P4
%  L1  |         p         |    
%        |                    |L3
%       P2----------P3
%                 L2

% m and c for L1
m1 = (P2(1,2)-P1(1,2))/(P2(1,1)-P1(1,1));
c1 = P1(1,2) - m1*P1(1,1);

% m and c for L2
m2 = (P3(1,2)-P2(1,2))/(P3(1,1)-P2(1,1));
c2 = P2(1,2) - m2*P2(1,1);

% m and c for L3
m3 = (P4(1,2)-P3(1,2))/(P4(1,1)-P3(1,1));
c3 = P3(1,2) - m3*P3(1,1);

% m and c for L4
m4 = (P1(1,2)-P4(1,2))/(P1(1,1)-P4(1,1));
c4 = P4(1,2) - m4*P4(1,1);

Y1 = m1*Img(:,1) + c1;
Y_diff1 = Y1-Img(:,2);
        
Y3 = m3*Img(:,1) + c3;
Y_diff3 = Y3-Img(:,2);
        
Img1 = Img(Y_diff1 <= 0 & Y_diff3 > 0,:);
if (size(Img1,1) < 1)
    Img1 = Img(Y_diff1 > 0 & Y_diff3 <= 0,:);
end
        
Y2 = m2*Img1(:,1) + c2;
Y_diff2 = Y2-Img1(:,2);

Y4 = m4*Img1(:,1) + c4;
Y_diff4 = Y4-Img1(:,2);

Img2 = Img1(Y_diff2 <= 0 & Y_diff4 > 0,:);
if (size(Img2,1) < 1)
    Img2 = Img1(Y_diff2 > 0 & Y_diff4 <= 0,:);
end

ret = size(Img2,1);
%if (size(Img2,1) < 1)
%    ret = 0;
%else
%    ret = size(Img2,1);
%end