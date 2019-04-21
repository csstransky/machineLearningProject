function ret = testInRectanglesAxisParallel(P1, P2, P3, P4, Img);

%                    L4
%        P2-----------P3
%  L1  |         p         |    
%        |                    |L3
%       P1----------P4
%                 L2

xMin = P1(1,1);
xMax = P3(1,1);
yMin = P1(1,2);
yMax = P2(1,2);

ret = 1;
for i = 1:size(Img,1)    
    if (Img(i,1) < xMin || Img(i,1) > xMax || Img(i,2) < yMin  || Img(i,2) > yMax)
        ret = 0;
        break;
    end        
end