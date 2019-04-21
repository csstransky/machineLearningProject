function [totalP numValP] = findImagePointCounts(I, P1, P2, P3, P4, pVal)
%             L4
%        P1---------P4
%    L1  |    p     |    
%        |          |L3
%       P2----------P3
%            L2

totalP = 0;
numValP = 0;

[ht wd] = size(I);

[minX, minY, maxX, maxY] = findMinMax(P1, P2, P3, P4, wd, ht);

for i = minX:maxX
    for j = minY:maxY
        P = [i j];
        if isInsideRectangle(P1, P2, P3, P4, P)
            totalP = totalP+1;
            if I(i,j) == pVal
                numValP = numValP+1;
            end
            %plot(j,i,'.b'); hold on;
        end
    end
end