function m = findImagePointMean(I, P1, P2, P3, P4)
%             L4
%        P1---------P4
%    L1  |    p     |    
%        |          |L3
%       P2----------P3
%            L2

s = 0;
nP = 0;

[ht wd] = size(I);

[minX, minY, maxX, maxY] = findMinMax(P1, P2, P3, P4, wd, ht);

for i = minX:maxX
    for j = minY:maxY
        P = [i j];
        if isInsideRectangle(P1, P2, P3, P4, P)
            s = s+I(i,j);            
            nP = nP+1;
            
            %plot(j,i,'.b'); hold on;
        end
    end
end

m = s/nP;