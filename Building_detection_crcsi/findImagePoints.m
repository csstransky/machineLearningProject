function Pixels = findImagePoints(I, P1, P2, P3, P4)
%             L4
%        P1---------P4
%    L1  |    p     |    
%        |          |L3
%       P2----------P3
%            L2

Pixels = [];

[ht wd band] = size(I);

[minX, minY, maxX, maxY] = findMinMax(P1, P2, P3, P4, wd, ht);

for i = minX:maxX
    for j = minY:maxY
        P = [i j];
        if isInsideRectangle(P1, P2, P3, P4, P)
            if band == 1
                Pixels = [Pixels; I(i,j,1)]; 
            elseif band == 2
                Pixels = [Pixels; [I(i,j,1) I(i,j,2)]]; 
            elseif band == 3
                Pixels = [Pixels; [I(i,j,1) I(i,j,2) I(i,j,3)]];                        
            end
            %plot(j,i,'.b'); hold on;
        end
    end
end
