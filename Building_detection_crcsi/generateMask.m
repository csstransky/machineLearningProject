function [pMask sMask Rm] = generateMask(A, R, Ad, Rd, ALS, LPDx, LPDy);

resolutionMask = 0.25;
%resolutionImage = abs(R(2,1));
%resolutionDEM = abs(Rd(2,1));

%assume image is wholey covered by both DEM and ALS data
%so LLpoint and URpoint of image will be LLpoint and URpoint for the masks
[height width bands] = size(A);
pMax = [1 width 1]*R; %URpoint
pMin = [height 1 1]*R; %LLpoint

%pMin = max([R(3,1); Rd(3,1); min(ALS(:,1))]);
%pMax = min([R(3,2); Rd(3,2); max(ALS(:,2))]);

%just to show
%figure; mapshow(A(:,:,1:3),R); hold on;
%plot(pMin(1,1),pMin(1,2), '+b'); hold on;
%plot(pMax(1,1),pMax(1,2), '+y'); hold off;

%TFW for masks
Rm(1,1) = 0;
Rm(2,1) = resolutionMask;
Rm(3,1) = pMin(1,1);
Rm(1,2) = -resolutionMask;
Rm(2,2) = 0;
Rm(3,2) = pMax(1,2);

%define masks
wMask = double(uint16(floor((pMax(1,1) - pMin(1,1)) / resolutionMask) + 1.0));
hMask = double(uint16(floor((pMax(1,2) - pMin(1,2)) / resolutionMask) + 1.0));
if ((wMask <= 0) || (hMask <= 0)) 
    return;
end
pMask = zeros(hMask,wMask) == 1;
sMask = zeros(hMask,wMask) == 0;

%just to show

%on object space
%pMaxMask = [1 wMask 1]*Rm; %URpoint
%pMinMask = [hMask 1 1]*Rm; %LLpoint
%figure; mapshow(pMask,Rm); hold on;
%plot(pMinMask(1,1),pMinMask(1,2), '+b'); hold on;
%plot(pMaxMask(1,1),pMaxMask(1,2), '+y'); hold off;

%on image space
%Pm = obj2obs(pMin,Rm);
%PM = obj2obs(pMax,Rm);
%figure; imshow(pMask); hold on;
%plot(Pm(1,2),Pm(1,1), '+b'); hold on;
%plot(PM(1,2),PM(1,1), '+y'); hold off;

%calculate neighbourhood size (dx or dy) in mask with respect to ALS data
resolution = floor (LPDx + LPDy + 0.5) * 0.5;
if (resolution < 0.5) 
    resolution = 0.5;
end
nbr	= double(resolution/resolutionMask)-1;

ALSabove = []; ALSg = [];

[sxd syd] = size(Ad);

%figure; mapshow(A(:,:,1:3),R); hold on;
for I = 1:size(ALS,1)%nTileX
    P = ALS(I,1:2);
    H = ALS(I,3);
    %find 2D position on DEM
    %y = floor((P(1,1)*Rd(1,2) - Rd(3,1)*Rd(1,2) + Rd(3,2)*Rd(1,1) - P(1,2)*Rd(1,1))/(Rd(2,1)*Rd(1,2)-Rd(2,2)*Rd(1,1)));
    %x = floor((P(1,2) - Rd(3,2) - y*Rd(2,2))/Rd(1,2));
    Pd = obj2obs(P,Rd);
    x = Pd(1,1); y = Pd(1,2);
    %ind = [x y]
    if (x > 0 && x <= sxd && y > 0 && y <= syd)
        Hd = Ad(x,y);

        nn = 0;
        while (Hd == -9999)
            nn = nn+1;
            stx = x-nn;
            enx = x+nn;
            sty = y-nn;
            eny = y+nn;
            if (stx < 1)
                stx = 1;
            end
            if (sty < 1)
                sty = 1;
            end
            if (enx > sxd)
                enx = sxd;
            end
            if (eny > syd)
                eny = syd;
            end
            Hd = myMean(Ad(stx:enx, sty:eny));
        end

        gThresh = Hd + 2.5;
        if (H >= gThresh)
            ALSabove = [ALSabove; ALS(I,:)];
            %plot(P(1,1),P(1,2),'+g'); hold on;
        else
            ALSg = [ALSg; ALS(I,:)];            
        end
    end
    %thresholding = I
end
%hold off;
%figure; mapshow(A,R);

for i = 1:size(ALSabove,1)
    P = ALSabove(i,1:2);
    %y = floor((P(1,1)*Rm(1,2) - Rm(3,1)*Rm(1,2) + Rm(3,2)*Rm(1,1) - P(1,2)*Rm(1,1))/(Rm(2,1)*Rm(1,2)-Rm(2,2)*Rm(1,1)));
    %x = floor((P(1,2) - Rm(3,2) - y*Rm(2,2))/Rm(1,2));
    Pm = obj2obs(P,Rm);
    x = Pm(1,1); y = Pm(1,2);
    
    if (x-nbr >0 && y-nbr >0 && x+nbr <= hMask && y+nbr <= wMask)
        sMask(x-nbr:x+nbr,y-nbr:y+nbr) = 0;
    end     
    
    
end
        
        
for i = 1:size(ALSg,1)
     P = ALSg(i,1:2);
     %y = floor((P(1,1)*Rm(1,2) - Rm(3,1)*Rm(1,2) + Rm(3,2)*Rm(1,1) - P(1,2)*Rm(1,1))/(Rm(2,1)*Rm(1,2)-Rm(2,2)*Rm(1,1)));
     %x = floor((P(1,2) - Rm(3,2) - y*Rm(2,2))/Rm(1,2));
     Pm = obj2obs(P,Rm);
     x = Pm(1,1); y = Pm(1,2);
    
     if (x-nbr >0 && y-nbr >0 && x+nbr <= hMask && y+nbr <= wMask)
        pMask(x-nbr:x+nbr,y-nbr:y+nbr) = 1;        
     end
            
end

%show
%figure; imshow(pMask);
%figure; imshow(sMask);

%write the masks
imwrite(sMask, 'SecondaryMask.jpg','JPG');
imwrite(pMask, 'PrimaryMask.jpg','JPG');
%here = 1;





function m = myMean(M);
[sx sy] = size(M);
s = 0; c = 0;
for i = 1:sx
    for j = 1:sy
        if (M(i,j) ~= -9999)
            s = s+M(i,j);
            c = c+1;
        end
    end
end
    
if (c > 0)
    m = s/c;
else
    m = -9999;
end