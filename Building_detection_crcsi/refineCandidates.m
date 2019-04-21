function refineCandidates(A, R, orinetation, edgeimage, fnamein, fnameout)

%path(path,'C:\Awrangjeb\Drive I\Sydney-ALS\Digital Ortho Photo');
%path(path,'C:\Awrangjeb\Drive I\Sydney-ALS\LAS Files');
%path(path,'C:\Develop\Knoch');
%path(path,'C:\Awrangjeb\Drive I\Hobert');
%path(path,'F:\Office Desktop - 22 June 2011\Develop\Knoch');
%path(path,'F:\Office Desktop - 22 June 2011\Awrangjeb\Drive I\Sydney-ALS\Matlab Code');

%read image 
%[A, R, bbox] = geotiffread('e3525n58030_RT.tif');
%[An, Rn, bboxn] = geotiffread('crop7_NDVI.tif');
%[sx sy sz] = size(A);

%read data
fp = fopen(fnamein,'r');
Buildings = fscanf(fp,'%f %f',[9 Inf]);
fclose(fp);
Buildings = Buildings';

%read ground mask
%maskg = imread('mask_up_522000e_5251000n_LB_01_dem1_pixels.jpg');

%{
%read NDVI_sigma data
fp = fopen('image_ndvi_crop7.txt', 'r');
size2 = fscanf(fp, '%f %f %f', [1 3]);
N = fscanf(fp, '%f %f %f', [3 inf]);
fclose(fp);
N = N';


%read texture info
fp = fopen('image_texture_crop7.txt','r');
T = fscanf(fp, '%f %f %d %d', [4 inf]);
fclose(fp);
T = T';

%}

%{
fid = fopen('gradient_e3525n58030_RT.txt','r');
grad = fscanf(fid,'%d\t%d\t%f\t%d',[4 inf]);
fclose(fid);
grad = grad';
%[grad_image grad_BW] = find_gradient_image(grad, sx, sy); % a pixel value of -9999 indicates no value
x = sx/2 - grad(:,2);
y = sy/2 + grad(:,1);
O = ones(size(x,1),1);
grad(:,1:2) = [x y O]*R;
%}

nB = 1;
%binDist = 5;
%Hist = [];
%H = zeros(1, 180/binDist);
%fid = fopen('Texture_building_info_522000e_5251000n_LB_01.txt','w');
fp = fopen(fnameout,'w');

minMax2MeanRatio = 4.0;
minMax2MeanRatio2 = 3.0;
minArea2TextureRatio = 45.0;
%minMax2TextureRatio = 0.20;

Resolution = abs(R(2,1)); %0.15m for FF and HB, 0.10m for MV & Knox
minAggregateBinHeight = 9/Resolution; %9m = 90 pixels for knox & MV, 60 pixels for FF
minAggregateBinHeight2 = 6/Resolution;%6m = 60 pixels for knox & MV, 40 pixels for FF
minBinHeight = 3/Resolution; %3m = 30 pixels for knox & MV, 20 pixels for FF
maxTotal = 90/Resolution; %90m = 900 pixels for knox & MV, 600 pixels for FF
%maxTotal2 = 600; %60m = 600 pixels for knox & MV, 400 pixels for FF


%Ig = rgb2gray(uint8(A(:,:,1:3)));
%BW = edge(Ig,'canny');
figure; mapshow(A(:,:,1:3),R); hold on;
%figure; mapshow(mat2gray(edgeimage),R); hold on;
%figure; mapshow(grad_image,R); hold on;
%figure; mapshow(maskg,R); hold on;
for i = 1:size(Buildings,1)
    %if (Buildings(i,9))
        P1 = Buildings(i,1:2);  P2 = Buildings(i,3:4);  P3 = Buildings(i,5:6);  P4 = Buildings(i,7:8); 
      
        %plot([P1(1,1) P2(1,1)], [P1(1,2) P2(1,2)],'-b','LineWidth', 2); hold on;
        %    plot([P2(1,1) P3(1,1)], [P2(1,2) P3(1,2)],'-b','LineWidth', 2); hold on;
        %    plot([P3(1,1) P4(1,1)], [P3(1,2) P4(1,2)],'-b','LineWidth', 2); hold on;
        %    plot([P4(1,1) P1(1,1)], [P4(1,2) P1(1,2)],'-b','LineWidth', 2); hold on;
            
        [ret retInd H] = validate_candidate_building(A,R,orinetation, edgeimage,P1,P2,P3,P4,minBinHeight,minAggregateBinHeight);
        
        Lp = sqrt((P1(1,1)-P2(1,1))*(P1(1,1)-P2(1,1)) + (P1(1,2)-P2(1,2))*(P1(1,2)-P2(1,2)))/Resolution;
        Wp = sqrt((P3(1,1)-P2(1,1))*(P3(1,1)-P2(1,1)) + (P3(1,2)-P2(1,2))*(P3(1,2)-P2(1,2)))/Resolution;
        area2TexturePixelRatio = (Lp*Wp)/ret(1,7);
        %perMaxima = sum(H(H>minBinHeight))/ret(1,7);
        
        
        %CP = [(P1(1,1)+P3(1,1))/2 (P1(1,2)+P3(1,2))/2];
        
        trueBuilding = 0;
        
        if (size(ret,1)>0 && size(retInd,1)==0)
            trueBuilding = ((area2TexturePixelRatio>= minArea2TextureRatio) || (sum(H>=minAggregateBinHeight2) > 1 && ret(1,4) >= minMax2MeanRatio2) || (sum(H>=minAggregateBinHeight) > 1 && ret(1,8) <= minAggregateBinHeight2) || (ret(1,2) >= minAggregateBinHeight && ret(1,7) <= maxTotal) || (ret(1,4) >= minMax2MeanRatio && ret(1,2) >= minBinHeight));
        elseif (size(ret,1)>0 && size(retInd,1)>0)
            trueBuilding = ((area2TexturePixelRatio>= minArea2TextureRatio) || (sum(H>=minAggregateBinHeight2) > 1 && ret(1,4) >= minMax2MeanRatio2) || (sum(H>=minAggregateBinHeight) > 1 && ret(1,8) <= minAggregateBinHeight2) || (ret(1,2) >= minAggregateBinHeight && ret(1,7) <= maxTotal) || (ret(1,4) >= minMax2MeanRatio && ret(1,2) >= minBinHeight) || (retInd(1,2) >= minBinHeight && ret(1,7) <= maxTotal));
        elseif (size(ret,1) == 0 && size(retInd,1)>0)
            trueBuilding = (retInd(1,2) >= minBinHeight);
        end
        %}
        %{
        if (size(ret,1)>0 && size(retInd,1)==0)
            trueBuilding = ((area2TexturePixelRatio>= minArea2TextureRatio) || (perMaxima >= minMax2TextureRatio && ret(1,7) <= maxTotal2) || (sum(H>=minAggregateBinHeight2) > 1 && ret(1,4) >= minMax2MeanRatio2) || (sum(H>=minAggregateBinHeight) > 1) || (ret(1,2) >= minAggregateBinHeight && ret(1,7) <= maxTotal) || (ret(1,4) >= minMax2MeanRatio && ret(1,2) >= minBinHeight));
        elseif (size(ret,1)>0 && size(retInd,1)>0)
            trueBuilding = ((area2TexturePixelRatio>= minArea2TextureRatio) || (perMaxima >= minMax2TextureRatio && ret(1,7) <= maxTotal2) || (sum(H>=minAggregateBinHeight2) > 1 && ret(1,4) >= minMax2MeanRatio2) || (sum(H>=minAggregateBinHeight) > 1) || (ret(1,2) >= minAggregateBinHeight && ret(1,7) <= maxTotal) || (ret(1,4) >= minMax2MeanRatio && ret(1,2) >= minBinHeight) || (retInd(1,2) >= minBinHeight && ret(1,7) <= maxTotal));
        elseif (size(ret,1) == 0 && size(retInd,1)>0)
            trueBuilding = (retInd(1,2) >= minBinHeight);
        end
        %}
        %{
        if (size(ret,1)>0 && size(retInd,1)==0)
            trueBuilding = ((sum(H>=minAggregateBinHeight2) > 1 && ret(1,4) >= minMax2MeanRatio2) || (sum(H>=minAggregateBinHeight) > 1) || (ret(1,2) >= minAggregateBinHeight && ret(1,7) <= maxTotal) || (ret(1,4) >= minMax2MeanRatio && ret(1,2) >= minBinHeight));
        elseif (size(ret,1)>0 && size(retInd,1)>0)
            trueBuilding = ((sum(H>=minAggregateBinHeight2) > 1 && ret(1,4) >= minMax2MeanRatio2) || (sum(H>=minAggregateBinHeight) > 1) || (ret(1,2) >= minAggregateBinHeight && ret(1,7) <= maxTotal) || (ret(1,4) >= minMax2MeanRatio && ret(1,2) >= minBinHeight) || (retInd(1,2) >= minBinHeight && ret(1,7) <= maxTotal));
        elseif (size(ret,1) == 0 && size(retInd,1)>0)
            trueBuilding = (retInd(1,2) >= minBinHeight);
        end
        %}
        if (trueBuilding)
            
            plot([P1(1,1) P2(1,1)], [P1(1,2) P2(1,2)],'-b','LineWidth', 2); hold on;
            plot([P2(1,1) P3(1,1)], [P2(1,2) P3(1,2)],'-b','LineWidth', 2); hold on;
            plot([P3(1,1) P4(1,1)], [P3(1,2) P4(1,2)],'-b','LineWidth', 2); hold on;
            plot([P4(1,1) P1(1,1)], [P4(1,2) P1(1,2)],'-b','LineWidth', 2); hold on;
            
            fprintf(fp,'%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %d\n', Buildings(i,1), Buildings(i,2), Buildings(i,3), Buildings(i,4), Buildings(i,5), Buildings(i,6), Buildings(i,7), Buildings(i,8), Buildings(i,9));
        else
            %{
            plot([P1(1,1) P2(1,1)], [P1(1,2) P2(1,2)],'-r','LineWidth', 2); hold on;
            plot([P2(1,1) P3(1,1)], [P2(1,2) P3(1,2)],'-r','LineWidth', 2); hold on;
            plot([P3(1,1) P4(1,1)], [P3(1,2) P4(1,2)],'-r','LineWidth', 2); hold on;
            plot([P4(1,1) P1(1,1)], [P4(1,2) P1(1,2)],'-r','LineWidth', 2); hold on;
            %}
        end
        %text(CP(1,1), CP(1,2), strcat(strcat(int2str(nB),','),int2str(ret*100)), 'Color', 'c');
        %text(CP(1,1), CP(1,2), int2str(nB), 'Color', 'c');
        
         %N_abv = find_ALS(A,R,N,P1,P2,P3,P4);
         %ma = mean(N_abv(:,3));
        %fprintf(fid,'%d\t%d\t%d\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%d\t%4.2f\n',nB,ret(1,1),ret(1,2),ret(1,3),ret(1,4),ret(1,5),ret(1,6),ret(1,7),ma);
        %fprintf(fid,'%d\t%4.2f\t%4.2f\t%d\n',i,area2TexturePixelRatio,perMaxima,ret(1,7));
        %{
        if (trueBuilding)
            [mx id] = max(ret(:,2));
            r = ret(id,:);
            fprintf(fid,'%d\t%d\t%d\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%d\n',nB,r(:,1),r(:,2),r(:,3),r(1,4),r(1,5),r(1,6),r(1,7));
        else
            fprintf(fid,'%d\n',nB);
        end
       
        Hist = [Hist;H];
         %}
        %{
        N_abv = find_ALS(A,R,N,P1,P2,P3,P4);
        ma = mean(N_abv(:,3));
        
        T_abv = find_ALS(A,R,T,P1,P2,P3,P4);
        id_edge = (T_abv(:,3) == 1);
        id_entr = (T_abv(:,4) == 1);
        m_edge = sum(id_edge)/size(T_abv,1);
        m_entr = sum(id_entr)/size(T_abv,1);
        
        G = find_ALS(A,R,grad,P1,P2,P3,P4);
        for j = 1:size(G,1)
            theta = 180*G(j,3)/pi;;
            %if (theta ~= -9999) % if defined
                %plot(G(j,1),G(j,2),'+');
                theta = theta + 90; % range: -90 to 90 is transferred to 0 to 180
                index = ceil(theta/binDist);
                if (index <= 0)
                    index = 1;
                elseif (index > 180/binDist)
                    index = 180/binDist;
                end
                Hist(nB,index) = Hist(nB,index)+1;
            %end
        end
        fprintf(fid,'%d\t%f\t%f\t%f\t',nB,ma,m_edge,m_entr);
        for k = 1:180/binDist
            fprintf(fid,'%d\t',Hist(nB,k));
        end
        fprintf(fid,'\n');
        %}
        nB= nB+1;
    %end
end
title('Refined building candidates');
hold off;
%fclose(fid);
fclose(fp);
%hold off;