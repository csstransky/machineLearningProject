function extendLines(A, R, bbox, mask, Rm, ndvi, ndviThresh, eMask, useEntropy, fnamein, fnameout)
%function detect_buildings02_mask05_04_02_hb_texture();
% add paths to the Matlab path list
%path(path,'C:\Awrangjeb\Drive I\Sydney-ALS\Digital Ortho Photo');
%path(path,'C:\Awrangjeb\Drive I\Sydney-ALS\LAS Files');
%path(path,'C:\Develop\Knoch');
%path(path,'C:\Awrangjeb\Drive I\Hobert');

% read orhto image
%[A, R, bbox] = geotiffread('524000e_5251000n_LB_02.tif');
extra = 0;
ImgRect = [bbox(1,1)-extra,bbox(1,2)-extra; bbox(1,1)-extra,bbox(2,2)+extra; bbox(2,1)+extra,bbox(2,2)+extra; bbox(2,1)+extra,bbox(1,2)-extra];

%read Line data
fp = fopen(fnamein, 'r');                  
L = fscanf(fp, '%f %f %f %f %f %f', [10 inf]);
fclose(fp);
L = L';

%{
%read Mask data
fp = fopen('mask_ground_524000e_5251000n_LB_02_dem1_pixels.txt', 'r');
size1 = fscanf(fp, '%f %f %f', [1 3]);
M = fscanf(fp, '%f %f %f', [3 inf]);
fclose(fp);
M = M';

%maskg = imread('mask_ground_crop7.jpg');
%figure; mapshow(maskg,R); hold on;
%for i=1:size(L,1)
%    P1 = L(i,1:2);
%    P2 = L(i,3:4);
%    plot([P1(1,1) P2(1,1)], [P1(1,2) P2(1,2)], '-c');
%end
%hold off;

%read NDVI_sigma data
fp = fopen('image_ndvi_524000e_5251000n_LB_02.txt', 'r');
size2 = fscanf(fp, '%f %f %f', [1 3]);
N = fscanf(fp, '%f %f %f', [3 inf]);
fclose(fp);
N = N';

%read texture info
fp = fopen('image_texture_524000e_5251000n_LB_02.txt','r');
T = fscanf(fp, '%f %f %d %d', [4 inf]);
fclose(fp);
T = T';

%read mask
mask = imread('mask_ground_524000e_5251000n_LB_02_dem1_pixels.jpg');
%}

% filter using ALS
dAl = 3; % radius of circle
%htThresh = 9.5;%7.5;
edgeThresh = 0.70; %nBin = 16; 
%Sin = 0.50; % Spacing_in_flight
%Sacross = 0.50; % Spacing_across_flight
%figure; mapshow(A(:,:,1:3),R); hold on; 
%ndviThresh = 10; 
entropyThresh = 0.30; %edgeShadowThresh = 0.05;

[srtLen idx] = sort(L(:,7),'descend');
 %seg_num = 0; Segs = []; dAc = 6; dAl = 3;
Flag = 1;
for i = 1:size(L,1)
    mat = L(idx(i,1),:);
    C1 = mat(1,1:2);
    C2 = mat(1,3:4);
    m = mat(1,5);
    c = mat(1,6);
    matP = mat(1,8:9);
    FlagT = L(idx(i,1),10);
    %plot([C1(1,1) C2(1,1)], [C1(1,2) C2(1,2)],'-g'); hold on;
     %plot(matP(1,1),matP(1,2),'*c'); hold on;
    %plot([C1(1,1) C2(1,1)], [C1(1,2) C2(1,2)],'-g','linewidth',1); hold on;
    
    [P1 P2 P3 P4] = find_rectVertices(m,c,dAl,C1,C2);

    Y = m*matP(1,1) + c;
    Y_diff = Y-matP(1,2);
                    
    Yp1 = m*P1(1,1) + c;
    Y_diffp1 = Yp1-P1(1,2);
                    
    if (sign(Y_diff) == sign(Y_diffp1))
        P_avb1 = P1;
        P_avb2 = P3;                        
    else
        P_avb1 = P2;
        P_avb2 = P4;                        
    end
                    
    [C1 P_avb1] = extend_line_segment(A,R,m,c,dAl,C1,P_avb1,C2,mask,Rm,ndvi,eMask,edgeThresh,ndviThresh, entropyThresh, ImgRect,Flag,FlagT, useEntropy); % extend left
    [C2 P_avb2] = extend_line_segment(A,R,m,c,dAl,C2,P_avb2,C1,mask,Rm,ndvi,eMask,edgeThresh,ndviThresh, entropyThresh, ImgRect,Flag,FlagT, useEntropy); %extend right
    

    
    %convert geographic to image
                    %for C1
                 c1y = floor((C1(1,1)*R(1,2) - R(3,1)*R(1,2) + R(3,2)*R(1,1) - C1(1,2)*R(1,1))/(R(2,1)*R(1,2)-R(2,2)*R(1,1)));
                 c1x = floor((C1(1,2) - R(3,2) - c1y*R(2,2))/R(1,2));
                    %for C2
                 c2y = floor((C2(1,1)*R(1,2) - R(3,1)*R(1,2) + R(3,2)*R(1,1) - C2(1,2)*R(1,1))/(R(2,1)*R(1,2)-R(2,2)*R(1,1)));
                 c2x = floor((C2(1,2) - R(3,2) - c2y*R(2,2))/R(1,2));
     
     %update
     mpxc = (C1(1,1)+C2(1,1))/2;
     mpyc = (C1(1,2)+C2(1,2))/2;
     mat(1,1:2) = C1;
     mat(1,3:4) = C2;
     %mat(1,7) = sqrt((c1x-c2x)*(c1x-c2x) + (c1y-c2y)*(c1y-c2y));
     mat(1,8:9) = find_midAbovePoint([mpxc mpyc],m,c,matP);
     mat(1,10:11) = [mpxc mpyc];
     L1(i,:) = [mat FlagT];
     
     %plot([C1(1,1) C2(1,1)], [C1(1,2) C2(1,2)],'-c','linewidth',2); hold on;
     %plot(mat(1,8),mat(1,9),'ow'); hold on;
     %plot(matP(1,1),matP(1,2),'ow'); hold on;
     %out = i
end
%hold off;

fp = fopen(fnameout,'w');
for i = 1:size(L1,1)
    fprintf(fp,'%f\t %f\t %f\t %f\t %f\t %f \t %f \t %f\t %f \t %f \t %f \t %d\n',L1(i,1), L1(i,2), L1(i,3),L1(i,4),L1(i,5), L1(i,6), L1(i,7),L1(i,8), L1(i,9), L1(i,10),L1(i,11),L1(i,12));
end
fclose(fp);
here = 1;