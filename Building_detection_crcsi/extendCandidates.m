function extendCandidates(A, R, pMask, sMask, Rm, useReduction, fnamein, fnameout)
%path(path,'C:\Awrangjeb\Drive I\Sydney-ALS\Digital Ortho Photo');
%path(path,'C:\Awrangjeb\Drive I\Sydney-ALS\LAS Files');
%path(path,'C:\Develop\Knoch');
%path(path,'C:\Awrangjeb\Drive I\Hobert');

% read orhto image
%[A, R, bbox] = geotiffread('524000e_5251000n_LB_02.tif');
%[sx sy sz] = size(A);
A = A(:,:,1:3);
YIQ = rgb2ntsc(A); H = YIQ(:,:,2); S = YIQ(:,:,3); I = YIQ(:,:,1);% I = uint8(floor((H + 0.6)*255/1.2)); 
%smooth
fH = fspecial('gaussian');
H = imfilter(H,fH,'replicate');
S = imfilter(S,fH,'replicate');
I = imfilter(I,fH,'replicate');
yiq3b(:,:,1) = H;
yiq3b(:,:,2) = S;
yiq3b(:,:,3) = I;

%mapshow(A, R); hold on; plot3(ALS(:,1), ALS(:,2), ALS(:,3),'.');
%Img = read_imagePixels_tri('image_hue_sat_int_ori_524000e_5251000n_LB_02.txt');
%Img = Img';
%A = decorrstretch(A1,'Tol',0.01);% figure; imshow(C)

%read ALS data
%fp = fopen('last_crop7dem.txt', 'r');
%ALS = fscanf(fp, '%f %f %f', [3 inf]);
%fclose(fp);
%ALS = ALS';

%{
%read up mask data
fp = fopen('mask_up_524000e_5251000n_LB_02_dem1_pixels.txt', 'r');
size2 = fscanf(fp, '%f %f %f', [1 3]);
Mup = fscanf(fp, '%f %f %f', [3 inf]);
fclose(fp);
Mup = Mup';

%read up mask data
fp = fopen('mask_ground_524000e_5251000n_LB_02_dem1_pixels.txt', 'r');
size2 = fscanf(fp, '%f %f %f', [1 3]);
Mg = fscanf(fp, '%f %f %f', [3 inf]);
fclose(fp);
Mg = Mg';
%}
%read image gradient for image 174_184
%{
fid = fopen('gradient.txt','r');
grad = fscanf(fid,'%d\t%d\t%f',[3 inf]);
fclose(fid);
grad = grad';
x = sx/2 - grad(:,2);
y = sy/2 + grad(:,1);
O = ones(size(x,1),1);
grad(:,1:2) = [x y O]*R;
%}

%read NDVI_sigma data
%fp = fopen('image_ndvi_sigma_crop7dem.txt', 'r');
%size2 = fscanf(fp, '%f %f %f', [1 3]);
%N = fscanf(fp, '%f %f %f', [3 inf]);
%fclose(fp);
%N = N';

%fp = fopen('buildings_crop7dem_up.txt','r');
%Bup = fscanf(fp, '%f %f %f %f %f %f %f %f', [8 inf]);
%Bup = Bup';
%fclose(fp);

fp = fopen(fnamein,'r');
Bground = fscanf(fp, '%f %f %f %f %f %f %f %f', [9 inf]);
Bground = Bground';
fclose(fp);

%fp = fopen('Segs_hue_crop7dem_img.txt','r');
%Segs = fscanf(fp, '%f %f %f %f %f %f ', [6 inf]);
%Segs = Segs';
%fclose(fp);

max_buildingLength = 50;
%[Bground CPoint] = remove_overlaps_ground(A,R,Bground,max_buildingLength);

%show straight line segments

%Ig = rgb2gray(A);
%figure; mapshow(Ig,R); hold on; 
%for i = 1:size(Bup,1)
%    P1 = Bup(i,1:2); P2 = Bup(i,3:4); P3 = Bup(i,5:6); P4 = Bup(i,7:8);
%    plot([P1(1,1) P2(1,1)], [P1(1,2) P2(1,2)],'-r','LineWidth', 2); hold on;
%    plot([P2(1,1) P3(1,1)], [P2(1,2) P3(1,2)],'-r','LineWidth', 2); hold on;
%    plot([P3(1,1) P4(1,1)], [P3(1,2) P4(1,2)],'-r','LineWidth', 2); hold on;
%    plot([P4(1,1) P1(1,1)], [P4(1,2) P1(1,2)],'-r','LineWidth', 2); hold on;
%    ImgInfo{i} = find_ALS(A,R,Img,P1,P2,P3,P4);
%end
%for i = 1:size(Bground,1)
    %P1 = Bground(i,1:2); P2 = Bground(i,3:4); P3 = Bground(i,5:6); P4 = Bground(i,7:8);
    %plot([P1(1,1) P2(1,1)], [P1(1,2) P2(1,2)],'-k','LineWidth', 3); hold on;
    %plot([P2(1,1) P3(1,1)], [P2(1,2) P3(1,2)],'-k','LineWidth', 3); hold on;
    %plot([P3(1,1) P4(1,1)], [P3(1,2) P4(1,2)],'-k','LineWidth', 3); hold on;
    %plot([P4(1,1) P1(1,1)], [P4(1,2) P1(1,2)],'-k','LineWidth', 3); hold on;
    
    %CP =  [(P1(1,1) + P3(1,1))/2 (P1(1,2) + P3(1,2))/2];
    %text(CP(1,1),CP(1,2),int2str(i), 'Color', 'cyan');
%end
%for i = 1:size(Segs,1)
%    P1 = Segs(i,1:2); P2 = Segs(i,3:4);
%    plot([P1(1,1) P2(1,1)], [P1(1,2) P2(1,2)],'-g','LineWidth', 2); hold on;    
%end
%hold off;

%calculate color histograms / slopes / centerpoints for ground buildings
per = 0.15;%0.30;
for i = 1:size(Bground,1)
    P1 = Bground(i,1:2); P2 = Bground(i,3:4); P3 = Bground(i,5:6); P4 = Bground(i,7:8);
    
    %centerpoint
    BGroundInfo(i,1:2) = [(P1(1,1) + P3(1,1))/2 (P1(1,2) + P3(1,2))/2];
    % m, c for P1P2 or P3P4
    BGroundInfo(i,3) = (P2(1,2) - P1(1,2))/(P2(1,1) - P1(1,1)); % m12 or m34
    BGroundInfo(i,4) = P1(1,2) - BGroundInfo(i,3)*P1(1,1); %c12
    BGroundInfo(i,5) = P3(1,2) - BGroundInfo(i,3)*P3(1,1); %c34
    % m, c for P2P3 or P4P1
    BGroundInfo(i,6) = (P2(1,2) - P3(1,2))/(P2(1,1) - P3(1,1)); % m23 or m41
    BGroundInfo(i,7) = P3(1,2) - BGroundInfo(i,6)*P3(1,1); %c23
    BGroundInfo(i,8) = P4(1,2) - BGroundInfo(i,6)*P4(1,1);%c41
    
    BGroundInfo(i,9:10) = P1;
    BGroundInfo(i,11:12) = P2;
    BGroundInfo(i,13:14) = P3;
    BGroundInfo(i,15:16) = P4;  
     
    if useReduction
    crop12 = per*sqrt((P1(1,1)-P2(1,1))*(P1(1,1)-P2(1,1)) + (P1(1,2)-P2(1,2))*(P1(1,2)-P2(1,2))); %crop34
    crop23 = per*sqrt((P3(1,1)-P2(1,1))*(P3(1,1)-P2(1,1)) + (P3(1,2)-P2(1,2))*(P3(1,2)-P2(1,2))); %crop41
    
    Cp = BGroundInfo(i,1:2);
    %plot(Cp(1,1), Cp(1,2), 'or');
    m12 = BGroundInfo(i,3); %m34
    c12 = BGroundInfo(i,4);
    c34 = BGroundInfo(i,5);
    m23 = BGroundInfo(i,6); %m41
    c23 = BGroundInfo(i,7);
    c41 = BGroundInfo(i,8);
    
     [Q1 Q2 Q3 Q4] = find_rectVertices(m12,c12,crop23,P1,P2);     
     Y = m12*[Cp(1,1); Q1(1,1)] + c12;
     Y_diff = Y-[Cp(1,2); Q1(1,2)];
     if ((Y_diff(1,1) <= 0 && Y_diff(2,1) <= 0) || (Y_diff(1,1)  > 0 && Y_diff(2,1) > 0)) % Cp & Q1 are on same side
        Q12L = Q1;
        Q12R = Q3;        
     else
        Q12L = Q2;
        Q12R = Q4;
     end                
     %plot(Q12L(1,1), Q12L(1,2), '+g');
     %plot(Q12R(1,1), Q12R(1,2), '*g');
     %plot([Q12L(1,1) Q12R(1,1)], [Q12L(1,2) Q12R(1,2)], '-g');
     
     [Q1 Q2 Q3 Q4] = find_rectVertices(m23,c23,crop12,P2,P3);     
%plot([Q1(1,1) Q2(1,1)], [Q1(1,2) Q2(1,2)],'-g'); hold on; 
%plot([Q2(1,1) Q4(1,1)], [Q2(1,2) Q4(1,2)],'-b'); hold on;
%plot([Q4(1,1) Q3(1,1)], [Q4(1,2) Q3(1,2)],'-y'); hold on;
%plot([Q1(1,1) Q3(1,1)], [Q1(1,2) Q3(1,2)],'-r'); hold off;

     Y = m23*[Cp(1,1); Q1(1,1)] + c23;
     Y_diff = Y-[Cp(1,2); Q1(1,2)];
     if ((Y_diff(1,1) <= 0 && Y_diff(2,1) <= 0) || (Y_diff(1,1)  > 0 && Y_diff(2,1) > 0)) % Cp & Q1 are on same side
        Q23U = Q1;
        Q23D = Q3;
     else
        Q23U = Q2;
        Q23D = Q4;
     end                
     %plot(Q23U(1,1), Q23U(1,2), '+b');
     %plot(Q23D(1,1), Q23D(1,2), '*b');
     %plot([Q23U(1,1) Q23D(1,1)], [Q23U(1,2) Q23D(1,2)], '-b');

     
     [Q1 Q2 Q3 Q4] = find_rectVertices(m12,c34,crop23,P4,P3);     
     Y = m12*[Cp(1,1); Q1(1,1)] + c34;
     Y_diff = Y-[Cp(1,2); Q1(1,2)];
     if ((Y_diff(1,1) <= 0 && Y_diff(2,1) <= 0) || (Y_diff(1,1)  > 0 && Y_diff(2,1) > 0)) % Cp & Q1 are on same side
        Q34L = Q1;
        Q34R = Q3;        
     else
        Q34L = Q2;
        Q34R = Q4;
     end                
     %plot(Q34L(1,1), Q34L(1,2), '+y');
     %plot(Q34R(1,1), Q34R(1,2), '*y');
      %plot([Q34L(1,1) Q34R(1,1)], [Q34L(1,2) Q34R(1,2)], '-y');
      
     [Q1 Q2 Q3 Q4] = find_rectVertices(m23,c41,crop12,P1,P4);     
     Y = m23*[Cp(1,1); Q1(1,1)] + c41;
     Y_diff = Y-[Cp(1,2); Q1(1,2)];
     if ((Y_diff(1,1) <= 0 && Y_diff(2,1) <= 0) || (Y_diff(1,1)  > 0 && Y_diff(2,1) > 0)) % Cp & Q1 are on same side
        Q41U = Q1;
        Q41D = Q3;
     else
        Q41U = Q2;
        Q41D = Q4;
     end                
     %plot(Q41U(1,1), Q41U(1,2), '+m');
     %plot(Q41D(1,1), Q41D(1,2), '*m');
     %plot([Q41U(1,1) Q41D(1,1)], [Q41U(1,2) Q41D(1,2)], '-m');
     
     c12 = Q12L(1,2) - m12*Q12L(1,1);
     c23 = Q23U(1,2) - m23*Q23U(1,1);
     c34 = Q34L(1,2) - m12*Q34L(1,1);
     c41 = Q41U(1,2) - m23*Q41U(1,1);
     BGroundInfo(i,4) = c12;
     BGroundInfo(i,5) = c34;    
     BGroundInfo(i,7) = c23;
     BGroundInfo(i,8) = c41;
     
     P1(1,1) = (c41-c12)/(m12-m23);
     P1(1,2) = m12*P1(1,1) + c12;
     %plot(P1(1,1), P1(1,2), 'ob');

     P2(1,1) = (c23-c12)/(m12-m23);
     P2(1,2) = m12*P2(1,1) + c12;
     %plot(P2(1,1), P2(1,2), 'ob');

     P3(1,1) = (c23-c34)/(m12-m23);
     P3(1,2) = m12*P3(1,1) + c34;
     %plot(P3(1,1), P3(1,2), 'ob');

     P4(1,1) = (c41-c34)/(m12-m23);
     P4(1,2) = m12*P4(1,1) + c34;
     %plot(P4(1,1), P4(1,2), 'ob');

     BGroundInfo(i,9:10) = P1;
     BGroundInfo(i,11:12) = P2;
     BGroundInfo(i,13:14) = P3;
     BGroundInfo(i,15:16) = P4;  
    %else%if useReduction
    end%if useReduction
end
%hold off;

%maskg = imread('mask_ground_524000e_5251000n_LB_02_dem1_pixels.jpg');
%masku = imread('mask_up_524000e_5251000n_LB_02_dem1_pixels.jpg');
figure; mapshow(A,R); hold on; 
dAl = 0.35; %ndviThresh = 35; 
edgeThresh = 0.40; %Buildings = [];
UpedgeThresh = 0.9;
UpedgeThreshPool = 0.7;
%YIQ = rgb2ntsc(A); H = YIQ(:,:,2); S = YIQ(:,:,3); I = YIQ(:,:,1);% I = uint8(floor((H + 0.6)*255/1.2)); 
%nDetected = 0; 
Buildings = []; n = 0;

for i = 1:size(Bground,1)        
        Flag = 1;
        %Cp = BGroundInfo(i,1:2);
        m12 = BGroundInfo(i,3); %m34
        %c12 = BGroundInfo(i,4);
        %c34 = BGroundInfo(i,5);
        m23 = BGroundInfo(i,6); %m41    
        %c23 = BGroundInfo(i,7);
        %c41 = BGroundInfo(i,8);
        P1 = BGroundInfo(i,9:10); P2 = BGroundInfo(i,11:12); P3 = BGroundInfo(i,13:14); P4 = BGroundInfo(i,15:16);  
        
        
        %len1 = sqrt((P1(1,1)-P4(1,1)) * (P1(1,1)-P4(1,1)) + (P1(1,2)-P4(1,2))*(P1(1,2)-P4(1,2)));
        %len2 = sqrt((P1(1,1)-P2(1,1)) * (P1(1,1)-P2(1,1)) + (P1(1,2)-P2(1,2))*(P1(1,2)-P2(1,2)));
        
        % if (i == 89)
        %plot([P1(1,1) P2(1,1)], [P1(1,2) P2(1,2)],'-b','LineWidth', 1); hold on;
        %plot([P2(1,1) P3(1,1)], [P2(1,2) P3(1,2)],'-b','LineWidth', 1); hold on;
        %plot([P3(1,1) P4(1,1)], [P3(1,2) P4(1,2)],'-b','LineWidth', 1); hold on;
        %plot([P4(1,1) P1(1,1)], [P4(1,2) P1(1,2)],'-b','LineWidth', 1); hold on;
        % end
        
        ids = check_overlaps(A,R,Buildings,max_buildingLength, P1, P2, P3, P4);
        
        %Array = [];
        if (size(ids,1)>0)
            for j = 1:size(ids,1)
                if (ids(j,2) == 4) % all 4 points are in
                    Flag = 0;
                    break;
                end
            end
        end
        
        %test if it is a pool, when the secondary mask will be dominantly
        %white within the initial building area
        P1m = obj2obs(P1,Rm);
        P2m = obj2obs(P2,Rm);
        P3m = obj2obs(P3,Rm);
        P4m = obj2obs(P4,Rm);
        [totalsP numValsP] = findImagePointCounts(sMask, P1m, P2m, P3m, P4m, 0);
        MupInfo = numValsP/totalsP;
        if MupInfo < UpedgeThreshPool
            Flag = 0;
        end
        
        if (Flag == 1)
            P1i = obj2obs(P1,R);
            P2i = obj2obs(P2,R);
            P3i = obj2obs(P3,R);
            P4i = obj2obs(P4,R);
            
            %HueSatIn = find_ALS(I,R,Img,P1,P2,P3,P4);            
            hsi = findImagePoints(yiq3b, P1i, P2i, P3i, P4i);
            h = hsi(:,1);
            s = hsi(:,2);
            int = hsi(:,3);
            %h = findImagePoints(H, P1, P2, P3, P4);
            %s = findImagePoints(S, P1, P2, P3, P4);
            %int = findImagePoints(I, P1, P2, P3, P4);
            
            %h = HueSatIn(:,3);
            %s = HueSatIn(:,4);
            %int = HueSatIn(:,5);
            if (size(hsi,1)>0)
                [HThresh Hhist] = find_histogram(h');
                [SThresh Shist] = find_histogram(s');
                [IThresh Ihist] = find_histogram(int');

                hThresh = find_Thresh(HThresh, Hhist);
                sThresh = find_Thresh(SThresh, Shist);
                iThresh = find_Thresh(IThresh, Ihist);
                %meanHue = mean(HueIn(:,3));
                %meanHue = 0.1;
                %{
                HueBinary = logical(zeros(size(H)));
                SatBinary = logical(zeros(size(S)));
                IntBinary = logical(zeros(size(I)));
                for j=1:size(hThresh,1)
                    BinH = (H >= hThresh(j,1)) & (H < hThresh(j,2)); 
                    HueBinary = HueBinary | BinH;
                 end

                for j=1:size(sThresh,1)
                    BinS = (S >= sThresh(j,1)) & (S < sThresh(j,2)); 
                    SatBinary = SatBinary | BinS;                
                end

                 for j=1:size(iThresh,1)
                    BinI = (I >= iThresh(j,1)) & (I < iThresh(j,2)); 
                    IntBinary = IntBinary | BinI;                
                end
                %}
               %HueBinary = H >= meanHue; 
               %figure; imshow(HueBinary); 
               %figure; imshow(SatBinary);
               %figure; mapshow(HueBinary,R); hold on; 

               %check in which Bups
                %[ret bupID] = testInBuidlings_ID(Bup,Cp); 
                %Imga = ImgInfo{bupID};

                %[P2 P3] = extend_segment_binary_image(HueBinary,R,m12,dAc,dAl,P2,P3,P1,Imga,edgeThresh,Sin, Sacross);
               %[P1 P2] = extend_segment_binary_image_thresh_huesatint(A,R,m23,dAl,P1,P2,P4,Img,edgeThresh,hThresh,sThresh,iThresh);
               %[P2 P3] = extend_segment_binary_image_thresh_huesatint(A,R,m12,dAl,P2,P3,P1,Img,edgeThresh,hThresh,sThresh,iThresh);
               %[P3 P4] = extend_segment_binary_image_thresh_huesatint(A,R,m23,dAl,P3,P4,P2,Img,edgeThresh,hThresh,sThresh,iThresh);
               %[P4 P1] = extend_segment_binary_image_thresh_huesatint(A,R,m12,dAl,P4,P1,P3,Img,edgeThresh,hThresh,sThresh,iThresh);
                %if (i == 6)
                %    here = 1;
                %end
               [P1_1 P2_1] = extend_building_side(A,R,m23,dAl,P1,P2,P4,yiq3b,pMask,sMask,Rm,UpedgeThresh,edgeThresh,hThresh,sThresh,iThresh);
               [P2_2 P3_1] = extend_building_side(A,R,m12,dAl,P2,P3,P1,yiq3b,pMask,sMask,Rm,UpedgeThresh,edgeThresh,hThresh,sThresh,iThresh);
               [P3_2 P4_1] = extend_building_side(A,R,m23,dAl,P3,P4,P2,yiq3b,pMask,sMask,Rm,UpedgeThresh,edgeThresh,hThresh,sThresh,iThresh);
               [P4_2 P1_2] = extend_building_side(A,R,m12,dAl,P4,P1,P3,yiq3b,pMask,sMask,Rm,UpedgeThresh,edgeThresh,hThresh,sThresh,iThresh);

               [P1 P2 P3 P4] = find_extended_rectangle(m12, m23, P1_1, P2_2, P3_2, P4_2);
                   
               %Len1 = sqrt((P1(1,1)-P4(1,1))*(P1(1,1)-P4(1,1)) + (P1(1,2)-P4(1,2))*(P1(1,2)-P4(1,2)));
               %Len2 = sqrt((P1(1,1)-P2(1,1))*(P1(1,1)-P2(1,1)) + (P1(1,2)-P2(1,2))*(P1(1,2)-P2(1,2)));

               %cut = 0;
               %if (Len1 > max_buildingLength || Len2 > max_buildingLength || Len1 > 3*len1 || Len2 > 3*len2)

              %      N_left = find_ALS(A,R,N,BGroundInfo(i,9:10),BGroundInfo(i,11:12),BGroundInfo(i,13:14),BGroundInfo(i,15:16));
              %      m = mean(N_left(:,3));
              %      if (m > ndviThresh)
              %          P1 = BGroundInfo(i,9:10); P2 = BGroundInfo(i,11:12); P3 = BGroundInfo(i,13:14); P4 = BGroundInfo(i,15:16);  
              %          cut = 1;
              %      end
              %else
                   %it is okay;
              % end
             %{
              if (Bground(i,9))
                  %validate
                  ret = valid_building(A,R,grad,P1,P2,P3,P4)
                  if (ret)
                    n = n+1;
                    Buildings(n,:) = [P1 P2 P3 P4 Bground(i,9)];% i Len1 Len2];
                  
                    plot([P1(1,1) P2(1,1)], [P1(1,2) P2(1,2)],'-c','LineWidth', 2); hold on;
                    plot([P2(1,1) P3(1,1)], [P2(1,2) P3(1,2)],'-c','LineWidth', 2); hold on;
                    plot([P3(1,1) P4(1,1)], [P3(1,2) P4(1,2)],'-c','LineWidth', 2); hold on;
                    plot([P4(1,1) P1(1,1)], [P4(1,2) P1(1,2)],'-c','LineWidth', 2); hold on;
                  end
              else
                plot([P1(1,1) P2(1,1)], [P1(1,2) P2(1,2)],'-b','LineWidth', 2); hold on;
                plot([P2(1,1) P3(1,1)], [P2(1,2) P3(1,2)],'-b','LineWidth', 2); hold on;
                plot([P3(1,1) P4(1,1)], [P3(1,2) P4(1,2)],'-b','LineWidth', 2); hold on;
                plot([P4(1,1) P1(1,1)], [P4(1,2) P1(1,2)],'-b','LineWidth', 2); hold on;
                
                n = n+1;
                Buildings(n,:) = [P1 P2 P3 P4 Bground(i,9)];% i Len1 Len2];
              end
              %}
                %nDetected = nDetected+1;
                %BGroundInfo(nDetected,9:10) = P1;
                %BGroundInfo(nDetected,11:12) = P2;
                %BGroundInfo(nDetected,13:14) = P3;
                %BGroundInfo(nDetected,15:16) = P4;  
                plot([P1(1,1) P2(1,1)], [P1(1,2) P2(1,2)],'-b','LineWidth', 2); hold on;
                plot([P2(1,1) P3(1,1)], [P2(1,2) P3(1,2)],'-b','LineWidth', 2); hold on;
                plot([P3(1,1) P4(1,1)], [P3(1,2) P4(1,2)],'-b','LineWidth', 2); hold on;
                plot([P4(1,1) P1(1,1)], [P4(1,2) P1(1,2)],'-b','LineWidth', 2); hold on;
                n = n+1;
                Buildings(n,:) = [P1 P2 P3 P4 Bground(i,9)];% i Len1 Len2];
            end
            %out = [i n]
           % if (i == 41)
            %    here = 1;
            %end
        end%Flag
end

title('Building candidates');
hold off;

fp = fopen(fnameout,'w');
for i = 1:size(Buildings,1)
    fprintf(fp,'%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %d\n', Buildings(i,1), Buildings(i,2), Buildings(i,3), Buildings(i,4), Buildings(i,5), Buildings(i,6), Buildings(i,7), Buildings(i,8), Buildings(i,9));
end
fclose(fp);

here = 1;