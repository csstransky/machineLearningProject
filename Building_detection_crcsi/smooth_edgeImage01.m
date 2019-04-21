function [edgeImage curve_sm] = smooth_edgeImage01(BW,Gap_size,minLength);

% read orhto image
%[A, R, bbox] = geotiffread('crop4.tif');

%read ALS data
%fp = fopen('Test 3_crop4.txt', 'r');
%ALS = fscanf(fp, '%f %f %f', [3 inf]);
%fclose(fp);
%ALS = ALS';

%[ght mask] = find_gHeightMask1(ALS,A,R);
%BW = edge(mask,'canny'); 

% Extract curves from the edge-image
%Gap_size = 1; minLength = 20;
[curve,curve_start,curve_end,curve_mode,curve_num,TJ,img1] = extract_curve(BW, Gap_size, minLength);  

[Szx Szy] = size(BW);
BW_edge=zeros(Szx,Szy);
seg_size = 50; sig = 3;
[gau w] = find_Gaussian(sig);
for i = 1:curve_num
    X = curve{i}(:,1);
    Y = curve{i}(:,2);
    L = size(X,1);
    curve_sm{i} = []; st = 1;
    while L>(seg_size + 2*w)
        x = X(st:st+seg_size-1,1);
        y = Y(st:st+seg_size-1,1);
        
        x1=[ones(w,1)*2*x(1)-x(w+1:-1:2);x;ones(w,1)*2*x(seg_size)-x(seg_size-1:-1:seg_size-w)];
        y1=[ones(w,1)*2*y(1)-y(w+1:-1:2);y;ones(w,1)*2*y(seg_size)-y(seg_size-1:-1:seg_size-w)];
            
        Fsx = conv(x1,gau,'same'); %xx = Fsx(w+1:seg_size+w);
        Fsy = conv(y1,gau,'same'); %yy = Fsy(w+1:seg_size+w);
        curve_sm{i}(st:st+seg_size-1,:) = [Fsx(w+1:seg_size+w,1) Fsy(w+1:seg_size+w,1)];
        st = st+seg_size;
        L = L-seg_size;
    end
    if L > w
        x = X(st:end,1);
        y = Y(st:end,1);        
        sz = size(x,1);
        
        x1=[ones(w,1)*2*x(1)-x(w+1:-1:2);x;ones(w,1)*2*x(L)-x(L-1:-1:L-w)];
        y1=[ones(w,1)*2*y(1)-y(w+1:-1:2);y;ones(w,1)*2*y(L)-y(L-1:-1:L-w)];
            
        Fsx = conv(x1,gau,'same'); %xx = Fsx(w+1:seg_size+w);
        Fsy = conv(y1,gau,'same'); %yy = Fsy(w+1:seg_size+w);
        curve_sm{i}(st:st+sz-1,:) = [Fsx(w+1:sz+w,1) Fsy(w+1:sz+w,1)];
    else
        curve_sm{i}(st:st+sz-1,:) = [X(st:end,1) Y(st:end,1)];
    end
    %if (i == 13)
    %    here = 1;
   % end
    Pts = round(curve_sm{i});
    for j = 1:size(Pts,1)
        BW_edge(Pts(j,1),Pts(j,2))=1;
    end
   % i1 = i
end

edgeImage = BW_edge;