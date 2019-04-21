function [ndvi ndvisig] = computeNDVI(A, show);

%path(path,'F:\Office Desktop - 22 June 2011\Awrangjeb\Drive I\Hobart\Imagery');

%[A, R, bbox] = geotiffread('524000e_5251000n_LB_02_sample.tif');

[height width bands] = size(A);

if bands > 3
    pidx1 = 1;
    pidx2 = 4;
elseif bands == 3
    pidx1 = 2;
	pidx2 = 1;
else
    return;
end

%get Gaussian noise
varP = 0.0;
varN = 0.0;

tileSize = 512;

%nTileV = 0;
%nTileH = 0;
%vn = [];
%vp = [];
i = 1; 
while (i < height)
    endx = i+tileSize-1;
    if (endx > height)
        endx = height;
    end
    %nTileV =  nTileV+1;
    
    j = 1;
    while (j < width)        
        endy = j+tileSize-1;
        if (endy > width)
            endy = width;
        end
        %nTileH =  nTileH+1;
        
        %if (nTileH == 10)
        %    here = 1;
        %end
        tile = A(i:endx, j:endy,:);
        sigmaN = getGaussianNoise(tile);
        
        sigN = sigmaN(1,pidx1);
		sigP = sigmaN(1,pidx2);
		if (sigN > varN) 
            varN = sigN;
        end
		if (sigP > varP) 
            varP = sigP;
        end
        j = j+tileSize;
        
        %vn = [vn sigN];
        %vp = [vp sigP];
    end
    i = i+tileSize;    
    
end
varN = varN*varN;
varP = varP*varP;

%{
fp = fopen('store.txt','w');
fprintf(fp,'%f\t%f',varN,varP);
fclose(fp);

clear all;

path(path,'F:\Office Desktop - 22 June 2011\Awrangjeb\Drive I\Hobart\Imagery');

[A, R, bbox] = geotiffread('524000e_5251000n_LB_02_sample.tif');

[height width bands] = size(A);
if bands > 3
    pidx1 = 1;
    pidx2 = 4;
elseif bands == 3
    pidx1 = 2;
	pidx2 = 1;
else
    return;
end

fp = fopen('store.txt','r');
fdata = fscanf(fp,'%f\t%f',[1 2]);
fclose(fp);
varN = fdata(1,1);
varP = fdata(1,2);
%}

%compute NDVI
FLT_EPSILON = 1.192092896e-07;
bandN = double(A(:,:,pidx1));
bandP = double(A(:,:,pidx2));
        
s = bandP + bandN;
vi = -101.0;
sigvi = -1.0;

%idx = s > FLT_EPSILON;
%sum = zeros(height,width);
sum = 1.0 ./ s;
vi = 100.0 * ((bandP - bandN) .* sum);
sigvi = 200.0 * sum .* sum .* sqrt(bandP .* bandP * varP + bandN .* bandN * varN);

idx = sigvi > 100.0;
sigvi(idx) = 100.0;
          
ndvi = vi;
ndvisig = sigvi;

if show
    
    mn = min(min(ndvi));
    mx = max(max(ndvi));
    figure; imshow(uint8((ndvi-mn)*255/(mx-mn))); title('NDVI image');
    
    mn = min(min(ndvisig));
    mx = max(max(ndvisig));
    figure; imshow(uint8((ndvisig-mn)*255/(mx-mn))); title('NDVI sigma image');
end


function sigmaN = getGaussianNoise(tile)

alpha1 = 0.15;
alpha2 = 0.35;

marginWidth = 0;
factor = 1;

[height width bands] = size(tile);

sigmaN = zeros(1,bands);
factorsq = factor * factor * 0.5;

histSize = 65536;

tileSizeX = height;
tileSizeY = width;

for b = 1:bands
    histogram = zeros(1,histSize);
    
    for i = 1:tileSizeX-1
        for j = 1:tileSizeY-1
            f00 = double(tile(i,j,b));
			f10 = double(tile(i,j+1,b));
			f01 = double(tile(i+1,j,b));
			f11 = double(tile(i+1,j+1,b));
            
            d1 = f11 - f00;
			d2 = f10 - f01;
    
			index = floor((d1 * d1 + d2 * d2) * factorsq) + 1;
			if (index <= histSize)
                histogram(1,index) = histogram(1,index)+1;
            end
        end
    end
    
    redSize = width * height - 2 * marginWidth * (width + height - 2* marginWidth);
  
	sigmaN(1,b) = getSigma(histogram, redSize, alpha1, alpha2) / factor;
end
    

function sigma = getSigma(histogram, numberOfPixels, alpha1, alpha2)
alphaCount = numberOfPixels * alpha1; 
count = 0;

i = 1;
while (count < alphaCount)
    count = count + histogram(1,i);	
    i = i+1;
end

alphaPoint = (i-1) - ((count - alphaCount)/histogram(1,i - 1));
theoretical = - 2.0 * log(1.0 - alpha1);
sigma = 2.0 * alphaPoint / theoretical;

histSize = size(histogram,2);

if (alpha2 > 0.0)
	if (sigma >= histSize) 
        sigma = histSize;
    end
	count = 0;
	for u = 1:(sigma + 1)
		count = count + histogram(1,u);
    end
	alphaCount = floor(count * alpha2 + 0.5);
        
	count = 0;
    u = 1;
	while (count < alphaCount)	 
        count = count + histogram(1,u);		
        u = u+1;
	end
        
	alphaPoint = (u-1) - ((count - alphaCount)/histogram(1,u - 1));
	theoretical =  - 2.0 * log(1 - (1 - exp(-1)) * alpha2);
	sigma = 2.0 * alphaPoint / theoretical;		    
end

sigma = sqrt(sigma/2);

    
    