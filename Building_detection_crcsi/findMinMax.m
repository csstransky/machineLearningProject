function [minX, minY, maxX, maxY] = findMinMax(P1, P2, P3, P4, wd, ht)

	 % find minX and maxX
	if (P1(1,1) <= P2(1,1))
	
		minX = floor(P1(1,1));
		maxX = ceil(P2(1,1));
	
	else
	
		minX = floor(P2(1,1));
		maxX = ceil(P1(1,1));
    end
    
	if (P3(1,1) < minX)
		minX = floor(P3(1,1));
    elseif (P3(1,1) > maxX)
		maxX = ceil(P3(1,1));
    end
    
	if (P4(1,1) < minX)
		minX = floor(P4(1,1));
    elseif (P4(1,1) > maxX)
		maxX = ceil(P4(1,1));
    end

	% find minY and maxY
	if (P1(1,2) <= P2(1,2))
	
		minY = floor(P1(1,2));
		maxY = ceil(P2(1,2));
	
	else
	
		minY = floor(P2(1,2));
		maxY = ceil(P1(1,2));
    end
    
	if (P3(1,2) < minY)
		minY = floor(P3(1,2));
    elseif (P3(1,2) > maxY)
		maxY = ceil(P3(1,2));
    end
    
	if (P4(1,2) < minY)
		minY = floor(P4(1,2));
    elseif (P4(1,2) > maxY)
		maxY = ceil(P4(1,2));
    end

	if (minX  <= 0)
		minX = 1;
    end
	if (minY <= 0)
		minY = 1;	
    end
	if (maxX <= 0)	
		maxX = 1;
    end
	if (maxY <= 0)
		maxY = 1;	
    end
	if (minX > ht)
		minX = ht;	
    end
	if (maxX > ht)
		maxX = ht;
    end
	if (minY > wd)
		minY = wd;
    end
	if (maxY > wd)
		maxY = wd;
    end
 