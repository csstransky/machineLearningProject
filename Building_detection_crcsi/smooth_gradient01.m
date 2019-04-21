function [smoothed_curve gard_curve c_intersect] = smooth_gradient01(curve, curve_mode, curve_num, sizex, sizey);
%function [smoothed_curve gard_curve c_intersect] = smooth_gradient01(curve, curve_mode, curve_num, R);

% grad = (y2-y1) / (x2-x1)
% c = yi - m.xi
s = 3.0;
[gau W] = find_Gaussian(s);
%figure;
for i = 1:curve_num    
    x = curve{i}(:,2) - sizey/2;
    y = sizex/2 - curve{i}(:,1);    
    %plot(x,y); hold on;
    %x = curve{i}(:,1);
    %y = curve{i}(:,2);    
    L = size(x,1);
    %o = ones(L,1);
    %M = [x y o]*R;
    %x = M(:,1);
    %y = M(:,2);
    if (L>W)    
        if curve_mode(i,:)=='loop'
            x1=[x(L-W+1:L);x;x(1:W)];
            y1=[y(L-W+1:L);y;y(1:W)];
        else
            x1=[ones(W,1)*2*x(1)-x(W+1:-1:2);x;ones(W,1)*2*x(L)-x(L-1:-1:L-W)];
            y1=[ones(W,1)*2*y(1)-y(W+1:-1:2);y;ones(W,1)*2*y(L)-y(L-1:-1:L-W)];
        end    
      
        xx=conv(x1,gau);
        xx=xx(W+1:L+3*W);
        yy=conv(y1,gau);
        yy=yy(W+1:L+3*W);

        Xu=[xx(2)-xx(1) ; (xx(3:L+2*W)-xx(1:L+2*W-2))/2 ; xx(L+2*W)-xx(L+2*W-1)];
        Yu=[yy(2)-yy(1) ; (yy(3:L+2*W)-yy(1:L+2*W-2))/2 ; yy(L+2*W)-yy(L+2*W-1)];
        m = Yu(W+1:L+W)./Xu(W+1:L+W);
        
        
    %smoothed_curve{i} = [sizex/2-yy(W+1:L+W)   sizey/2+xx(W+1:L+W)];
    smoothed_curve{i} = [xx(W+1:L+W)   yy(W+1:L+W)];
    gard_curve{i} =  atan(m);
    c_intersect{i} = yy(W+1:L+W) - gard_curve{i}.*xx(W+1:L+W);
    else
        smoothed_curve{i} = [];
        %smoothed{i} = [];
        gard_curve{i} = [];
        c_intersect{i} = [];
    end
end