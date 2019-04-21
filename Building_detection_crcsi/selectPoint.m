%%%%%%%%%%% select points from extracted curves
function [extremum] = selectPoint(curve, curve_mode, curve_num, sizex, sizey);

extremum = [];
%S = 1.5*s;
s = 3.0;
S = 4.0;
[gau w] = find_Gaussian(s);
[Gau W] = find_Gaussian(S);
extra = W-w;
gau1 = [zeros(1,extra) gau zeros(1,extra)];
DoG = Gau-gau1;
K_thresh = 0.05;
%t = 0.1;

%smoothed_curve = cell(curve_num);
%point_selected = cell;
%Width = w; Sig = s;
%for i = 1:curve_num    
    x = curve(:,2) - sizey/2;
    y = sizex/2 - curve(:,1);
    L = size(x,1);
    if (L>W)
    % Calculate curvature
        if curve_mode=='loop'
            x1=[x(L-W+1:L);x;x(1:W)];
            y1=[y(L-W+1:L);y;y(1:W)];
        else
            x1=[ones(W,1)*2*x(1)-x(W+1:-1:2);x;ones(W,1)*2*x(L)-x(L-1:-1:L-W)];
            y1=[ones(W,1)*2*y(1)-y(W+1:-1:2);y;ones(W,1)*2*y(L)-y(L-1:-1:L-W)];
        end
       
        xx=conv(x1,DoG);
        xx=xx(2*W+1:L+2*W);
        yy=conv(y1,DoG);
        yy=yy(2*W+1:L+2*W);
        K = xx.^2 + yy.^2;
       
        % Find curvature local maxima as corner candidates
        extremum=[];
        N=size(K,1);
        n=0;
        Search=1;
        
        for j=1:N-1
            if (K(j+1)-K(j))*Search>0
                n=n+1;
                extremum(n)=j;  % In extremum, odd points is minima and even points is maxima
                Search=-Search;
            end
        end
        if mod(size(extremum,2),2)==0
            n=n+1;
            extremum(n)=N;
        end
    
        %n=size(extremum,2);
        %flag=ones(size(extremum));
  
        % Compare with adaptive local threshold to remove round corners
        %for j=2:2:n
        %    if K(extremum(j))<t
        %        flag(j)=0;
        %    end
        %end
        
        %plot(curve(:,1),curve(:,2),'-');
        %hold on
        %for k=1:size(extremum,2)
        %    plot(curve(extremum(1,k),1),curve(extremum(1,k),2),'+r');
        %end

        extremum= [extremum(2:2:n)];
        corners = K(extremum)>K_thresh;
        extremum = extremum(find(corners == 1));
        extremum = [1 extremum N];
        %figure; plot(curve(:,1),curve(:,2),'-');
        %hold on
        %for k=1:size(extremum,2)
        %    plot(curve(extremum(1,k),1),curve(extremum(1,k),2),'+r');
        %end
        
        %flag=flag(2:2:n);
        %extremum=extremum(find(flag==1));        
        %extremum=extremum-W;
        %extremum=extremum(find(extremum>0 & extremum<=L));
        
        %xx=conv(x1,gau);
        %xx=xx(W+1:L+3*W);
        %yy=conv(y1,gau);
        %yy=yy(W+1:L+3*W);
        
    %smoothed_curve{i} = [xx yy];
    %point_selected{i} = extremum;
    %Width = [Width W];
    %Sig = [Sig s];
    %else
     %   smoothed_curve{i} = [];
      %  point_selected{i} = [];
        %Width = [Width 0];
        %Sig = [Sig 0];
    end
    %here = 1;
%end