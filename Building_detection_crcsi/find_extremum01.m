function ext = find_extremum01(K);

% Find curvature local maxima as corner candidates
        extremum=[];
        N=size(K,1);
        %if (K(2) > K(1))
            Search=1;
        %else
        %    Search=-1;
        %end
        n=0;           
        
        
        for j=1:N-1
            if (K(j+1)-K(j))*Search>0
                n=n+1;
                extremum(n)=j;  % In extremum, odd points are minima and even points are maxima
                Search=-Search;
            end
        end
        
        ext = [];
        if (size(extremum,1)>0)
            if (extremum(1,1) > 1)
                ext = [ext 1];
            end
            ext = [ext extremum(2:2:n)];

            if (extremum(1,size(extremum,2)) < N && Search == -1)
                ext = [ext N];
            end
        end
%        if (size(extremum,2) == 0)
%            extremum(1,1) = 1;
%            extremum(1,2) = N-1;
%            extremum(1,3) = N;
%        end
        
%        if (size(extremum,2) == 1) % uphill only
%            extremum(1,2) = N-1;
%            extremum(1,3) = N;
%        end
       %{  
         if mod(size(extremum,2),2)==0            
            extremum = [extremum N];
         end
         
        
        if (extremum(1,size(extremum,2)) < N)
            extremum = [extremum N-1 N];
         end
         
        if (extremum(1,1) > 1)
            extremum = [1 2 extremum];
        end

        for i=1:floor(size(extremum,2)/2)
            %st = 
        end
         %}