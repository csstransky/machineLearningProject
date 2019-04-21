function Thresh = find_Thresh(hThresh, Hhist);

T = sum(Hhist,2);
PThresh = 0.97; %percentage of pixels to be counted
PThresh_t = 0.9;
extremum = find_extremum(Hhist');
ThSize = floor(size(extremum,2)/2);

for i = 1:ThSize
    S{i} = Hhist(1,extremum(1,2*i-1):extremum(1,2*i+1));
    STh{i} = hThresh(1,extremum(1,2*i-1):extremum(1,2*i+1)+1);
    Ssum(i,1) = sum(S{i},2);
end

if (ThSize == 0)
        Thresh = hThresh;
else

        
    [stSum indices] = sort(Ssum,'descend');
    Tsum = 0;
    for i=1:ThSize
        id = indices(i,1);
        sm = stSum(i,1);
        smHist = S{id};
        smInd = STh{id};
        [stHist indx] = sort(smHist,'descend');
        tSum = 0;
        for j=1:size(stHist,2)
            tSum = tSum+stHist(1,j);
            per = tSum/sm;
            if (per>PThresh)
                break;
            end        
        end    
        indMat = indx(1,1:j);
        mn = min(indMat);
        mx = max(indMat);
        Thresh(i,1) = smInd(1,mn);
        Thresh(i,2) = smInd(1,mx+1);
        Tsum = Tsum+sm;
        Perc = Tsum/T;
        if (Perc>PThresh_t)
            break;
        end
    end

end
here = 1;