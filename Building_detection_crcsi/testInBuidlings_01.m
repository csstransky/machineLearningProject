function ret = testInBuidlings_01(Buildings,Seg);
ret = 1;
if (size(Buildings,1)>0)
    for i=1:size(Buildings,1)
        ret = testInRectangles(Buildings(i,1:2), Buildings(i,3:4), Buildings(i,5:6), Buildings(i,7:8), [Seg(1,1:2);Seg(1,3:4);Seg(1,8:9);Seg(1,10:11)]);        
        if (ret < 1)
            break;
        end
    end
end