function decision = check_vertical_alignment(pointA, pointB, windowSize, threshold, tolerance)

    horzDifference = abs(pointA(1)-pointB(1));
    horzCheck = horzDifference <= threshold;
    vertDistance = abs((pointA(2)+pointB(2))/2);
    if horzCheck
        thisTheta = abs((pointA(1)+pointB(1))/2);
        temp = (vertDistance*cosd(thisTheta) - windowSize/2) / tand(thisTheta);
        geoCheck = abs(temp - (windowSize/2 - vertDistance*sind(thisTheta))) / temp <= tolerance;
        if geoCheck && (vertDistance > windowSize/2) && (abs(pointA(2)-pointB(2)) > 7)
            decision = 1;
        else
            decision = 0;
        end
    else
        decision = 0;
    end
end