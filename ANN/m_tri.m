function Output = m_tri(idx, Data, image)
% Calculate the mean of the tristimulus values of Data..

red_value = (Data == idx) .* image(:, :, 1);
green_value = (Data == idx) .* image(:, :, 2);
blue_value = (Data == idx) .* image(:, :, 3);

RED = mean(mean(red_value(find(red_value ~= 0))));
GREEN = mean(mean(green_value(find(green_value ~= 0))));
BLUE = mean(mean(blue_value(find(blue_value ~= 0))));

Output = [RED, GREEN, BLUE];
end