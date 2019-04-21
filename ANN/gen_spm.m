function Output = gen_spm(center, radius)
% find the super patches fom the given radius
Output = zeros(length(center));

for i = 1:length(center(:, 1))
    temp1 = center(i, :);
    
    for j = 1:length(center(:, 1))
        temp2 = center(j, :);
        if (temp2(1)-temp1(1))^2+(temp2(2)-temp1(2))^2<=radius^2
            Output(i, j) = 1;
        end
    end
end

end

