function Output = patMat_init(X)

Output = zeros(size(X));
Y = randperm(max(max(X)));

for i = 1:length(Y)
    Output(find(X == i)) = Y(i);
end

end

