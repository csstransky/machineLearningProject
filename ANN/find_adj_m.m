function Output = find_adj_m(X)

[row, col]=size(X);
Output=zeros(row, col);

for i=1:row
    for j=1:col
        if abs(i-j)==1
            Output(i,j)=1;
        end
    end
end

end