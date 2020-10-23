function amatrix = gmisMN(M, N, a, b) % - horizontal
% MxN - размерность матрицы
% ax^2 + bx + c ; 
for i=1:M
    for j=1:N
        amatrix(i,j)=0;
    end
end

for i=1:M
    for j=1:N
        % NO SYMMETRY
  %         amatrix(i,j)=amatrix(i,j)+( a*(j-1)^2 + b*(j-1) )/100;
        % NICE SYMMETRY
            amatrix(i,j)=amatrix(i,j)+ ( a*(j - (N+1)/2 )^2 + b*(j-1)  )/100;
    end
end
end
