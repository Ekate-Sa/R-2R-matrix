function amatrix = vmisMN(M, N, a, b) % - vertical
% MxN - размерность матрицы
% ax^2 + bx + c ; x is percent per step
for i=1:M
    for j=1:N
        amatrix(i,j)=0;
    end
end

for i=1:M
    for j=1:N
        % NO SYMMETRY
%            amatrix(i,j)=amatrix(i,j)+ ( a*(i-1)^2 + b*(i-1)  )/100;
        % NICE SYMMETRY
            amatrix(i,j)=amatrix(i,j)+ ( a*(i - (M+1)/2 )^2 + b*(i-1)  )/100;
    end
end
end