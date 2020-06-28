function a = gmisMN(M, N, min, max) % - горизонтальная
% MxN - размерность матрицы
%min, max - диапазон слева направо
for i=1:M
    for j=1:N
        a(i,j)=0;
    end
end

step = ( max - min )/(N-1);

for i=1:M
    for j=1:N
        a(i,j)=a(i,j)+(1 *(min/100 + (j-1)*step/100));
        
    end
end
end
