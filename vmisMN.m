function a = vmisMN(M, N, min, max) % - горизонтальная
% MxN - размерность матрицы
%min, max - в каком диапазоне меняется ошибка
for i=1:M
    for j=1:N
        a(i,j)=0;
    end
end

step = ( max - min )/(M-1);

for i=1:M
    for j=1:N
        a(i,j)=a(i,j)+(1 * (min/100 + (i-1)*step/100));
        
    end
end
end