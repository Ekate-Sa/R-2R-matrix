function a = vmisMN(M, N, min, max) % - ��������������
% MxN - ����������� �������
%min, max - � ����� ��������� �������� ������
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