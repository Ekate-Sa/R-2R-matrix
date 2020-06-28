%������ ������ ������� � ���������� � ������� ������ ����-�
% "-1 ��� ������ ������, ��������� - ����� ���������
%____���� ��������� �������

Rmatrix = [16 4 6 8 10 12 16 ;
14 5 7 9 11 13 14 ;
15 4 6 8 10 12 -1 ; 
-1 2 3 2 1 0 -1 ;
-1 -1 -1 -1 -1 -1 -1 ;
-1 0 1 2 3 2 -1 ;
-1 12 10 8 6 4 15 ;
14 13 11 9 7 5 14 ;
16 12 10 8 6 4 16 ;] ;

% 0 - R0, .. , 2 2 - � ����� R2
siz = size(Rmatrix);
M = siz(1); %����� �����
N = siz(2); %����� ��������
razr = 8; %����� ��������_____���� �������

min_g = 4; % min - �����, max - ������; _____���� �������
max_g = 18;
min_v = 10; % ��� ������������ ������ ����
max_v = -6;
%�������� ������� � ���������� ++++++++
res_real = ones(M, N) + gmisMN(M, N, min_g, max_g) + vmisMN(M,N,min_v,max_v) 

R0 = 0; % �������� R0, ����� ���� ����� � ���������� ++++
for i = 1 : M
    for j = 1 : N
        if ( Rmatrix(i,j) == 0 ) 
            R0 = R0 + res_real(i,j) ;
        end
    end
end
% �������� ������ R1 - R16 ++++++
for k = 1:razr*2
    R(k) = 0;
    for i = 1 : M
        for j = 1 : N
            if ( Rmatrix(i,j) == k )
                R(k) = R(k) + res_real(i,j) ;
            end
        end
    end
end

 % Rp - ������ "R��" ++++++
Rp(1) = R0 + R(1);
for k = 2:razr
    Rp(k) = 0;
    a=Rp(k-1); b=R(2*(k-1)); c=R(2*k - 1);
    %Rp(k) = par( Rp(k-1) , R(2*(k-1)) ) + R(2*k - 1) ;
    Rp(k) = par( a , b ) + c ;
end


 % Raf - ������ "R�����" +++++
for k = 1:(razr-1)
    Raf(k) = 0;
end
Raf(razr - 1) = R(2*razr) + R(2*razr - 1);
for k = (razr-2):-1:1
    Raf(k) = 0;
    a=Raf(k+1); b=R(2*(k+1)); c=R(2*k+1);
    Raf(k) = par(a , b ) + c ;
end

%�������� ������ ����-� ++++++++ (��� ��������� R ���������� ���������)
for k = 1:razr
    koef(k) = 0;
end

koef(razr) = del(R(2*razr) , Rp(razr) ) ;
for i = 1:(razr-1)
    koef(i) = del( par(Rp(i),Raf(i)) , R(2*i))*del(R(2*razr),R(2*razr-1)) ;
    q = 1;
    for k = (i+1):(razr-1)
        a=Raf(k); b=R(2*k); c=R(2*k-1);
        q = q * del( par(a,b) , c );
    end
    koef(i) = koef(i) * q ;
end

% �������� Vout
Vsup = 1.8 ; % V ���. "1" _____ ���� �������
[Vout, Vout_perf] = voltageR(koef,Vsup);
lsb = Vout_perf(1);

%�������� DNL
DNL(1) = Vout(1)/lsb - 1;
for k = 2:length(Vout)
    DNL(k) = (Vout(k)-Vout(k-1))/lsb - 1;
end
%DNL_perf = DNL; ___ ���������� 
DNL = DNL-DNL_perf;

figure;

plot(DNL);
legend('DNL');
grid on;

% �������� INL
INL = abs(Vout - Vout_perf)/lsb ;
%INL_perf = INL; ___ ���������� 
INL = INL - INL_perf;
figure;

plot(INL);
legend('INL');
grid on;

