%скрипт делает матрицу с градиентом и считает массив коэф-в
% "-1 это пустая ячейка, остальное - номер резистора
%____Надо Заполнить Матрицу

Rmatrix = [16 4 6 8 10 12 16 ;
14 5 7 9 11 13 14 ;
15 4 6 8 10 12 -1 ; 
-1 2 3 2 1 0 -1 ;
-1 -1 -1 -1 -1 -1 -1 ;
-1 0 1 2 3 2 -1 ;
-1 12 10 8 6 4 15 ;
14 13 11 9 7 5 14 ;
16 12 10 8 6 4 16 ;] ;

%без компенсации
% Rmatrix = [0 1 3 5 7 9 11 13 15 ; 
% -1 2 4 6 8 10 12 14 16 ;
% -1 2 4 6 8 10 12 14 16 ; ] ;

% 0 - R0, .. , 2 2 - в сумме R2
siz = size(Rmatrix);
M = siz(1); %число строк // raws
N = siz(2); %число столбцов // columns
razr = 8; %число разрядов_____надо вводить // bits

% quadratic gradient
quad_max_g = 0 ;
quad_max_v = 0 ;
% UNCOMMENT OUR FIGHTER (+fix vmisMN, gmisMN)
    % No symmetry
%     ag = quad_max_g / ( N-1 )^2 ;
%     av = quad_max_v / ( M-1 )^2 ;
    % Nice symmetry
    ag = quad_max_g / ( (N+1)/2 - 1 )^2 ;
    av = quad_max_v / ( (M+1)/2 - 1 )^2 ;

% linear gradient w/symmetry
lin_max_g = 8 ;
lin_max_v = 8 ;
cg = (-1)* lin_max_g ; 
cv = (-1)* lin_max_v ;
bg = 2 * lin_max_g / (N-1) ;
bv = 2 * lin_max_v / (M-1) ;

% matrix of nominals w/gradients ++++++++
res_real = ones(M, N) + gmisMN(M, N, ag, bg) + vmisMN(M,N,av,bv) + (cg + cv)/100 
% visual of gradients
figure;
mesh(res_real);

R0 = 0; % separate R0 to simplify numeration ++++
for i = 1 : M
    for j = 1 : N
        if ( Rmatrix(i,j) == 0 ) 
            R0 = R0 + res_real(i,j) ;
        end
    end
end
% get array R1 - R[2razr] (no R0) ++++++
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

 % Rp - массив "Rдо" ++++++
Rp(1) = R0 + R(1);
for k = 2:razr
    Rp(k) = 0;
    a=Rp(k-1); b=R(2*(k-1)); c=R(2*k - 1);
    %Rp(k) = par( Rp(k-1) , R(2*(k-1)) ) + R(2*k - 1) ;
    Rp(k) = par( a , b ) + c ;
end


 % Raf - массив "Rпосле" +++++
for k = 1:(razr-1)
    Raf(k) = 0;
end
Raf(razr - 1) = R(2*razr) + R(2*razr - 1);
for k = (razr-2):-1:1
    Raf(k) = 0;
    a=Raf(k+1); b=R(2*(k+1)); c=R(2*k+1);
    Raf(k) = par(a , b ) + c ;
end

%получаем массив коэф-в ++++++++ (ok w/ perfect R's)
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

% получаем Vout
Vsup = 1.8 ; % V лог. "1" _____ надо вводить
[Vout, Vout_perf] = voltageR(koef,Vsup);
lsb = Vout_perf(1);

%получаем DNL
DNL(1) = Vout(1)/lsb - 1;
for k = 2:length(Vout)
    DNL(k) = (Vout(k)-Vout(k-1))/lsb - 1;
end
%DNL_perf = DNL; use if perfect needed

figure;

plot(DNL);
legend('DNL');
grid on;

% получаем INL
INL = abs(Vout - Vout_perf)/lsb ;
%INL_perf = INL; ___ установить 
%INL = INL - INL_perf;
figure;

plot(INL);
legend('INL');
grid on;

