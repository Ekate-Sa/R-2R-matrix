%получает 2 массива Vout и строит графики
function [Vout,Vout_perf] = voltageR(koef,Vsup)

razr = length(koef);
N = 2^razr - 1;

for k = 1:razr
    koef_perf(k) = 1/2^k;
end
koef_perf = koef_perf(end:-1:1);

for i=1:N 
    Vout(i) = 0;
    Vout_perf(i) = 0;
    a = i;
    for k = 1:razr
        Vout(i) = Vout(i) + rem(a,2)*Vsup*koef(k);
        Vout_perf(i) = Vout_perf(i) + rem(a,2)*Vsup*koef_perf(k);
        a = floor(a/2);
    end
end

figure;
title('Выходная характеристика');
plot(Vout_perf, 'r--');
hold on;
plot(Vout);
grid on;
hold off;
end

