function [sol] = solution(n)

L =     2; # 2m
h = L / n;

w =   0.3; # 30cm = 0.3m

d =  0.03; # 3cm = 0.03m
E = 1.3e10;
I = w * d^3 / 12;

g =         9.81; # Definisikan sebagai konstanta gravitasi

f = (-480 * w * d * g);

correct = zeros(n, 1);

for k=1:n
    x = h * k;
    sol(k,1) = f/(24*E*I) * x^2 * (x^2 - 8*x + 24);
end