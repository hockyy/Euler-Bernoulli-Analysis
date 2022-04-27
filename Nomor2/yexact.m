function c = yexact(n)
%YEXACT Calculates correct score
%   Detailed explanation goes here
  L =     2; % 2m
  h = L / n;
  w =   0.3; % 30cm = 0.3m
  d =  0.03; % 3cm = 0.03m
  
  p =  480; % 480kg/m^3
  E =       1.3e10; % N/m^2
  I = w * d^3 / 12;
  g =         9.81; % Definisikan sebagai konstanta gravitasi
  
  f = -1*g*p*w*d;
  
  yx = @(x) f/24/E/I .* x.^2 .* (x.^2 - 4*L.*x + 6*L^2);
  ins = linspace(0.2, 2, n)';
  c = yx(ins)';
end
