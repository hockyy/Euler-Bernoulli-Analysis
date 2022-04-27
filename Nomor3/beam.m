function [A, x] = beam(n)
  A = getBeam(n);
  f = papanLoncat(n);
  b = sparse(n, 1);
  for i=1:n
    b(i) = f;
  end
  x = A\b;
endfunction

function [A] = getBeam(n)
  A = sparse(n, n);
  A(1, 1) =   16;
  A(1, 2) =   -9;
  A(1, 3) =  8/3;
  A(1, 4) = -1/4;
  A(2, 1) =   -4;
  A(2, 2) =    6;
  A(2, 3) =   -4;
  A(2, 4) =    1;
  for i=3:n-2
    A(i, i - 2) =  1;
    A(i, i - 1) = -4;
    A(i, i)     =  6;
    A(i, i + 1) = -4;
    A(i, i + 2) =  1;
  end
  A(n-1,n-3) =   16/17;
  A(n-1,n-2) =  -60/17;
  A(n-1,n-1) =   72/17;
  A(n-1,n)   =  -28/17;
  A(n,n-3)   =  -12/17;
  A(n,n-2)   =   96/17;
  A(n,n-1)   = -156/17;
  A(n,n)     =   72/17;
endfunction

function [f] = papanLoncat(n)

  L =     2; # 2m
  h = L / n;
  w =   0.3; # 30cm = 0.3m
  d =  0.03; # 3cm = 0.03m

  kerapatan =  480; # 480kg/m^3
  E =       1.3e10; # N/m^2
  I = w * d^3 / 12;
  g =         9.81; # Definisikan sebagai konstanta gravitasi

  f = (-480 * w * d * g);
  f =    h^4/(E * I) * f;

endfunction
