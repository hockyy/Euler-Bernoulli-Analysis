function [A, x] = beamSparse(n)
  A = getEulerBernoulliBeamSparse(n);
  f = papanLoncat(n);
  b = sparse(n, 1);
  for i=1:n
    b(i) = f;
  end
  x = gaussianBeamSparse(A, b);
  % norm(x - A \ b)
endfunction

function [A] = getEulerBernoulliBeamSparse(n)
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

  L =     2; % 2m
  h = L / n;
  w =   0.3; % 30cm = 0.3m
  d =  0.03; % 3cm = 0.03m

  kerapatan =  480; % 480kg/m^3
  E =       1.3e10; % N/m^2
  I = w * d^3 / 12;
  g =         9.81; % Definisikan sebagai konstanta gravitasi

  f = (-480 * w * d * g);
  f =    h^4/(E * I) * f;

endfunction

function [b] = gaussianBeamSparse(A, b)
  [n n] = size(A);
  A = [A b];
  A(n, :) += 0.75 * A(n - 1, :);
  for j=1:n
    pos = j;
    if(j + 1 <= n && abs(A(j+1,j)) > abs(A(pos, j))) pos = j+1;
    elseif(j + 2 <= n && abs(A(j+2, j)) > abs(A(pos, j))) pos = j+2;
    endif
    tmp = A(j, :);
    A(j, :) = A(pos, :);
    A(pos, :) = tmp;
    for i=j+1:min(n, j+2)
      m = A(i, j)/A(j, j);
      A(i, :) -= m * A(j, :);
    end
  end

  for i=n:-1:1
    for j=i+1:min(i+4, n)
      A(i, n+1) -= A(j, n + 1) * A(i, j);
      A(i, j) = 0;
    end
    A(i, n+1) /= A(i, i);
    A(i, 2) = 1;
  end
  b = A(:, n+1);
endfunction
