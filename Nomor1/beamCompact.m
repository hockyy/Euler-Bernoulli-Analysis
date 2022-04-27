
function [x] = beamCompact(n)
  A = getEulerBernoulliBeamCompact(n);
  b = papanLoncat(n) * ones(n, 1);
  x = gaussianBeamCompact(A, b);
endfunction

function [A] = getEulerBernoulliBeamCompact(n)
  A = zeros(n, 7);
  A(1, 3:6) = [16, -9, 8/3, -1/4];
  A(2, 2:5) = [-4, 6, -4, 1];
  for i=3:n-2
    A(i, 1:5) = [1, -4, 6, -4, 1];
  end
  A(n-1, 1:4) = [16/17, -60/17, 72/17, -28/17];
  A(n, 1:3) = [3, -6, 3];
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

function [b] = gaussianBeamCompact(A, b)
  [n m] = size(A);
  b(n) += 0.75 * b(n - 1);
  for j=1:n
    pos = [j, 3];
    if(j + 1 <= n && abs(A(pos(1), pos(2))) < abs(A(j+1, 2)))
      pos = [j+1, 2];
    endif
    if(j + 2 <= n && abs(A(pos(1), pos(2))) < abs(A(j+2, 1)))
      pos = [j+2, 1];
    endif
    
    
    if(pos(2) != 3)
      A(j,:) = [A(j, 4-pos(2):7) zeros(1, 3-pos(2))];
      A(pos(1),:) = [zeros(1, 3-pos(2)) A(pos(1), 1:4 + pos(2))];
      
      tmp = A(j, :);
      A(j, :) = A(pos(1), :);
      A(pos(1), :) = tmp;
      
      tmp = b(j, :);
      b(j, :) = b(pos(1), :);
      b(pos(1), :) = tmp;
    endif
    
    
    for i=1:min(2, n-j)
      m = A(j + i, 3 - i)/A(j, 3);
      for k=1:(7-i)
        A(j + i, k) -= m * A(j, k + i);
      end
      b(j + i) -=  m * b(j);
    end
    
  end
  
  for i=n:-1:1
    for j=4:min(7, n - i + 3)
      b(i) -= b(i + j - 3) * A(i, j);
      A(i, j) = 0;
    end
    b(i) /= A(i, 3);
    A(i, 3) = 1;
  end
endfunction
