function [error, cond] = nomor3(n)
  error = zeros(n, 1);
  cond = zeros(n, 1);
  for k=1:n
    val = 10 * 2^k;
    [A, x] = beamSparse(val);
    err = x - solution(val);
    error(k) = abs(err(val));
    cond(k) = condest(A, 2);
    
    disp(val);
    disp(error(k));
    disp(cond(k));
  endfor
endfunction
