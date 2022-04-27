x = 1;
hs = 2.^(-(2:12))';

t = table(hs);

u = -9.81*480*0.3*0.03/(1.3e10*(0.3*(0.03)^3/12));
t.u = repelem(u, numel(hs), 1);

t

t.uh = zeros(size(t.u));

% loop and populate values
for i = 1:numel(t.uh)
    h = t.hs(i);
    n = 2/h;
    idx = x/2*n;
    ys = beamCompact(n);
    t.uh(i) = vpa(ys(idx - 2)) - 4*ys(idx - 1) + 6*ys(idx) - 4*ys(idx + 1) + ys(idx+2);
    t.uh(i) = t.uh(i) / h^4;
end

t

for i = 1:numel(t.uh) - 1
    k1 = (t.uh(i) - t.u(i));
    t.uhh(i) = k1;
end

format long
t

% Fungsi Quadratic
fx = @(x) x^4 + x^3 - 5*x^2 + 3*x;

fxx = @(x) 4*x^3 + 3*x^2 - 10*x + 3;


hs = 2.^(-(0:20))';
v = table(hs);
v.u = zeros(size(v.hs));
v.uh1 = zeros(size(v.hs));
v.uh2 = zeros(size(v.hs));
v.uh3 = zeros(size(v.hs));

x = 4; % f(x) = 3*16+4+1 = 53 , f'(x) = 26
for i = 1:numel(v.u)
    h = v.hs(i);
    f_approx1 = @(x) (fx(x) - fx(x - h))./h; % backward difference
    f_approx2 = @(x) (fx(x+h) - fx(x))./h; % forward difference
    f_approx3 = @(x) (fx(x+h) - fx(x - h))./(2.*h); % central difference
    
    v.u(i) = fxx(x);
    v.uh1(i) = f_approx1(x);
    v.uh2(i) = f_approx2(x);
    v.uh3(i) = f_approx3(x);
end

v

% Calculate Order

for i = 1:(numel(v.u) - 1)
    v.uhh1(i) = (v.uh1(i) - v.u(i)) / (v.uh1(i + 1) - v.u(i));
    v.uhh2(i) = (v.uh2(i) - v.u(i)) / (v.uh2(i + 1) - v.u(i));
    v.uhh3(i) = (v.uh3(i) - v.u(i)) / (v.uh3(i + 1) - v.u(i));
end

v.order1 = log2(v.uhh1);
v.order2 = log2(v.uhh2);
v.order3 = log2(v.uhh3);
v

p1 = loglog(hs, abs(v.u - v.uh1), 'r-', hs, hs, 'b--*');

legend(p1, '|u - u_h|', 'h', 'location', 'best');
xlabel('h');
title('Error for Backward Difference')

p2 = loglog(hs, abs(v.u - v.uh2), 'r-', hs, hs, 'b--*');

legend(p2, '|u - u_h|', 'h', 'location', 'best');
xlabel('h');
title('Error for Forward Difference')

p3 = loglog(hs, abs(v.u - v.uh3), 'r-', hs, hs, 'b--*', hs, hs.^2, 'g-.d');

legend(p3, '|u - u_h|', 'h', 'h^2', 'location', 'best');
xlabel('h');
title('Error for Central Difference');
xlim([10^-4 10^0])
