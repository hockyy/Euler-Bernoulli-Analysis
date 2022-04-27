ye = @(x) 9*x^2;

x = 4;
fp = @(x) 2.^(-x);
hs = fp(0:30)';
rs = ones(numel(hs), 1);
for i = 1:numel(hs)
    h = hs(i);
    ya_approx1 = @(x) (y(x + h)-y(x))/h;
    ya_approx2 = @(x) (y(x + h)-y(x -h))/(2*h);
    p1 = ya_approx1(x);
    p2 = ya_approx2(x);
    p3 = ye(x);
    rs(i) = log2(abs(p1 - p3));
end

format ShortEng
rs

f_uh = @(h) (y(4 + h)-y(4))./h;

t = table(hs);
t.uh = f_uh(t.hs);
t.uhh = zeros(size(t.uh));

for i = 1:numel(t.uh) - 1
    t.uhh(i) = t.uh(i) - t.uh(i + 1);
end

t.uhhh = zeros(size(t.uh));

for i = 1:numel(t.uh) - 2
    t.uhhh(i) = (t.uh(i) - t.uh(i + 1)) / (t.uh(i + 1) - t.uh(i + 2));
end

t.order = log2(t.uhhh)

x = 0.1
y = @(x) -9.81*480*0.3*0.03/(24*1.3e10*(0.3*(0.03)^3/12)).*x.^2.*(x.^2-8*x+24);
f_uh = @(h) (16*y(x) - 9*y(x+h) + 8/3*y(x + 2*h) - 1/4*y(x+3*h))./(h.^4);

%f_uh = @(h) (y(x - 2*h) - 4*y(x-h) + 6*y(x) - 4*y(x+h) + y(x+2*h))./(h.^4);

hs = fp(0:30)';
t = table(hs);

t.uh = f_uh(t.hs);
t.u = zeros(size(t.uh));
t.u = repelem(-9.81*480*0.3*0.03/(1.3e10*(0.3*(0.03)^3/12)), numel(t.uh),1);
t.uhh = zeros(size(t.uh));

for i = 1:numel(t.uh) - 1
    k1 = (t.uh(i) - t.u(i));
    k2 = (t.uh(i + 1) - t.u(i));
    t.uhh(i) = k1/k2;
end

format long
t.order = log2(t.uhh)
