ns = [10 20 40 80 100]

for i = 1:(numel(ns) - 1)
    x = linspace(0.2, 2, ns(i))';
    y1 = yexact(ns(i));
    y2 = beamCompact(ns(i));
    subplot(3, 4, [i+(i - 1), i+(i - 1) + 1])
    plot(x, y1,'r-o', x, y2, 'b-.*')
    legend('exact', 'approx')
    title(sprintf('Plot using n=%d steps', ns(i)))
    xlabel('length (m)')
    ylabel('height (m)')
end

i = 5;
x = linspace(0.2, 2, ns(i))';
y1 = yexact(ns(i));
y2 = beamCompact(ns(i));
subplot(3, 4, [10 11])
plot(x, y1,'r-o', x, y2, 'b-.*')
legend('exact', 'approx')
title(sprintf('Plot using n=%d steps', ns(i)))
xlabel('length (m)')
ylabel('height (m)')
