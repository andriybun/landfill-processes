function laplace_use
    clear, close all
    tmax = 100;
    [ft1, t1] = nilt(@flt, tmax);
    [t2, ft2] = invlap(@flt, 0.01, tmax);
    figure(3)
    set(3,'color','white')
    plot(t1, ft1), grid on, zoom on
    xlabel('t [s]'), ylabel('f(t)')
    hold on
    plot(t2, ft2, 'g'), grid on, zoom on
    
    hold off

    return
    
    function res = flt(s)
        eta = 500;
        res = 1 ./ s .* exp(-s.^(1/2) .* tanh(eta .* s.^(1/2)));
    end

end