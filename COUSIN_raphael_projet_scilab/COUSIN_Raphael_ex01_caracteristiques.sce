//ex1Q1
    t = linspace(-1, 1, 100);
    x0 = linspace (-2, 2, 5);

    clf();
    subplot(1, 2, 1);
    xtitle( 'cas a(x) = x','t','x')
    for i = 1 : 5
        plot2d(t, x0(i) * exp(t) ,[i]);
    end;
    subplot(1,2,2);
    xtitle( 'cas a(x) = -x','t','x')
    for i = 1 : 5
        plot2d(t, x0(i) * exp(-t), [i]);
    end;
    legends(['x0 = -2';'x0 = -1';'x0 = 0';'x0 = 1';'x0 = 2'],[1,2,3,4,5],opt="below");
    //xs2pdf(gcf(),"Q1");

