//ex1 Q11 convergence au temps T
    xmin = -2;
    xmax = 2;
    T = 0.85;

function[y] = u_ini(x)
    if length(x) == 0
        y = 'u_ini';
    else
        y = bool2s((abs(x)<=1/2));
    end
endfunction

function[y] = u_lin(x)
    if length(x) == 0
        y = 'u_lin';
    else
        y = bool2s(abs(x)<= 1/2) .* ones(1,length(x)) + 2*x.*bool2s((x>= -1/2) & (x<=0)) - 2*x.*bool2s((x<= 1/2) & (x>=0));
    end
endfunction

function[erreur1, erreur8] = erreur_schema6(u0, a, dx, dt)
    N = 0.85 / dt;
    lambda = dt / dx;
    x = [xmin : dx : xmax];
    L = length(x);
    ap = (max(a, zeros(1, length(a))));//a^+
    am = (max(-a, zeros(1, length(a))));//a^-
    u = u0(x);//solution exacte
    v = u0(x);//solution approchee

    for n = 1 : N
        v(2 : L - 1) = lambda * ap(1 : L -2) .* v(1 : L - 2) + (1 - lambda * (ap(2 : L - 1) + am(1 : L - 2))) .* v(2 : L - 1) + lambda * am(2 : L - 1) .* v(3 : L);
        u = u0(x * exp(sign(a(1)) * n * dt)) * exp(sign(a(1)) * n * dt);//solution analytique
    end

    erreur1 = dx * sum (abs(u - v));//norme L1
    erreur8 = max (abs (u - v));// norme infini
endfunction

function[erreur1, erreur8] = erreur_schema7(u0, a, dx, dt)
    N = 0.85 / dt;
    lambda = dt / dx;
    x = [xmin : dx : xmax];
    L = length(x);
    ap = (max(a, zeros(1, length(a))));//a^+
    am = (max(-a, zeros(1, length(a))));//a^-
    u = u0(x);//solution exacte
    v = u0(x);//solution approchee

    for n = 1 : N
        v(2 : L - 1) = lambda * ap(1 : L -2) .* v(1 : L - 2) + (1 - lambda * (am(2 : L - 1) + ap(1 : L - 2))) .* v(2 : L - 1) + lambda * am(2 : L - 1) .* v(3 : L);
        u = u0(x * exp(sign(a(1)) * n * dt));
    end

    erreur1 = dx * sum (abs(u - v));//norme L1
    erreur8 = max (abs (u - v));// norme infini
endfunction

function[] = print_erreur(u0)
    dx = 0.2;
    dt = 0.05;
    lambda = dt / dx;
    if(2 * lambda * 2 <= 1)//pour 2 * lambda * max(a(x))
        disp('condition CFL OK, lambda = ' + string(lambda));
    else
        disp('pas de condition CFL, lambda = ' + string(lambda));
    end
    xabs = linspace(-1, -1, 10);
    erreur1 = -ones(4,10)
    erreur8 = -ones(4,10);
    //cas a(x) = x
    for i = 1 : 10
        x = [xmin : dx : xmax];
        [erreur1(1,i), erreur8(1,i)] = erreur_schema6(u0, x, dx, dt);
        [erreur1(2,i), erreur8(2,i)] = erreur_schema6(u0, -x, dx, dt);
        [erreur1(3,i), erreur8(3,i)] = erreur_schema7(u0, x, dx, dt);
        [erreur1(4,i), erreur8(4,i)] = erreur_schema7(u0, -x, dx, dt);
        xabs(i) = dx + dt;
        dx = dx / 2;
        dt = dt / 2;
   end
    clf();
    plot2d(xabs', [erreur1(1,:)' erreur8(1,:)' erreur1(2,:)' erreur8(2,:)'], [1, 2, 3, 4]);
    legends(['cas a(x) = x norme L1'; 'cas a(x) = x norme infini'; 'cas a(x) = -x norme L1'; 'cas a(x) = -x norme infini'],[1, 2, 3, 4],opt = "below");
    xtitle('Erreur solution analytique et schema6 au temps T = 0.85 pour ' + u0([]), 'dx + dt');
//    xs2pdf(gcf(),"ex1_conv_s6'+u0([]));
    figure("BackgroundColor",[1,1,1]);
    plot2d(xabs', [erreur1(3,:)' erreur8(3,:)' erreur1(4,:)' erreur8(4,:)'], [1, 2, 3, 4]);
    legends(['cas a(x) = x norme L1'; 'cas a(x) = x norme infini'; 'cas a(x) = -x norme L1'; 'cas a(x) = -x norme infini'],[1, 2, 3, 4],opt = "below");
    xtitle('Erreur solution analytique et schema7 au temps T = 0.85 pour ' + u0([]), 'dx + dt');
//    xs2pdf(gcf(),"ex1_conv_s7'+u0([]));
endfunction

print_erreur(u_ini);// Q11
print_erreur(u_lin);//Q11
