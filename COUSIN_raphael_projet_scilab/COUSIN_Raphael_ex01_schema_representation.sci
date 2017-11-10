//exercice 1 Q7-Q11

    xmin = -2;
    xmax = 2;
    T = 0.85;
    
//choix dx et dt

    dx = 0.2;//pas d'espace
    M = 1 / dx;//nombre de pas d'espace
    dt = 0.05;//pas de temps
    N = 0.85 / dt;//nombre de pas de temps
    lambda = dt / dx;
    if(2 * lambda * 2 <= 1)//pour 2 * lambda * max(a(x))
        disp('condition CFL OK, lambda = ' + string(lambda));
    else
        disp('pas de condition CFL, lambda = ' + string(lambda));
    end

//initialisation des variables et fonctions utiles
    x = [xmin : dx : xmax];
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
    xbis=[xmin : 0.01 : xmax]; //pour lisser en cas de pas grand
    L = length(x);

function[] = print_schema6(u0, a)
    ap = (max(a, zeros(1, length(a))));//a^+
    am = (max(-a, zeros(1, length(a))));//a^-
    u = u0(xbis);//solution exacte
    v = u0(x);//solution approchee

    if(a(1) < 0)//pour le cadre du plot
        R = [-2, -0.1, 2, 1];
    else
        R =[-1, -0.1, 1, 3];
    end

    clf;
    plot2d(x', v', 2, rect = R);
    plot2d(xbis', u', 1, rect = R);
    legends(['solution apporchee'; 'solution analytique'],[2,1],opt="ur");

//    xtitle('t = ' + string(0));
//    xs2pdf(gcf(),"ex1_s6_a" + string(sign(a(1))) + string(u0([])) + '_1');
    for n = 1 : N
//        clf;
        v(2 : L - 1) = lambda * ap(1 : L -2) .* v(1 : L - 2) + (1 - lambda * (ap(2 : L - 1) + am(1 : L - 2))) .* v(2 : L - 1) + lambda * am(2 : L - 1) .* v(3 : L);
        u = u0(xbis * exp(sign(a(1)) * n * dt)) * exp(sign(a(1)) * n * dt);//solution analytique
        plot2d(x', v', 2, rect = R);
        plot2d(xbis', u', 1, rect = R);
//        xtitle('t = ' + string(n * dt));
//        legends(['solution apporchee'; 'solution analytique'],[2,1],opt="ur");
//        xs2pdf(gcf(),"ex1_s6_a" + string(sign(a(1))) + string(u0([])) + '_' + string(1 + n));
    end
    xtitle('cas a(x) = ' + string(sign(-a(1)))+ 'x , ' + string(u0([])) + ' et dx = ' + string(dx) + ' , dt = ' + string(dt), 'x');
//    xs2pdf(gcf(),"ex1_s6_a" + string(sign(a(1))) + string(u0([])) + "_0");
endfunction

function[] = print_schema7(u0, a)
    ap = (max(a, zeros(1, length(a))));
    am = (max(-a, zeros(1, length(a))));
    u = u0(xbis);//solution exacte
    v = u0(x);//solution approchee

    if(a(1) < 0)//pour le cadre du plot
        R = [-2, -0.1, 2, 1];
    else
        R =[-1, -0.1, 1, 1];
    end
    clf;
    plot2d(x', v', 2, rect = R);
    plot2d(xbis', u', 1, rect = R);
    legends(['solution apporchee'; 'solution analytique'],[2,1],opt="ur");

//    xtitle('t = ' + string(0));
//    xs2pdf(gcf(),"ex1_s7_a" + string(sign(a(1))) + string(u0([])) + '_1');
    for n = 1 : N
//        clf;
        v(2 : L - 1) = lambda * ap(1 : L -2) .* v(1 : L - 2) + (1 - lambda * (am(2 : L - 1) + ap(1 : L - 2))) .* v(2 : L - 1) + lambda * am(2 : L - 1) .* v(3 : L);
        u = u0(xbis * exp(sign(a(1)) * n * dt));
        plot2d(x', v', 2, rect = R);
        plot2d(xbis', u', 1, rect = R);
//        xtitle('t = ' + string(n * dt));
//        legends(['solution apporchee'; 'solution analytique'],[2,1],opt="ur");
//        xs2pdf(gcf(),"ex1_s7_a" + string(sign(a(1))) + string(u0([])) + '_' + string(1 + n));
    end
    xtitle('cas a(x) ='+string(sign(-a(1)))+'x , '+ string(u0([])) + ' et dx = ' + string(dx)+' , dt = ' + string(dt), 'x');
//    xs2pdf(gcf(),"ex1_s7_a" + string(sign(a(1))) + string(u0([])) + "_0");
endfunction

print_schema6(u_ini, x); //Q7-8
print_schema6(u_ini, -x);//Q7-8
print_schema7(u_ini, x);//Q9-10
print_schema7(u_ini, -x);//Q9-10

print_schema6(u_lin, x); 
print_schema6(u_lin, -x);
print_schema7(u_lin, x);
print_schema7(u_lin, -x);
