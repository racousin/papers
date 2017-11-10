//ex02 Q3

//initialisation matrice
function[A] = init_matA(N, h, lambda) //stockage plein
    for i = 1 : N - 1
        A(i, i) = 2 + 2 * lambda * h^2 / 3;
        A(i, i + 1) = -1 + lambda * h^2 / 6;
        A(i + 1, i) = -1 + lambda * h^2 / 6;
    end
    A(N,N) = 2 + 2 * lambda * h^2 / 3;
    A = 1 / h * A;
endfunction

function[b] = init_matb(N, h)
    for i = 1 : N - 1
        b(i) = 2 * sin(%pi * h * i) - sin(%pi * h * (i - 1)) - sin(%pi * h * (i + 1));
    end
    b(N) = 2 * sin(%pi * h * N) - sin(%pi * h * (N - 1));
    b = ((1 + %pi^2) / (h * %pi^2)) * b;
endfunction

//calcul solution approchee
function[u] = elem_fini(N)
    lambda = 1;
    h = 1 / (N + 1);
    x = [0 : h : 1];

    b = init_matb(N, h);
    A = init_matA(N, h, lambda);
    u = A \ b; //resolution du systeme
    u = [0; u ; 0];
endfunction

//erreurs
function[err2] = erreur_elem_L2(u)
    N = length(u) - 2;
    h = 1 / (N + 1);
    x = [0 : h : 1];
    uex = sin(%pi * x);

    err2 = 0;
    for k = 1 : N
        err2 = err2 + h * ((u(k + 1) - uex(k + 1))^2 + (u(k) - uex(k))^2)/2; //formule des trapezes
    end
endfunction

//affichage
function[] = print_elem(u)
    h = 1/(length(u) - 1);
    x = [0 : h : 1];
    clf;
    plot2d(x, u, 2); //solution approchee 
    xbis = [0 : 0.001 : 1]; // pour lisser la solution exacte
    plot2d(xbis, sin (%pi * xbis) , 1); // sol exacte
    legends(['solution apporchee'; 'solution exacte'],[2,1],opt="ur");
    xtitle('discretisation uniforme elements finis de pas h = ' + string(h),'x');
//    xs2pdf(gcf(),"ex02" + string(N));
endfunction

function[] = print_erreur

    err2 = linspace(-1,-1,100);
    for N = 1 : 100
        err2(N) = erreur_elem_L2(elem_fini(N));
    end
    plot2d(log([1:100]), log(err2), 1); //solution approchee 
    legends(['erreur L2'],[1],opt="ur");
    xtitle('courbe d''erreur en fonction de N, echelle log-log','N');
//    xs2pdf(gcf(),"ex02_erreur");
endfunction
print_erreur
//print_elem(elem_fini(1));
//print_elem(elem_fini(4));
//print_elem(elem_fini(7));
//print_elem(elem_fini(10));


