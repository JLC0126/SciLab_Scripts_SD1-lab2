function [] = retrato_de_fase(A)

    t = linspace(-2, 2, 101);  // variación del tiempo
    C1 = [-1 -1 -1 0 0 0 1 1 1];  // valores de las constantes
    C2 = [-1 0 1 -1 0 1 -1 0 1];
    
    [V, D] = spec(A);  // cálculo de los valores y vectores propios
    lambda1 = D(1, 1);
    lambda2 = D(2, 2);
    v1 = V(:, 1);
    v2 = V(:, 2);

    T = A(1, 1) + A(2, 2);
    D = A(1, 1) * A(2, 2) - A(2, 1) * A(1, 2);

    //////////////////////////////////////////////////////////
    if lambda1 ~= lambda2
        // caso 1: valores propios reales distintos, t^2-4d=0
        X1 = [exp(lambda1 * t) .* v1(1); exp(lambda2 * t) .* v1(2)];
        X2 = [exp(lambda2 * t) .* v2(1); exp(lambda2 * t) .* v2(2)];
    else
        if A(1, 2) == 0 && A(2, 1) == 0
            // caso valores propios reales repetidos, t^2-4d=0
            // si se obtienen dos vectores propios linealmente independientes,
            // en este caso A(1,2)==0 y A(2,1)==0:
            X1 = [exp(lambda1 * t) .* v1(1); exp(lambda1 * t) .* v1(2)];
            X2 = [exp(lambda1 * t) .* v2(1); exp(lambda1 * t) .* v2(2)];
        else
//si se obtiene un solo vector propio linealmente independiente,
//tenemos que encontrar otro vector P si A(1,2)==0 es distinto de cero
//P=[1;(v1(1)-(A(1,1)-lambda1))/A(1,2)], en caso contrario podemos
// definir P=[(v1(2)-(A(2,2)-lambda1))/A(2,1);1]
//solución:
            X1 = [exp(lambda1 * t) .* v1(1); exp(lambda1 * t) .* v1(2)];
            P = [1; (v1(1) - (A(1, 1) - lambda1)) / A(1, 2)];
            X2 = [exp(lambda1 * t) .* v1(1) + exp(lambda1 * t) .* P(1); exp(lambda1 * t) .* v1(2) + exp(lambda1 * t) .* P(2)];
        end
    end
    // caso 2: valores propios complejos
    alpha = real(lambda1);
    beta = imag(lambda1);
    x1 = real(v1(1));
    x2 = imag(v1(1));
    y1 = real(v1(2));
    y2 = imag(v1(2));
    X1 = [x1 * exp(alpha * t) .* cos(beta * t) - x2 * exp(alpha * t) .* sin(beta * t); 
          y1 * exp(alpha * t) .* cos(beta * t) - y2 * exp(alpha * t) .* sin(beta * t)];
    X2 = [x2 * exp(alpha * t) .* cos(beta * t) + x1 * exp(alpha * t) .* sin(beta * t); 
          y2 * exp(alpha * t) .* cos(beta * t) + y1 * exp(alpha * t) .* sin(beta * t)];
    // representación gráfica de soluciones
    ind = 1;
    while ind < length(C1)
        X = C1(ind) .* X1 + C2(ind) .* X2;
        x = X(1, :);
        y = X(2, :);
        plot(x, y, 'r->');
        xgrid();
        ind = ind + 1;
    end
    //campo de direcciones del sistema de ecuaciones diferenciales
    function [xdot] = SL(t, xy)
        xd1 = A(1,1)*xy(1) + A(1,2)*xy(2);
        xd2 = A(2,1)*xy(1) + A(2,2)*xy(2);
        xdot = [xd1; xd2];
    endfunction
    
    xf = -8:0.5:8;
    yf = -8:0.5:8;
    fchamp(SL, 0, xf, yf);
    xgrid();
    grafico = get('hdl');
    grafico.colored = 'on';
    ejes.isoview = 'on';

    xlabel('x(t)');
    ylabel('y(x)');
    title('Campo de direcciones y retrato de fase');

endfunction
