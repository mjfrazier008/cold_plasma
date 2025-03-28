function H = discH_kplus1(ky, op, Om, N, L, periodic)
    dx = 2*L/N;
    c = floor(0.1*N);
    kplus = op/sqrt(1+ op/Om);
    x = linspace(-L, L, N);
    %opx = linspace(0, 2*op0, N);
    kzx = 0.9*kplus + 0.2*kplus./(1+ exp(-x));
    %opx = linspace(1.5, 3, N);
    if periodic == true
        % opp1 = opx(c)*exp(-x(c))/(1+ exp(-x(c)));
        % opp2 = opx(N-c)*exp(-x(N-c))/(1+ exp(-x(N-c)));
        % x1 = linspace(x(N-c), L+(L-x(N-c)), 2*c);
        % p = spline([x1(1), x1(2*c)], [opp2 [opx(N-c), opx(c)] opp1]);
        % opx(1:c) = ppval(p, x1(c+1:2*c));
        % opx(N-c+1:N) = ppval(p, x1(1:c));
        kzx1 = kzx(c) + (kzx(N-c)-kzx(c))./(1+exp(linspace(-L, L, 2*c)));
        kzx(1:c) = kzx1(c+1:2*c);
        kzx(N-c+1:N) = kzx1(1:c);
    end
    f = [0, -1j*Om, 0; 1j*Om, 0, 0; 0, 0, 0];
    I = eye(3)*1j;
    H = zeros(9*N);
    for n = 2:N-1
        kz = kzx(n);
        kcr1 = [0, -kz/2, ky/2; kz/2, 0, -1j/dx; -ky/2, 1j/dx, 0];
        kcr2 = [0, kz/2, -ky/2; -kz/2, 0, -1j/dx; ky/2, 1j/dx, 0];
        kcr3 = [0, -kz/2, ky/2; kz/2, 0, 1j/dx; -ky/2, -1j/dx, 0];
        kcr4 = [0, kz/2, -ky/2; -kz/2, 0, 1j/dx; ky/2, -1j/dx, 0];
        H(9*n-8:9*n, 9*n-11:9*n+6) = [zeros(3), f, -op*I, zeros(3),...
            zeros(3), zeros(3);
            kcr4, op*I, zeros(3), kcr2, zeros(3), zeros(3);
            zeros(3), zeros(3), kcr1, zeros(3), zeros(3) kcr3];
    end
    if periodic == true
        kz = kzx(1);
        kcr1 = [0, -kz/2, ky/2; kz/2, 0, -1j/dx; -ky/2, 1j/dx, 0];
        kcr2 = [0, kz/2, -ky/2; -kz/2, 0, -1j/dx; ky/2, 1j/dx, 0];
        kcr3 = [0, -kz/2, ky/2; kz/2, 0, 1j/dx; -ky/2, -1j/dx, 0];
        kcr4 = [0, kz/2, -ky/2; -kz/2, 0, 1j/dx; ky/2, -1j/dx, 0];
        H(1:9, 1:9) = [f, -op*I, zeros(3); op*I, zeros(3), kcr2;
            zeros(3), kcr1, zeros(3)];
        H(4:6, 9*N-2:9*N) = kcr4;
        H(7:9, 13:15) = kcr3;
        kz = kzx(N);
        kcr1 = [0, -kz/2, ky/2; kz/2, 0, -1j/dx; -ky/2, 1j/dx, 0];
        kcr2 = [0, kz/2, -ky/2; -kz/2, 0, -1j/dx; ky/2, 1j/dx, 0];
        kcr3 = [0, -kz/2, ky/2; kz/2, 0, 1j/dx; -ky/2, -1j/dx, 0];
        kcr4 = [0, kz/2, -ky/2; -kz/2, 0, 1j/dx; ky/2, -1j/dx, 0];
        H(9*N-8:9*N, 9*N-8:9*N) = [f, -op*I, zeros(3); op*I, zeros(3), kcr2;
            zeros(3), kcr1, zeros(3)];
        H(9*N-5:9*N-3, 9*N-11:9*(N-1)) = kcr4;
        H(9*N-2:9*N, 4:6) = kcr3;
    else
        kz = kzx(1);
        kcr2 = [0, kz/2, -ky/2; -kz/2, 0, 0; ky/2, 0, 0];
        kcr1 = [0, -kz/2, ky/2; kz/2, 0, -1j/dx; -ky/2, 1j/dx, 0];
        kcr3 = [0, -kz/2, ky/2; kz/2, 0, 1j/dx; -ky/2, -1j/dx, 0];
        kcr4 = [0, kz/2, -ky/2; -kz/2, 0, 1j/dx; ky/2, -1j/dx, 0];
        H(1:9, 1:9) = [f, -opx(1)*I, zeros(3); opx(1)*I, zeros(3), kcr2;
            zeros(3), kcr1, zeros(3)];
        H(7:9, 13:15) = kcr3;
        kz = kzx(N);
        kcr2 = [0, kz/2, -ky/2; -kz/2, 0, -1j/dx; ky/2, 1j/dx, 0];
        kcr3 = [0, -kz/2, ky/2; kz/2, 0, 1j/dx; -ky/2, -1j/dx, 0];
        kcr4 = [0, kz/2, -ky/2; -kz/2, 0, 1j/dx; ky/2, -1j/dx, 0];
        kcr1 = [0, -kz/2, ky/2; kz/2, 0, 0; -ky/2, 0, 0];
        H(9*N-8:9*N, 9*N-8:9*N) = [f, -opx(N)*I, zeros(3); opx(N)*I, zeros(3), kcr2;
            zeros(3), kcr1, zeros(3)];
        H(9*N-5:9*N-3, 9*N-11:9*(N-1)) = kcr4;
    end