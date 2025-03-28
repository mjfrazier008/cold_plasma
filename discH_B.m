function H = discH_B(ky, kz, op, N, L)
    dx = 2*L/N;
    c = floor(0.1*N);
    x = linspace(-L, L, N);
    periodic = true;
    if op < kz
        Om0= op/(1-(op/kz)^2);
    else
        Om0 = op/((op/kz)^2-1);
    end
    Omx =  0.75*Om0- 1.5*Om0./(1+exp(-x));
    % Omx = 1.2*Om0 - 0.4*Om0./(1+exp(-2*x));
    % Omx = linspace(-1.2*Om0, 1.2*Om0, N);
    % Omx = (heaviside(x)-1/2)*2.2*Om0;
    if periodic == true
        % opp1 = Omx(c)*exp(-x(c))/(1+ exp(-x(c)));
        % opp2 = Omx(N-c)*exp(-x(N-c))/(1+ exp(-x(N-c)));
        % x1 = linspace(x(N-c), L+(L-x(N-c)), 2*c);
        % p = spline([x1(1), x1(2*c)], [opp2 [Omx(N-c), Omx(c)] opp1]);
        % Omx(1:c) = ppval(p, x1(c+1:2*c));
        % Omx(N-c+1:N) = ppval(p, x1(1:c));
        Om1 = Omx(c) + (Omx(N-c)-Omx(c))./(1+exp(2*linspace(-L, L, 2*c)));
        Omx(1:c) = Om1(c+1:2*c);
        Omx(N-c+1:N) = Om1(1:c);
    end
    kcr1 = [0, -kz/2, ky/2; kz/2, 0, -1j/dx; -ky/2, 1j/dx, 0];
    kcr2 = [0, kz/2, -ky/2; -kz/2, 0, -1j/dx; ky/2, 1j/dx, 0];
    kcr3 = [0, -kz/2, ky/2; kz/2, 0, 1j/dx; -ky/2, -1j/dx, 0];
    kcr4 = [0, kz/2, -ky/2; -kz/2, 0, 1j/dx; ky/2, -1j/dx, 0];
    I = eye(3)*1j;
    H = zeros(9*N);
    for n = 2:N-1
        Om = Omx(n);
        f = [0, -1j*Om, 0; 1j*Om, 0, 0; 0, 0, 0];
        H(9*n-8:9*n, 9*n-11:9*n+6) = [zeros(3), f, -op*I, zeros(3),...
            zeros(3), zeros(3);
            kcr4, op*I, zeros(3), kcr2, zeros(3), zeros(3);
            zeros(3), zeros(3), kcr1, zeros(3), zeros(3) kcr3];
    end
    if periodic == true
        f = [0, -1j*Omx(1), 0; 1j*Omx(1), 0, 0; 0, 0, 0];
        H(1:9, 1:9) = [f, -op*I, zeros(3); op*I, zeros(3), kcr2;
            zeros(3), kcr1, zeros(3)];
        H(4:6, 9*N-2:9*N) = kcr4;
        H(7:9, 13:15) = kcr3;
        f = [0, -1j*Omx(N), 0; 1j*Omx(N), 0, 0; 0, 0, 0];
        H(9*N-8:9*N, 9*N-8:9*N) = [f, -op*I, zeros(3); op*I, zeros(3), kcr2;
            zeros(3), kcr1, zeros(3)];
        H(9*N-5:9*N-3, 9*N-11:9*(N-1)) = kcr4;
        H(9*N-2:9*N, 4:6) = kcr3;
    else
        f = [0, -1j*Omx(1), 0; 1j*Omx(1), 0, 0; 0, 0, 0];
        kcr2 = [0, kz/2, -ky/2; -kz/2, 0, 0; ky/2, 0, 0];
        H(1:9, 1:9) = [f, -op*I, zeros(3); op*I, zeros(3), kcr2;
            zeros(3), kcr1, zeros(3)];
        H(7:9, 13:15) = kcr3;
        f = [0, -1j*Omx(n), 0; 1j*Omx(n), 0, 0; 0, 0, 0];
        kcr2 = [0, kz/2, -ky/2; -kz/2, 0, -1j/dx; ky/2, 1j/dx, 0];
        kcr1 = [0, -kz/2, ky/2; kz/2, 0, 0; -ky/2, 0, 0];
        H(9*N-8:9*N, 9*N-8:9*N) = [f, -op*I, zeros(3); op*I, zeros(3), kcr2;
            zeros(3), kcr1, zeros(3)];
        H(9*N-5:9*N-3, 9*N-11:9*(N-1)) = kcr4;
    end