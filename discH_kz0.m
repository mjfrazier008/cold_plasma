function H = discH_kz0(ky, op, N, L)
    dx = 2*L/N;
    c = floor(0.1*N);
    x = linspace(-L, L, N);

    % Defines how Om varies, change this line if desired:
    Omx = 1 - 2./(1+exp(-x));
    % Omx = heaviside(x)-0.5;
    % Omx = linspace(0.5, -0.5, N);

    % Ensures periodic boundary conditions:
    Om1 = Omx(c) + (Omx(N-c)-Omx(c))./(1+exp(3*linspace(-L, L, 2*c)));
    Omx(1:c) = Om1(c+1:2*c);
    Omx(N-c+1:N) = Om1(1:c);
    H = zeros(5*N);
    for n = 2:N-1
        Om = Omx(n);
        H(5*n-4:5*n, 5*(n-1):5*n+4) = [[0, 0, -1j*Om, 1j*op, 0, 0, 0, 0, 0, 0];
            [0, 1j*Om, 0, 0, 1j*op, 0, 0, 0, 0, 0];
            [-ky/2, -1j*op, 0, 0, 0, -ky/2,0, 0, 0, 0];
            [1j/dx, 0, -1j*op, 0, 0, -1j/dx, 0, 0, 0, 0];
            [0, 0, 0, -ky/2, 1j/dx, 0, 0, 0, -ky/2, -1j/dx]];
    end
    H(1:5, 1:9) = [[0, -1j*Omx(1), 1j*op, 0, 0, 0, 0, 0, 0];
            [1j*Omx(1), 0, 0, 1j*op, 0, 0, 0, 0, 0];
            [-1j*op, 0, 0, 0, -ky/2, 0, 0, 0, 0];
            [0, -1j*op, 0, 0, -1j/dx, 0, 0, 0, 0];
            [0, 0, -ky/2, 1j/dx, 0, 0, 0, -ky/2, -1j/dx]];
    H(3:4, 5*N) = [-ky/2; 1j/dx];
    H(5*N-4:5*N, 5*(N-1):5*N) = [[0, 0, -1j*Omx(N), 1j*op, 0, 0];
            [0, 1j*Omx(N), 0, 0, 1j*op, 0];
            [-ky/2, -1j*op, 0, 0, 0, -ky/2];
            [1j/dx, 0, -1j*op, 0, 0, -1j/dx];
            [0, 0, 0, -ky/2, 1j/dx, 0]];
    H(5*N, [3, 4]) = [-ky/2, -1j/dx];