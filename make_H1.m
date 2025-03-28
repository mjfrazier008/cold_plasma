function [H, e, v] = make_H1(kp, op, Om, kz)
    kx = kp*cos(0);   
    ky = kp*sin(0);
    f = [0, -1j*Om, 0; 1j*Om, 0, 0; 0, 0, 0];
    kcr = [0, -kz, ky; kz, 0, -kx; -ky, kx, 0];
    I = eye(3)*1j;
    H = [f, op*I, zeros(3); -op*I, zeros(3), -kcr; zeros(3), kcr, zeros(3)];
    [v, e] = eig(H);