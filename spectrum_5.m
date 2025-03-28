%% Identical to spectrum.m but only for 5x5 transverse magnetic modes (kz = 0)

neigs = 100;
L = 15;
eigcenter = L*0.55;
Om = 1;
op = 1;
N = 2000;
Nk = 150;
c = floor(.1*N);
kleft = -1.5;
kright = 1.5;

xk = linspace(kleft, kright, Nk);
E1 = zeros([Nk, neigs]);
E2 = zeros(size(E1));
E3 = zeros(size(E1));
E4 = zeros(size(E1));
E5 = zeros(size(E1));
E6 = zeros(size(E1));
x = linspace(-L, L, 5*N);
z = [0, 0, 0, 0, 0, 0, 0, 0, 1];
B = [0, 0, 0, 0, 0, 0, 1, 1, 1];
B = repmat(B, 1, N);
Z = repmat(z, 1, N);
parfor n = 1:Nk
    ky = xk(n);
    H = sparse(discH_kz0(ky, op, N, L));
    [v, e] = eigs(L*H, neigs, eigcenter);
        e = diag(e);
        for m = 1:neigs
            A = abs(v(:, m));
            mn = x*A/sum(A);
            std = sqrt(x.^2*A/sum(A) - mn^2);
            if norm(v(1:5*c, m))^2 + norm(v(5*(N-c):5*N, m))^2 > .7*norm(v(:, m))^2
            %if std < 0.25*L && (mn < -.75*L || mn > .75*L)
                E1(n, m) = e(m);
            % else
            %     E2(n, m) = e(m);
            elseif std < .25*L
                E2(n, m) = e(m);
            elseif std < .5*L && mn < 0
                E3(n, m) = e(m);
            elseif std < .5*L && mn > 0
                E4(n, m) = e(m);
            % elseif norm(Z.*v(:, m)') < 0.1*(norm(B.*v(:, m)'))
            %     E6(n, m) = e(m);
            else
                E5(n, m) = e(m);
            end
        end
end

%% 

x1 = zeros([1, Nk*neigs]);
for n = 1:Nk
    x1(1+(n-1)*neigs:n*neigs) = ones(1, neigs)*xk(n);
end
E1 = reshape(E1', 1, Nk*neigs);
E2 = reshape(E2', 1, Nk*neigs);
E3 = reshape(E3', 1, Nk*neigs);
E4 = reshape(E4', 1, Nk*neigs);
E5 = reshape(E5', 1, Nk*neigs);
% E6 = reshape(E6', 1, Nk*neigs);
%if periodic == false
    figure();
    hold on
    scatter(x1, E1, 15, 'k', ".");
    scatter(x1, E2, 15, 'r', ".");
    scatter(x1, E3, 15, 'b', ".");
    scatter(x1, E4, 15, 'p', ".");
    scatter(x1, E5, 15, 'c', ".");
    % scatter(x1, E6, 15, 'yellow', ".");
    legend({'spurious', 'edge mode', 'bulk left', 'bulk right', 'bulk both', 'TE'});
    hold off
    figure();
    hold on
    scatter(x1, E2, 15, 'r', ".");
    scatter(x1, E3, 15, 'r', ".");
    scatter(x1, E4, 15, 'r', ".");
    scatter(x1, E5, 15, 'r', ".");
    %scatter(x1, E6, 15, 'yellow', ".");
    legend({'edge mode', 'bulk left', 'bulk right', 'bulk both', 'TE'});
    hold off
% else
    % figure();
    % hold on
    % scatter(x1, E1, 15, 'k', ".");
    % scatter(x1, E2, 15, 'r', ".");
    % legend({'spurious', 'eigenvalue'});
    % hold off
    % figure();
    % scatter(x1, E2, 15, 'r', ".");
% end