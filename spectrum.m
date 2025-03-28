    neigs = 200;
    Om = 1;
    kz = 2;
    op = 1;
    Om0 = op/(1-(op/kz)^2);
    op0minus = abs(Om)/2*(sqrt((kz/Om)^4 + 4*(kz/Om)^2)-(kz/Om)^2);
    op0plus =  abs(Om0)/2*(sqrt((kz/Om0)^4 + 4*(kz/Om0)^2)+(kz/Om0)^2);
    L = 20;
    N = 700;
    Nk = 150;
    c = floor(.1*N);
    periodic = true;
    kleft = -15;
    kright = 15;
    xk = linspace(kleft, kright, Nk);
    E1 = zeros([Nk, neigs]);
    E2 = zeros(size(E1));
    E3 = zeros(size(E1));
    E4 = zeros(size(E1));
    E5 = zeros(size(E1));
    E6 = zeros(size(E1));
    x = linspace(-L, L, 9*N);
    z = [0, 0, 0, 0, 0, 0, 0, 0, 1];
    B = [0, 0, 0, 0, 0, 0, 1, 1, 1];
    B = repmat(B, 1, N);
    Z = repmat(z, 1, N);
    parfor n = 1:Nk
        ky = xk(n);
        H = sparse(discH_B(ky, kz, Om, N, L));
        [v, e] = eigs(H, neigs, 0.8);
        %if periodic == false
            e = diag(e);
            for m = 1:neigs
                A = abs(v(:, m));
                mn = x*A/sum(A);
                std = sqrt(x.^2*A/sum(A) - mn^2);
                if norm(v(1:9*c, m))^2 + norm(v(9*(N-c):9*N, m))^2 > .5*norm(v(:, m))^2
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
        % else
        %     e = abs(diag(e));
        %     for m = 1:neigs
        %         if norm(v(1:9*c, m))^2 + norm(v(9*(N-c):9*N, m))^2 > .5*norm(v(:, m))^2
        %             E1((n-1)*neigs + m) = e(m);
        %         else
        %             E2((n-1)*neigs + m) = e(m);
        %         end
        %     end
        % end
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
    scatter(x1, E4, 15, 'g', ".");
    scatter(x1, E5, 15, 'cyan', ".");
    % scatter(x1, E6, 15, 'yellow', ".");
    legend({'spurious', 'edge mode', 'bulk left', 'bulk right', 'bulk both', 'TE'});
    hold off
    figure();
    hold on
    scatter(x1, E2, 15, 'r', ".");
    scatter(x1, E3, 15, 'r', ".");
    scatter(x1, E4, 15, 'r', ".");
    scatter(x1, E5, 15, 'r', ".");
    % scatter(x1, E6, 15, 'yellow', ".");
    legend({'edge mode', 'bulk left', 'bulk right', 'bulk both', 'TE'});
    hold off
% else
%     figure();
%     hold on
%     scatter(x1, E1, 15, 'k', ".");
%     scatter(x1, E2, 15, 'r', ".");
%     legend({'spurious', 'eigenvalue'});
%     hold off
%     figure();
%     scatter(x1, E2, 15, 'r', ".");
% % end