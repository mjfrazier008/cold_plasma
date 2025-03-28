%% Script for calculating spectrum of interface Hamiltonians for cold plasma
%
% 1. Choose discrete Hamiltonian function (discH_..., line 53) based on which
% parameter you would like to vary. B = vary magnetic field/cyclotron
% freqeuncy, opplus = vary plasma frequency around higher critical value,
% TCLW = vary plasma frequency around lower critical value, both = vary
% both op and Om. To change how parameters vary alter discH functions.
%
% 2. Choose parameter values (see below). Script will calculate interface
% spectrum for a range of k values, e.g. if only 1 k value is desired set
% kleft = kright and Nk = 1.
% 
% parameters:
%
% neigs: number of eigenvalues you want to calculate. suggest only
% ~10^2. without parallel computing even 50 may take an hour or more on my machine.
%
% eigcenter: script will calculate neigs nearest eigenvalues to eigcenter.
% Many times the best choice is the critical values om0plus or om0minus but
% may need to be adjusted based on density of states on either side of band
% gap.
%
% Om, kz, op: system parameters, see paper for details. 
%
% kleft, kright: lower and upper bound of perpendicular wave number
%
% Nk: Number of linearly spaced points between kleft and kright to diagonalize
%
% L: 1/2 length of domain centered around 0. 
%
% N: Number of discretization points.
%
% c: defines proportion of domain considered the spurious edge. recommend
% keeping at 10%. 

    neigs = 200;
    eigcenter = 0.8;
    Om = 1;
    kz = 2;
    op = 1;
    L = 20;
    N = 700;
    Nk = 150;
    c = floor(.1*N);
    kleft = -1.5;
    kright = 1.5;

    % You may find these values useful (critical values of Om and op):
    % Om0 = op/(1-(op/kz)^2);
    % op0minus = abs(Om)/2*(sqrt((kz/Om)^4 + 4*(kz/Om)^2)-(kz/Om)^2);
    % op0plus =  abs(Om)/2*(sqrt((kz/Om)^4 + 4*(kz/Om0)^2)+(kz/Om)^2);

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
        [v, e] = eigs(H, neigs, eigcenter);
            e = diag(e);
            for m = 1:neigs
                A = abs(v(:, m));

                % mean and std weight of eigenvectors used to categorize by
                % bulk modes, edge modes, and spurious edge modes:
                mn = x*A/sum(A);
                std = sqrt(x.^2*A/sum(A) - mn^2);
                
                % Spurious modes:
                if norm(v(1:9*c, m))^2 + norm(v(9*(N-c):9*N, m))^2 > .5*norm(v(:, m))^2
                %if std < 0.25*L && (mn < -.75*L || mn > .75*L)
                    E1(n, m) = e(m);
                % else
                %     E2(n, m) = e(m);
                % Edge modes:
                elseif std < .25*L
                    E2(n, m) = e(m);
                % Bulk mode left:
                elseif std < .5*L && mn < 0
                    E3(n, m) = e(m);
                % Bulk mode right:
                elseif std < .5*L && mn > 0
                    E4(n, m) = e(m);
                % Transverse Electric mode:
                % elseif norm(Z.*v(:, m)') < 0.1*(norm(B.*v(:, m)'))
                %     E6(n, m) = e(m);
                % Bulk both left and right:
                else
                    E5(n, m) = e(m);
                end
            end
    end
%% Plot your results:
% E1: spurious, E2: edge, E3: bulk left, E4: bulk right, E5: bulk both, E6:
% transverse electric. 

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