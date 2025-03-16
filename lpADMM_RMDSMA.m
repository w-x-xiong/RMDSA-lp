function [D,E,fail] = lpADMM_RMDSMA(D_hat,rho,p,gamma,H,maxiter)

sz = size(D_hat);
L = sz(1);
%initialization
D = D_hat;
D(H+1:end,:) = 0;
E = D - D_hat;
Delta = zeros(sz);
fail = false;
%k: counter for iterations
k = 0;
while 1
    k = k + 1;
    if k > maxiter
        fprintf('It cannot reach the convergence criterion in %d iterations\n', maxiter)
        fail = true;
        break
    end
    %ADMM section
    %update D
    N = D_hat+E-(1/rho)*Delta;
    [U_t,S_t,V_t] = svds(N,H);
    U_new = U_t*(S_t^(1/2));
    V_new = (S_t^(1/2))*V_t';
    D_new = U_new*V_new;
    
    %update E
    E_new = E;
    IMME_new = (1/rho)*Delta + D_new - D_hat;
    if (p==1)
        for i = 1:L
            for j = 1:L
                E_new(i,j) = max(IMME_new(i,j)-(1/rho),0) - max(-IMME_new(i,j)-(1/rho),0);
            end
        end
    elseif (p==2)
        for i = 1:L
            for j = 1:L
                E_new(i,j) = IMME_new(i,j)/(1+2*(1/rho));
            end
        end
    elseif ((p<2) && (p>1))
        for i = 1:L
            for j = 1:L
                E_new(i,j) = prox_lp(1/rho, IMME_new(i,j), p);
            end
        end
    end

    
    
    %update Delta
    Delta = Delta+rho*(D_new-D_hat-E_new);
    
    D = D_new;
    E = E_new;
    
    if ((norm(D_new-E_new-D_hat, 'fro')/norm(D_hat,'fro'))<gamma)
        D = D_new;
        break
    end
    
end
fprintf('The ADMM process involves iterating %d times\n', k)
end

