function [prox] = prox_lp(tau, b, p)
%proximal mapping of lp-norm (for 1<p<2)

if (b >= 0)

    fun1 = @(a) (a-b)/tau+p*a^(p-1);
    a0 = [0 b];
    aplus = fzero(fun1,a0);
    objaplus = abs(aplus)^p+(1/2/tau)*(aplus-b)^2;
    if (objaplus <= (b^2/2/tau))
        prox = aplus;
    else
        prox = 0;
    end

else

    fun2 = @(a) (a-b)/tau-p*(-a)^(p-1);
    a0 = [b 0];
    aminus = fzero(fun2,a0);
    objaminus = abs(aminus)^p+(1/2/tau)*(aminus-b)^2;
    if (objaminus <= (b^2/2/tau))
        prox = aminus;
    else
        prox = 0;
    end

end

end

