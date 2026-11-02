function [logL] = loglik_2(PARAM,DATA)

yt = [DATA.yt]';
nt = [DATA.nt]';
lag_nt = [DATA.lag_nt]';

n = length(yt);
alpha1 = PARAM(1);
alpha2 = PARAM(2);
alpha3 = PARAM(3);
beta = PARAM(4);
sigma = PARAM(5);
    
    logL = 0;
for i=1:n
 logL = logL + log(normpdf(yt(i), (alpha1)*( exp(beta*nt(i)) - exp(beta*lag_nt(i)) ) ... 
     + (alpha2)*( exp(2*beta*nt(i)) - exp(2*beta*lag_nt(i)) ) ... 
     + (alpha3)*( exp(3*beta*nt(i)) - exp(3*beta*lag_nt(i)) ) ... 
     , sigma ));

end
logL= - logL;
end
