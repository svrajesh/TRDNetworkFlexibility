function val = logdet(A)

A = (A+A')/2;


[L,p] = chol(A,'lower');


if p == 0

    val = 2*sum(log(diag(L)));

else

    % Eigenvalue fallback

    eigVals = eig(A);

    eigVals(eigVals<=0) = eps;

    val = sum(log(eigVals));

end


end
