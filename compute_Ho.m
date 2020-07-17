function Ho = compute_Ho(C,epsilon,sigma,NP) 
A = linspace(C(1), C(2), NP)';
B = linspace(epsilon(1), epsilon(2), NP)';
C = linspace(sigma(1), sigma(2), NP)';
Ho = [A(randperm(NP)) B(randperm(NP)) C(randperm(NP))];
end