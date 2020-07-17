function Ht = compute_Ht(NP,m,Ho,Fm,CR)

Hm = rand1(Ho,NP,Fm,m);

Ho_max = max(Ho);
Ho_min = min(Ho);

for n = 1:3
    if Hm(n) > Ho_max(n)
        Hm(n) = Ho_max(n);
    else
        if Hm(n) < Ho_min(n)
            Hm(n) = Ho_min(n);
        end
    end
end
   
jrand = floor(rand()*3+1);
for n = 1:3
    j = rand();
    if (j<=CR || n == jrand)
        Ht(n) = Hm(n);
    else
        Ht(n) = Ho(m,n);
    end
end
end
