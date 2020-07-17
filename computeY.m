function Y = computeY(n,ymin,ymax)
Q=length(ymin);
L=lhsdesign(n,Q);
Y=zeros(n,Q);
for i = 1:n
    for j = 1:Q
        Y(i,j)=(ymax(j)-ymin(j))*L(i,j) + ymin(j);
    end
end
end