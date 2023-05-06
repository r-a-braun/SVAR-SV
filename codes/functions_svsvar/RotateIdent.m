function Q = RotateIdent(A,r)
K = size(A,1);
GA = A';
Qp = eye(K);

for i = r+1:K-1
    for j = K:-1:i+1
        G = eye(K);
        rho = sign(GA(i,i)).*sqrt(GA(i,i)^2+GA(j,i)^2);
        c = GA(i,i)/rho;
        s = GA(j,i)/rho;
        G([i,j],[i,j]) = [c,s;-s,c];
        GA = G*GA;
        Qp = G*Qp;
    end
end
Q = Qp';
        