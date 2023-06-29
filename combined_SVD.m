function [U_dash, S_dash, V_dash] =combined_SVD(U, S, V, a)
    [m,n] = size(S);
    if m<=n
        u1 = U(1,1:m-1);
        M = eye(m-1)-((u1'*u1)/(1+U(1,m)));
        [Q,Omega,W] = svd([M -u1']*S(1:m,1:m));
        U = (U(2:m,1:m-1)*M-U(2:m,m)*u1)*Q;
        S = [Omega(1:m-1,1:m-1) zeros(m-1,n-m+1)];
        V = [V(:,1:m)*W V(:,m+1:n)];
    else
        U1 = U(1:m,1:n);
        U2 = U(1:m,n+1:m);
        u1 = U1(1,1:n)';
        U12 = U2(2:m,1:m-n);
        u2 = U2(1,1:m-n)';
        reflectionVector = u2 / norm(u2) - [1; zeros(m-n-1, 1)];
        if reflectionVector == 0
            temp_matrix = U12;
        else
            temp_matrix = U12 - 2*(U12*reflectionVector)*reflectionVector'/(reflectionVector'*reflectionVector);
        end
        my = norm(u2);                        
        [Q,Omega,W] = svd((eye(length(u1)) - 1/(1+my)*(u1*u1'))*S(1:n,1:n));            
        U = [(U1(2:m,1:n)*(eye(length(u1)) - (u1*u1')/(1+my))-temp_matrix(:,1)*u1')*Q, temp_matrix(:,2:end)];
        S = [Omega;zeros(m-n-1,n)];
        V = V*W; 
    end    
    m = m-1;
    if m<=n
        [u, s, Q] = svd([S(1:m, 1:m); a'*U]);
        V1_base = [V(:, 1:m), zeros(n, 1); zeros(1, m), 1];
        V2_temp = [ V(:, m+1:n); zeros(1, n-m)];
        V2_dash = [V1_base * u(:, m+1), V2_temp]; 
        V_dash = [V1_base * u(:, 1:m), V2_dash];
        S_dash = [s(1:m, 1:m); zeros(n-m+1, m)]';
        U_dash = U * Q;
    else
        U1 = U(:, 1:n);
        U2 = U(:, n+1:m);
        z2 = U2' * a;
        reflectionVector = z2 / norm(z2) - [1; zeros(m-n-1, 1)];
        [W, Omega_dash, Q] = svd([S(1:n, 1:n), zeros(n, 1); a'*U1, norm(z2)]);
        V_dash = [V, zeros(n, 1); zeros(1, n), 1]* W;
        U_dash2 = U2-2*((U2*reflectionVector)*reflectionVector')/(reflectionVector' * reflectionVector);
        U1_dash = [U1, U_dash2(1:m, 1)] * Q;
        U_dash = [U1_dash, U_dash2(1:m, 2:m-n)];
        S_dash = [Omega_dash, zeros(n+1, m-n-1)]';
    end
end