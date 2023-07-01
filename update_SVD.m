function [U_dash, S_dash, V_dash] = update_SVD(U, S, V, a)
    [m, n] = size(S);
    if m<=n
        V1 = V(:, 1:m);
        V2 = V(:, m+1:n);
        D = S(1:m, 1:m);
        
        z = a'*U;
        
        [u, s, Q] = svd([D; z]);
        W = u(:, 1:m);
        w = u(:, m+1);
        Omega = s(1:m, 1:m);
        V1_base = [V1, zeros(n, 1); zeros(1, m), 1];
        V1_dash = V1_base * W;
        w_dash = V1_base * w;
        V2_temp = [V2; zeros(1, n-m)];
        V2_dash = [w_dash, V2_temp];
         
        U_dash = U * Q;
        S_dash = [Omega; zeros(n-m+1, m)]';
        V_dash = [V1_dash, V2_dash];
    else
        U1 = U(:, 1:n);
        U2 = U(:, n+1:m);
        D = S(1:n, 1:n);
        
        z1 = U1' * a;
        z2 = U2' * a;
        
        
        z2_normalized = z2 / norm(z2);
        
        reflectionVector = z2_normalized - [1; zeros(m-n-1, 1)];
        
        [W, Omega_dash, Q] = svd([D, zeros(n, 1); z1', norm(z2)]);
       

        U_dash2 = U2-2*((U2*reflectionVector)*reflectionVector')/(reflectionVector' * reflectionVector);
        v = U_dash2(1:m, 1);
        U2_dash = U_dash2(1:m, 2:m-n);
        
        U1_dash = [U1, v] * Q;
        
        U_dash = [U1_dash, U2_dash];
        S_dash = [Omega_dash, zeros(n+1, m-n-1)]';
        V_dash = [V, zeros(n, 1); zeros(1, n), 1]* W;
    end
end