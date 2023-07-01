function [U_dash,S_dash,V_dash] = downdate_SVD(U,S,V)
    [m,n] = size(S);

    if m<=n
            V1 = V(:,1:m);
            V2 = V(:,m+1:n);
            D = S(1:m,1:m);
            
            U11 = U(2:m,1:m-1);
            u1 = U(1,1:m-1);
            x = U(2:m,m);
            my = U(1,m);
            
            M = eye(m-1)-((u1'*u1)/(1+my));
            X=U11*M-x*u1;
            
            C = [M -u1']*D;
            
            [Q,Omega,W] = svd(C);
            D_dash = Omega(1:m-1,1:m-1);
            
            U_dash = X*Q;
            S_dash = [D_dash zeros(m-1,n-m+1)];
            V_dash = [V1*W V2];
    
    else
            U1 = U(1:m,1:n);
            U2 = U(1:m,n+1:m);
            D = S(1:n,1:n);
            
            U11 = U1(2:m,1:n);
            u1 = U1(1,1:n)';
            U12 = U2(2:m,1:m-n);
            u2 = U2(1,1:m-n)';
            
            u2_normalized = u2 / norm(u2);
            reflectionVector = u2_normalized - [1; zeros(m-n-1, 1)];

            %%Code added to solve the problem when rotationVector is 0
            if reflectionVector == 0
                rotation_result = U12;
            else
                rotation_result = U12 - 2*(U12*reflectionVector)*reflectionVector'/(reflectionVector'*reflectionVector);
            end

            x=rotation_result(:,1);
            X12= rotation_result(:,2:end);
            my = norm(u2);
            
            U2_dash = X12;
            
            X = U11*(eye(length(u1)) - (u1*u1')/(1+my))-x*u1';
            
            C = (eye(length(u1)) - 1/(1+my)*(u1*u1'))*D;
            [Q,Omega,W] = svd(C);
            
            U1_dash = X*Q;
            D_dash = Omega;
            
            U_dash = [U1_dash, U2_dash];
            S_dash = [D_dash;zeros(m-n-1,n)];
            V_dash = V*W; 
    end    
end