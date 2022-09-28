function [OUT] = MBJK( T,rho,delta)
     
   
    mu=[0 0]; 
    sigma =[1 delta; delta 1]; %Covariance matrix

    
    uv= mvnrnd(mu,sigma,T);  
    u=uv(:,1); %random variable u with mean zero and varience 1
    v=uv(:,2); %random variable v with mean zero and varience 1

    x=zeros();
    r=zeros();
    beta = 0;
    alpha =0;
    shi =0;
    
    L = 0.3*T;
    k = T - L + 1;
    
    x(1)=normrnd(0,1); 
    
    for t= 2:T %simulate value for x & r

        x(t) = shi + rho .* x(t-1) + v(t);
        r(t)=alpha + beta .* x(t-1) + u(t);
 
    end
    
    R=r(2:length(r))';
    X=vertcat(ones(1,length(x)-1),x(1:length(x)-1))';
    
    
    OLS_estimates=((inv(X'*X)) * X'*R); %calculate beta hat for each full sample
        
    OLS_beta= OLS_estimates(2);         % calulate beta hat
 
    
    %X & R split to sub samples 
    R_N=zeros();
    X_N=zeros();
    
    R_N = R;
    X_N = X(:,2);
    
   
    f=zeros();
    g=zeros();
    
    arr1=zeros();
     
    for j=1:k % allocate indexeses
        for n=1:L+1
            arr1(j,n)=j+n-1;
        end
    end
    
    for j=1:k  % assign valuses for the above indexeses 
        for n=1:L+1
           if arr1(j,n)<T
                f(j,n)=R_N(arr1(j,n),1);
                g(j,n)= X_N(arr1(j,n),1);
           else
               f(j,n)=NaN;
               g(j,n)=NaN;
           end
        end
       
    end
    
    f_new=zeros();
    g_new =zeros();
    
    for j=1:k % remove NaN valueses 
        for n=1:L+1
           if isnan(f(j,n))==false
                f_new(j,n)=f(j,n);
                g_new(j,n)=g(j,n);
                
           end
        end
       
    end 
  % subsample end  
  
    X_iwant=g_new';
    R_iwant=f_new';
    
    beta_hat_new=zeros();
  
   
    for col=1:k  %calculate beta hat for each sub samples
        
      X_new_column = X_iwant(:,col);
      X_trial = horzcat(ones(length(X_new_column),1),X_new_column);
      R_new = R_iwant(:,col);

      OLS_Estimates_new =((inv(X_trial'*X_trial)) * X_trial'*R_new);
      beta_hat_subsample = OLS_Estimates_new(2);

      beta_hat_new_col=(T*OLS_beta - L*beta_hat_subsample)/(T-L);
      beta_hat_new(col)=beta_hat_new_col;
       MBJK_sum= sum(beta_hat_new);
                    

     end
            
         MBJK = MBJK_sum/k;

        
        variances =(diag(inv(X'*X))); %Obtaining Variances  for estimators
        Var_T=variances(2);        %calculate beta varience
        

        OUT = [MBJK ,Var_T];
    
     
 
end
