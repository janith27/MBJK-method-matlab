function [out] = PLUGINS(T,alpha,rho,shi,delta)

         beta=0;
         mu=[0 0];
         sigma =[1 delta; delta 1]; %Covariance matrix

         uv= mvnrnd(mu,sigma,T); 
         u=uv(:,1); %random variable u with mean zero and varience 1
         v=uv(:,2); %random variable v with mean zero and varience 1

         x=zeros();
         r=zeros();
         x(1)=normrnd(0,1);
         
        
         for t= 2:T %simulate value for x & r

              x(t) = shi +rho.* x(t-1) + v(t);
              r(t)= alpha + beta .* x(t-1) + u(t);
 
         end
        
         
         R=r(2:length(r))';
         X=vertcat(ones(1,length(x)-1),x(1:length(x)-1))';

         OLS_Estimates_i=((inv(X'*X)) * X'*R); %Obtaining OLS  estimators
         beta_hat_T= OLS_Estimates_i(2);    % calulate beta hat 
         
         
         
         Y=x(2:length(x))';


         OLS_Estimates_i=((inv(X'*X)) * X'*Y); %Obtaining OLS  estimators
         rho_hat_OLS= OLS_Estimates_i(2);
         Rho_OLS= rho_hat_OLS + (1+3*rho_hat_OLS)/T;
         

         
         variances =(diag(inv(X'*X)));%Obtaining Variances  for estimators
         variance = variances(2);     % calculate beta varience
         VAR_T=variance;   
        
       
       
    %end
    
     % calculate efcilent by evaluating residual

         mdl_u_t=fitlm(x(1:T-1),r(2:T));
         U_R(:,1)= mdl_u_t.Residuals.Raw;
        
         mdl_v_t=fitlm(x(1:T-1),x(2:T));
         V_R(:,1)= mdl_v_t.Residuals.Raw;
         var_v_t=var(V_R(:,1));
         
         new_sigma = cov(U_R(:,1),V_R(:,1));
         efcilent = new_sigma(1,2)/var_v_t;
       %  residual end
       
   
 
        plugin= -efcilent *(1+3* Rho_OLS)/T;
        out=[plugin,VAR_T,beta_hat_T];
          

end


  
