function [out] = OLS( T,beta,alpha,rho,shi,delta )

     
        
         mu=[0 0];                  %Mean equal to zero
         sigma =[1 delta; delta 1]; %Covariance matrix

    
         uv= mvnrnd(mu,sigma,T);  
         u=uv(:,1); %random variable u with mean zero and varience 1
         v=uv(:,2); %random variable v with mean zero and varience 1

         x=zeros();
         r=zeros();
         x(1)=normrnd(0,1);

         for t= 2:T %simulate value for x & r


              x(t) = shi + rho .* x(t-1) + v(t);
              r(t)=alpha + beta .* x(t-1) + u(t);
 
         end

         R=r(2:length(r))';
         X=vertcat(ones(1,length(x)-1),x(1:length(x)-1))';

         OLS_Estimates_i=((inv(X'*X)) * X'*R); %Obtaining OLS  estimators
         variances =(diag(inv(X'*X)));         %Obtaining Variances  for estimators


         OLS_Estimate= OLS_Estimates_i(2);  % calulate beta hat
         variance = variances(2);           % calculate beta varience

          
         out=[ OLS_Estimate,variance];
    

end
