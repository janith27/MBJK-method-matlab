

clc
clear

alpha=0;
shi=0;
beta=0.5;


times=[60,120,360];
d=[-0.8,-0.9,-0.95];
rho_aray=[0.95,0.99,0.999];

% OLS

fprintf('_____________________________________________________________________________________OLS_________________________________________________________________________________________________\n');
fprintf('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
fprintf('%4s%3s%35s%23s%36s%24s%35s%25s\n','T','|','delta= -0.8','|','delta= -0.9','|','delta= -0.95','|');
fprintf('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
fprintf('%7s\t','|');

for i=1:3
    fprintf('%11s%6s%14s%6s%14s%6s\t','rho=0.95','|','rho=0.99','|','rho=0.999','|');
end
fprintf('\n');
fprintf('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
for p=1:length(times)
    
  T=times(p); 
  fprintf('\n');
  fprintf('%4d  |\t',T);
  
  for m=1:length(d)
     delta=d(m);
     
     for n=1:length(rho_aray)
        
         rho=rho_aray(n);
          for k=1:10000
                 [a(k,:)] = OLS( T,beta,alpha,rho,shi,delta ); %call OLS function
          end
  
        OLS_bias = mean([a(:,1)])-beta;
        OLS_RMSE = sqrt( OLS_bias^2 + mean([a(:,2)]));
        fprintf('%1.3f  (%.3f)  |\t',OLS_bias,OLS_RMSE);
   
    end
  end
end
fprintf('\n');
fprintf('\n');
fprintf('\n');
fprintf('\n');

% PLUGINS
fprintf('_____________________________________________________________________________________PLUGINS_____________________________________________________________________________________________\n');
fprintf('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
fprintf('%4s%3s%35s%23s%36s%24s%35s%25s\n','T','|','delta= -0.8','|','delta= -0.9','|','delta= -0.95','|');
fprintf('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
fprintf('%7s\t','|');

for i=1:3
    fprintf('%11s%6s%14s%6s%14s%6s\t','rho=0.95','|','rho=0.99','|','rho=0.999','|');
end
fprintf('\n');
fprintf('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
 
for p=1:length(times)
  T=times(p); 
  fprintf('\n');
  fprintf('%4d  |\t',T);
  for m=1:length(d)
     delta=d(m);
     
     for n=1:length(rho_aray)
         
         rho=rho_aray(n);
         for k=1:10000
            [b(k,:)] = PLUGINS(T,alpha,rho,shi,delta); 
         end
          
        Plugin_bias = mean([b(:,3)])-mean([b(:,1)]);
        Plugin_RMSE = sqrt( Plugin_bias^2 + mean([b(:,2)]));

        fprintf('%1.3f  (%.3f)  |\t', Plugin_bias,Plugin_RMSE);
         
   
    end
  end
end
fprintf('\n');
fprintf('\n');
fprintf('\n');
fprintf('\n');


% MBJK
fprintf('_____________________________________________________________________________________MBJK________________________________________________________________________________________________\n');
fprintf('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
fprintf('%4s%3s%35s%23s%36s%24s%35s%25s\n','T','|','delta= -0.8','|','delta= -0.9','|','delta= -0.95','|');
fprintf('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
fprintf('%7s\t','|');

for i=1:3
    fprintf('%11s%6s%14s%6s%14s%6s\t','rho=0.95','|','rho=0.99','|','rho=0.999','|');
end
fprintf('\n');
fprintf('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
for p=1:length(times)
  T=times(p); 
  fprintf('\n');
  fprintf('%4d  |\t',T);
  for m=1:length(d)
     delta=d(m);
     
     for n=1:length(rho_aray)
        
         rho=rho_aray(n);
         for k=1:10000
            [c(k,:)] =MBJK( T,rho,delta );
          end
        
        MBJK_bias=  mean([c(:,1)])  ;
        MBJK_RMSE = sqrt( MBJK_bias^2 + mean([c(:,2)]));
      
        fprintf('%1.3f  (%.3f)  |\t',MBJK_bias,MBJK_RMSE);
   
    end
  end
end
fprintf('\n');

 