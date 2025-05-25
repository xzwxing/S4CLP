function [Z,E,F] = S4CLP(X,Omega,I_S,Y,lambda,alpha,beta,DEBUG)
% 
% This routine solves the following l1-norm optimization problem by using ADMM 
%------------------------------
% min |Z|_1+lambda*|E|_1 + 0.5*alpha* sum_sum
%                  \|F_i-F_j\|^2(|z_ij|+|z_ji|)+beta*sum_{i \in Omega}\|F_i-Y_i\|^2
% s.t., X = XZ+E, Z_ij=0, (i,j) in \Omega
%      
%  By introduce new variable: A=Z-P_Omega(Z)
%
%  min |Z|_1+lambda*|E|_1+ alpha* sum_ij \|F_i-F_j\|^2 |z_ij| + beta*sum_{i \in Omega}\|F_i-Y_i\|^2
% s.t., X = XA+E,A=Z-P_Omega(Z)
%        
% inputs:
%        X -- D*N data matrix, D is the data dimension, and N is the number
%             of data vectors.
%        Z--N*N;  E--D*N;   F--k*N, k the number of classes
%----------by Zhiwei Xing
%
  tol_recErr = 2e-6;   tol_coe = 2e-6;
  maxIter =1000; 

  [d, n] = size(X); 
  [row,col]=find(Omega==1);
  Ind=sub2ind([n,n],row,col);

  rho = 1.1; % rho \in [1,1.5];
  eta =1.1;
  eta_max =8e20;%eta_max =4e10;
  normfX = norm(X,'fro');

%% Initializing optimization variables
  Z = zeros(n,n);
  A = zeros(n,n);
  E = zeros(d,n);
  F = zeros(size(Y));
  Delta1=zeros(d,n);
  Delta2=zeros(n,n);
  
%% Start main loop

  Inv= X'*X +eye(n) ;
  iter =0; 
  while iter<maxIter
     iter = iter + 1;
     Zk=Z;Ak=A;Ek=E;Fk=F;

     %% update Z  %%% 
     % THET=£¨thet_ij)_n*n with thet_ij=\|F_i-F_j\|^2 
     NormF= sum(F.*F,2); % F in R^{n*c}   
     THET=bsxfun(@plus,NormF',NormF) - 2*(F*F'); 
     
     temp=A-Delta2/eta;
     J=wthresh(temp,'s',(alpha*THET+1)/eta);
     P_omega_J=J(Ind);
     Z=J-sparse(row,col,P_omega_J,n,n);    
            
     %% update A  %%% 
     temp=Delta2/eta +X'*(X-E-Delta1/eta) +Z;  
     A=Inv\temp;
     P_omega_A=A(Ind);
     A=A-sparse(row,col,P_omega_A,n,n);
      
     %% update E %%%
     temp = X- X*A -Delta1/eta;
     E=wthresh(temp,'s',lambda/eta);
     
         %% update F %%%
     WW= BuildAdjacency(Z,10); % graph construction     
     L=diag(sum(WW,2))-WW;
     F=(L+ (beta/alpha)*I_S+0.00001*eye(n))\( (beta/alpha)*I_S*Y);
      
   %% stopping criteria  %%%
    %if iter==1 || mod(iter,10)==0 || iter>=50
    relChg_Z = max(max( abs(Z - Zk) ) );
    relChg_AZ = max(max( abs(A - Ak) ) );
    relChg_AP = max(max( abs(F - Fk) ) );
    relChg_E = max(max( abs( E - Ek) ) );
        
    relChg= max([relChg_Z,relChg_E,relChg_AZ,relChg_AP] );
    recErr=norm( X*A+E-X,'fro')/normfX;
    convergenced = recErr <tol_recErr && relChg < tol_coe;
    %end
    
    if DEBUG
       if iter==1 || mod(iter,50)==0 || convergenced
          disp(['iter ' num2str(iter) ',eta=' num2str(eta) ...
            ',recErr=' num2str(recErr) ',relChg=' num2str(relChg)]);
       end
    end
    if convergenced
        break;
    end  
     eta = min(eta_max, eta*rho);
     Delta1 = Delta1 + eta*( X*A+E-X );   
     Delta2 = Delta2 + eta*( Z-A ); 
  end

   
   