%% ORL, 
  clear, clc, close all
  addpath(genpath('./'));
  
  load('ORL.mat'); % Load your data
  nX=size(fea,2); 
  class=max(gnd);
  fea= Norm2colum(fea); 
  lambda=1; alpha=0.1; beta=100; % for ORL
  
  %% 20 replicates of each experiment 
   Kpercent=[0.05,0.1,0.15,0.2,0.3,0.4];   
        
   repeat=20;
   Acc=zeros(repeat,length(Kpercent)); nMI=Acc; 
for i=1:length(Kpercent) 
   kpercent=Kpercent(i); 
   parfor k=1:repeat
      %%% Generate label constraint information: matrix Omega,
      Sel = Selction_sample(gnd,kpercent); % Randomly selected labeled samples
      Omega =Contr_matrix(nX,Sel,class); 
      
      tt=cell2mat(Sel); tt=sort( tt(:) );
      I_S=zeros(nX,1);  I_S(tt)=1; I_S=diag(I_S);
  
      Y0=zeros(nX,class);% label supervise information
      for j=1:class,   Y0(Sel{j},j)=1;  end  
    
      [Z,E,F] = S4CLP(fea,Omega,I_S,Y0,lambda,alpha,beta,1);
    
      [~,Grps] = max(F,[],2); % Obtain predicted labels           
      Acc(k,i) = sum(gnd(:) == Grps(:))/length(gnd);
      Nmi(k,i) = compute_nmi(double(gnd(:)),Grps(:));
      disp(['ACC(' num2str(k) ',' num2str(i) '),Nmi(' num2str(k) ',' num2str(i) '):=' num2str(Acc(k,i)) ',' num2str(Nmi(k,i))]);
   end
end

