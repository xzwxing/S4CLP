function X = Norm2colum(Y)
% Unit the 2-norm of each column of Y
   Yvec= sqrt(sum(Y.^2));
   Yvec=repmat(Yvec,size(Y,1),1);
   X=Y./Yvec; 