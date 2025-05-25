 
function Sel = Selction_sample(gnd,percent)
%%% Randomly select some samples as supervisory information
%   gnd£ºgroundtruth£¬
%   percent: ratio of selected samples

    class=max(gnd);
    Ind = cell(class,1);
    for j=1:class
      Ind{j} = find(gnd==j); 
    end
  
      Sel=cell(class,1);
   for j=1:class
       num_j=numel(Ind{j});
       KLabel=round( num_j *percent );
       temp= Ind{j}( sort( randperm(num_j,KLabel ) ) );
       Sel{j}=temp(:);
   end
   

