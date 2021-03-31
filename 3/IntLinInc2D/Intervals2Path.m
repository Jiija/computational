function [P] = Intervals2Path(S)
%
% The function  Intervals2Path  generates a closed path  P  
% from the nonempty matrix  S  of boundary intervals. 
% The input matrix  S  should not have unbounded elements. 

    [P]=zeros(0,2);

    if size(S,1)==0
      return;
    end

    bp = [S(1,1) S(1,2)];
    P = bp;
    bs = bp;

    while size(S,1)>0

       for k=1:size(S,1);
%         if isequal(bs,[S(k,1) S(k,2)])
          if max(abs(bs-[S(k,1) S(k,2)]))<1.e-8
             i=k;
             break;
          end
       end
      
       es = [S(i,3) S(i,4)];
%      if ~isequal(es,bs)
       if max(abs(bs-es))>1.e-8
          P = [P; es];
%         if isequal(es,bp)
          if max(abs(bp-es))<1.e-8
             return;  % the path is closed
          end
          bs = es;
       end
       S(i,:)=[];

    end
    
end