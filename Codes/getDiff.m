% Created by Jianqiu Zhang for calculating the differential of a vector XIC
% with length ll.
function diff=getDiff(XIC,ll,direction)
diff=zeros(ll,1);
if direction==-1
    for i=1:ll
        if i==1
           diff(i)=XIC(i);
        else
           diff(i)=XIC(i)-XIC(i-1);
        end   
    end 
else   
    for i=1:ll
        if i==ll
           diff(i)=XIC(i);
        else
           diff(i)=XIC(i+1)-XIC(i);
        end   
    end 
end   