function  log_KL_Value=KL_calculate(Vector01,Vector02)

% Vector01=XIC_normal_Vector01;
% Vector02=XIC_normal_Vector0102;

%                     %%%%%%%%%%%% calculate KL based on non zeros position
%                     ID_notzeros01=find(Vector01~=0);
%                     ID_notzeros02=find(Vector02~=0);
%                     
%                     [C,IA,IB]=intersect(ID_notzeros01,ID_notzeros02);
%                     Vector01v1=Vector01(C)./sum(Vector01(C));
%                     Vector02v1=Vector02(C)./sum(Vector02(C));
%                     if ~isempty(C)
%                         log_KL_Value=log(sum(Vector01v1.*log(Vector01v1./Vector02v1)));
%                     else
%                         log_KL_Value=100;
%                     end
%                     %%%%%%%%%%%%%
                    
                    %%%%%%%%%%%% calculate KL based on add small values to zeros position
                    ID_zeros01=find(Vector01==0);
                    ID_zeros02=find(Vector02==0);
                    Vector01v1=Vector01;
                    Vector02v1=Vector02;
                    if ~isempty(ID_zeros01)
                        Vector01v1(ID_zeros01)=max(Vector01)/100000;
                    end
                    if ~isempty(ID_zeros02)
                        Vector02v1(ID_zeros02)=max(Vector02)/100000;
                    end
                    Vector01v2=Vector01v1./sum(Vector01v1);
                    Vector02v2=Vector02v1./sum(Vector02v1);
                    log_KL_Value=log2(sum(Vector01v2.*log(Vector01v2./Vector02v2)));
                    %%%%%%%%%%%%%
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    