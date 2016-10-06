function Matrix_w=Matrix_Weight_func(Matrix)

% figure
% plot(Matrix)
if sum(sum(Matrix,1))~=0
    Weight_vector_C=sum(Matrix,1)./sum(sum(Matrix,1));%%%%% Column weighted    
    Matrix_C=Matrix*diag(Weight_vector_C);
    Weight_vector_R=sum(Matrix_C,2)./sum(sum(Matrix_C,2));%%%%%% Row weighted
%     figure
%     plot(Weight_vector_R.^1.5)
    Weight_vector_R_n=Weight_vector_R.^0.75;
    Weight_vector_R_n=Weight_vector_R_n./sum(Weight_vector_R_n);
    Matrix_w=(diag(Weight_vector_R_n)*Matrix);
else
    Matrix_w=Matrix;    
end





