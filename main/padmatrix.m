function [paddata] = padmatrix(matdata,L,val,dim)
%% padmatrix: pad a matrix with specified values in set dimensions
%
%   INPUT:
%       matdata  	:   input data (matrix)
%       L           :   pad length on each side
%       val         :   pad value
%       dim         :   dimension (along this)
%
%   OUTPUT:
%       paddata     :   padded data
%

sz_mat = size(matdata);
sz_L = size(L);
if all(sz_L==1)
    L(2) = L(1);
    L(1) = 0;
end
mat1 = val*ones(sz_mat(dim), L(1));
mat2 = val*ones(sz_mat(dim), L(2));

if dim==2
    paddata = [mat1' ; matdata ; mat2'];
elseif dim==1
    paddata = [mat1 ; matdata ; mat2];
else
    error('Must be along 1st or 2nd dimension')
end

end