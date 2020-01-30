function [paddata] = padmat(matdata,L,val,dim)
%% padmat: pad a matrix with specified value on both sides
%
%   INPUT:
%       data     	:   input data (matrix)
%       L           :   pad length on each side
%       val         :   pad value
%       dim         :   dimension
%   OUTPUT:
%       paddata     :   padded data
%

if dim==1
    paddata = [val*ones(L(1),1) ; matdata ; val*ones(L(2),1)];
elseif dim==2
    paddata = [val*ones(1,L(1)) , matdata , val*ones(1,L(2))];
else
    error('Must be along 1st or 2nd dimension')
end

end