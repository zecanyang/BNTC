function binary_tensor = sparse_miss(Is, rho, miss_method)
% sparse_miss: Missing tensor data
%
% %[Syntax]%: 
%   binary_tensor = sparse_miss(Is, rho, miss_method)
%
% %[Inputs]%:
%   Is            :       tensor size; vector
%   rho           :       missing data rate;
%   miss_method   :       miss method; string, 'rand', 'fiber', 'segment'
%
% %[Outputs]%:
%   binary_tensor :       size Is binary tenosr, element->0 missing data;
%
%
% %[Author Notes]%   
%   Author:        Zecan Yang
%   Email :        zecanyang@gmail.com
%   Release date:  Nov 03, 2021
%
%

binary_tensor = zeros(Is);
% randn('state', 5); rand('state', 5); %#ok<RAND>
randn('state', 2022); rand('state', 2022); %#ok<RAND>
switch miss_method
    case 'rand'
        omega = rand(Is)>=rho; % missing data index
        binary_tensor(omega) = 1;
    case 'fiber'
        omega = find(rand(Is(1)*Is(2), 1)>=rho); % missing data index
        Pi = zeros(Is(1), Is(2));
        Pi(omega) = 1;
        binary_tensor = repmat( Pi, [1 1 Is(3)]);
    otherwise
        disp('miss method err!!!')
end

end