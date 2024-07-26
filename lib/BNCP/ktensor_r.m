function tensor = ktensor_r( Factors, r )
%
% %[Author Notes]%
%   Author:        Zecan Yang
%   Email :        zecanyang@gmail.com
%   Release date:  April 08, 2022
%


N = length(Factors);
One_Factors = cell(N, 1);

for i=1:N
    One_Factors{i} = Factors{i}(:, r);
end

tensor =  double(ktensor(One_Factors));
end





