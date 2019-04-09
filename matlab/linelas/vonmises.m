% genearl formula for von miss's effective stress
%
%   Inputs:
%       s   1x6     stress tensor in Voigt notation
%
function vm = vonmises(s)
    vm = sqrt((s(1)-s(2))^2 + (s(2)-s(3))^2 + (s(3)-s(1))^2 + ...
        6*(s(4)^2+s(5)^2+s(6)^2)) / sqrt(2);
end
% function vm = vonmises(s)
%     st = [
%         s(1) s(6) s(5)
%         s(6) s(2) s(4)
%         s(5) s(4) s(3)
%     ];
%     vm = sqrt(sum(st.^2,'all'));
% end



