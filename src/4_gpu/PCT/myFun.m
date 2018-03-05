function [D, E] = myFun(A, B, C)
    D = A.* B+ C;
    E = A + B .* C + 50.45;
end
