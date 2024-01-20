function [loadfile] = problem_select(i)
if i==1
    loadfile='control1.dat-s';
elseif i==2
    loadfile='control2.dat-s';
elseif i==3
    loadfile='control3.dat-s';
elseif i==4
    loadfile='truss1.dat-s';
elseif i==5
    loadfile='truss2.dat-s';
elseif i==6
    loadfile='truss3.dat-s';
end

end

