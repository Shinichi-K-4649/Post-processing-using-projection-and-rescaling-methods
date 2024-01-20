function [dimacs] = cal_dimacsError(A,b,c,K,sol)

x = sol.x; y = sol.y; z = sol.z;

dimacs.err1 = norm(A*x-b)/(1+max(abs(b)));

tmp = eigK(x,K); tmp = min(tmp);
tmp = -tmp/(1+max(abs(b))); tmp = [0,tmp];
dimacs.err2 = max(tmp);

dimacs.err3 = norm(z-c+A.'*y)/(1+max(abs(c)));

tmp = eigK(z,K); tmp = min(tmp);
tmp = -tmp/(1+max(abs(c))); tmp = [0,tmp];
dimacs.err4 = max(tmp);

dimacs.err5 = (c.'*x - b.'*y)/(1+abs(c.'*x)+abs(b.'*y));
dimacs.err6 = x.'*z/(1+abs(c.'*x)+abs(b.'*y));

end

