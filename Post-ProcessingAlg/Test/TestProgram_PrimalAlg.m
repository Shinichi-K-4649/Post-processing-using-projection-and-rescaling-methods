clear

%% Select a problem
num =3;
[loadfile] = problem_select(num);
[At,b,c,K] = fromsdpa(loadfile); A = At.';

%% Select an input solution
Solver_List = ["Mosek","SDPT3","SDPA"]; solver_num = 1;
Solver = Solver_List(solver_num);
Tol = "default";
loadfile = erase(loadfile,".dat-s");
loadfile = Solver+"_para="+Tol+"_"+loadfile+".mat";
load(loadfile);
if Solver == "Mosek"
    [psol,dsol] = Sol_Mos2Sed(prob,res);
    sol.x = psol; sol.y = dsol.y; sol.z = dsol.z;
end

%% Set parameters
para.Maxtime = 1800;
para.Tol = 1e-12;
para.PraTerCon = 30;

%% Start post-processing algorithm with primal model
fprintf('\n Starting problem #i= %i  from file   %s \n',num,loadfile)
[bound,TmpSol,PPinfo] = PostProcessingAlg(A,b,c,K,sol,para);
fprintf('\n Finish post-processing \n')

%% Comparison of output and input solutions
[dimacs_orig] = cal_dimacsError(A,b,c,K,sol);
fprintf('\n Dimacs Error %7s   %7s   %7s   %7s   %7s   %7s', ...
    "err1","err2","err3","err4","err5","err6")
fprintf('\n Input Sol    %3.3E %3.3E %3.3E %3.3E %3.3E %3.3E', ...
    dimacs_orig.err1, dimacs_orig.err2, dimacs_orig.err3, dimacs_orig.err4, dimacs_orig.err5, dimacs_orig.err6)
dimacs = cal_dimacsError(A,b,c,K,TmpSol);
fprintf('\n Tmp   Sol    %3.3E %3.3E %3.3E %3.3E %3.3E %3.3E', ...
    dimacs.err1, dimacs.err2, dimacs.err3, dimacs.err4, dimacs.err5, dimacs.err6)
fprintf('\n')