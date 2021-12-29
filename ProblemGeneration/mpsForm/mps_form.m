function mps_form(experiment_name, data_dir)

cd ~/gurobi911/linux64/matlab
gurobi_setup
cd ~/ProblemGeneration/mpsForm

experiment = [data_dir, '/', experiment_name, '/SVMtoLP'];
[nb,nf,mi,me,m,n,Ae,Ai,be,bi,c] = dataRead(experiment);

lb = [zeros(nb,1); -inf(nf,1)]';
ub = [];
% Build Gurobi model
model.modelsense = 'min';
model.obj = c';
model.A = [Ai; Ae]; % A must be sparse
model.sense = [repmat('<',size(Ai,1),1); repmat('=',size(Ae,1),1)]';
model.rhs = [bi; be]'; % rhs must be dense

if ~isempty(lb)
    model.lb = lb;
else
    model.lb = -inf(size(model.A,2),1)'; % default lb for MATLAB is -inf
end
if ~isempty(ub)
    model.ub = ub;
end

experiment_name_mps = [data_dir, '/', experiment_name,'.mps'];
gurobi_write(model, experiment_name_mps);

end
