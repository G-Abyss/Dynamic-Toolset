clear
clc
%%
gravity_vec = [0, 0, -9.8015];

rbt_config = Demo7DOF;
rbt = struct();
rbt.rbt_df = DefineRobot('Demo7DOF',rbt_config);
rbt.geom = GeometryCalculation(rbt.rbt_df);
rbt.dyn = Dynamics(rbt.rbt_df, rbt.geom, gravity_vec);
rbt.base = DynamicBaseParamCalc(rbt.rbt_df, rbt.dyn);


%% save output
[filepath,name,ext] = fileparts(mfilename("fullpath"));
out_path = fullfile(filepath,"output");
rbt_data_path = fullfile(out_path,rbt.rbt_df.name+"_rbt.mat");
[status, msg, msgID] = mkdir(out_path);
addpath(out_path);

save(rbt_data_path,"rbt");

%% H matrix for prime parameter
make_H_function(rbt);

%% Y matrix for base parameter 
make_Y_function(rbt);

%% MCG matrix : M(q, P)  C(q, dq, P)  G(q, P) 
make_MCG_function(rbt); 