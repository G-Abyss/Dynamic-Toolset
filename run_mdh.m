clear
clc
%%
gravity_vec = [0; 0; 9.80665];

rbt_config = Demo_MDH;
rbt = struct();
rbt.rbt_df = DefineRobot_MDH('Demo_MDH',rbt_config);
rbt.geom = GeometryCalculation_MDH(rbt.rbt_df);
rbt.dyn = Dynamics_MDH(rbt.rbt_df, rbt.geom, gravity_vec);
rbt.base = DynamicBaseParamCalc_MDH(rbt.rbt_df, rbt.dyn);


%% save output
[filepath,name,ext] = fileparts(mfilename("fullpath"));
out_path = fullfile(filepath,"output");
rbt.rbt_df.out_path = out_path;
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