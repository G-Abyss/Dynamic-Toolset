function func_handle = make_Y_function(rbt)

mkdir(rbt.rbt_df.out_path);
tic;
out_path = rbt.rbt_df.out_path + "/";
Y = rbt.base.H_b;
disp("Creating Y Function ...")
func_handle = matlabFunction(Y, 'File',out_path+rbt.rbt_df.name+'_Y_func', 'Vars', {rbt.rbt_df.coordinates', rbt.rbt_df.d_coordinates', rbt.rbt_df.dd_coordinates'});
disp('function writen into '+ out_path+rbt.rbt_df.name + '_Y_func.m,time:'+string(toc))

tic;
disp("Creating prime2beta Function ...")
matlabFunction(vpa(rbt.base.prime2beta), 'File',out_path+rbt.rbt_df.name+'_prime_param2beta_param');
disp('prime param to beta param mapping matrix writen into '+ out_path + rbt.rbt_df.name + '_prime_param2beta_param.m,time:'+string(toc))

end