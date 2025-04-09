function func_handle = make_H_function(rbt)

mkdir(rbt.rbt_df.out_path);
tic;

out_path = rbt.rbt_df.out_path + "/";
H = rbt.base.H;
disp("Begin Creating H Function ...")
func_handle = matlabFunction(H, 'File',out_path+rbt.rbt_df.name+'_H_func', 'Vars', {rbt.rbt_df.coordinates', rbt.rbt_df.d_coordinates', rbt.rbt_df.dd_coordinates'});
disp('function writen into '+ out_path+rbt.rbt_df.name + '_H_func.m,time:'+string(toc))

end