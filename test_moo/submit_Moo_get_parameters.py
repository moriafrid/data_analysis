from open_pickle import read_from_pickle
import os
import sys
dataset_jobs_folder='/running'
spine_type='mouse_spine'
in_parallel=False

shrinkage_by=1.2#round(1.0/0.7,4)
resize_diam_by=1
# in_parallel = sys.argv[1]
for resize_diam_by in [1,1.2,round(1.0/0.7,4)]:
    for shrinkage_by in [1,1.2, round(1.0/0.7,4)]:
        if shrinkage_by!=1 and resize_diam_by!=1:
            name2run='dend*'+str(round(resize_diam_by,2))+' &shrinkage by '+str(round(shrinkage_by,2))
        elif shrinkage_by!=1:
            name2run='F_shrinkage='+str(round(shrinkage_by,2))
        elif resize_diam_by!=1:
            name2run='dend*'+str(round(resize_diam_by,2))+' &shrinkage by '+str(round(shrinkage_by,2))
        else:
            name2run='no change'

        passive_val_name='RA_initial' #or RA_const
        try:dict_1=read_from_pickle('../data/fit/'+spine_type+'/'+name2run+'/different_initial_conditions/RA0_100:300:2_RA0_50:100:0.5/RA0_10_minimums.p')
        except:continue

        for i,key in enumerate(dict_1.keys()):
            if i !=1: continue
            if in_parallel:
                command="sbatch -p ss.q runs_change_passive_val_parallel.sh"
                send_command = " ".join([command, '30',str(dict_1[key]['RM']),str(round(dict_1[key]['RA'],4)),str(dict_1[key]['CM']),str(shrinkage_by),str(resize_diam_by),passive_val_name])
            else:
                command="sbatch -p ss.q runs_change_passive_val.sh"
                send_command = " ".join([command,"1",str(dict_1[key]['RM']),str(round(dict_1[key]['RA'],4)),str(dict_1[key]['CM']),str(shrinkage_by),str(resize_diam_by),passive_val_name])
            print(send_command)
            os.system(send_command)
        #
        passive_val_name='RA_const'
        dict_1=read_from_pickle('../data/fit/'+spine_type+'/'+name2run+'/const_param/RA/analysis/RA_const_10_minimums.p')
        for i,key in enumerate(dict_1.keys()):
            if i !=1: continue
            if in_parallel:
                command="sbatch -p ss.q runs_change_passive_val_parallel.sh"
                send_command = " ".join([command, '30',str(dict_1[key]['params']['RM']),str(round(dict_1[key]['params']['RA'],4)),str(dict_1[key]['params']['CM']),str(shrinkage_by),str(resize_diam_by),passive_val_name])
            else:
                command="sbatch -p ss.q runs_change_passive_val.sh"
                send_command = " ".join([command,"1",str(dict_1[key]['params']['RM']),str(round(dict_1[key]['params']['RA'],4)),str(dict_1[key]['params']['CM']),str(shrinkage_by),str(resize_diam_by),passive_val_name])

            print(send_command)
            os.system(send_command)
        #
