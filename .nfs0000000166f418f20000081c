from split_data import split2phenomena
from correct_noise import clear_noise
from open_pickle import bring_phenomena
from open_pickle import read_phenomena
from find_Rinput import find_Rinput
from check_dynamics import check_dynamics
if __name__=='__main__':
    #### creat the data
    folder_ = '/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis'
    split2phenomena(folder_)
    #### clear the data
    folder_data = folder_+'/data/'
    correct_syn,mean_syn = read_phenomena('syn', folder_data, 1)
    correct_short_pulse,mean_short_pulse = read_phenomena('short_pulse', folder_data, 1)
    correct_spike, mean_spike= read_phenomena('spike', folder_data, 1)

    #### find the fit for Rm
    folder_iv_curves = folder_data+'traces_img/2017_05_08_A_0006/'
    R, Rinput_list = find_Rinput(folder_)

    #check if the model is dynamic or not
    folder_data_short_pulse = folder_data+'traces_img/2017_05_08_A_0006/short_pulse/'
    check_dynamics(mean_short_pulse,folder_data_short_pulse)


