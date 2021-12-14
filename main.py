#from correct_noise import clear_noise
from clear_noises import clear_noise
from split_data import split2phenomena
from find_Rinput import find_Rinput
from check_dynamics import check_dynamics
from open_pickle import read_from_pickle

if __name__=='__main__':
    #### creat the data
    folder_ = '/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/'
    split2phenomena(folder_)

    #### clear the data
    folder_data = folder_+'data/'
    # mean_short_pulse, x_mean_short_pulse, correct_short_pulse, x_correct_short_pulse = clear_noise('short_pulse', folder_data,STDtime2clear=2.5)
    # clear_noise('syn', folder_data,STDtime2clear=1.8)
    # clear_noise('spike', folder_data)

    #correct_syn,mean_syn,hz = read_phenomena('syn', folder_data, 1)
    #print('beginning short_pulse',flush=True)
    #correct_short_pulse,mean_short_pulse,hz = read_phenomena('short_pulse', folder_data, 1)
    #correct_spike, mean_spike,hz= read_phenomena('spike', folder_data, 1)

    #### find the fit for Rm
    folder_iv_curves = folder_data+'traces_img/2017_05_08_A_0006/'
    # R, Rinput_list = find_Rinput(folder_iv_curves)

    #check if the model is dynamic or not
    # folder_data_check_dynamic = folder_data+'check_dynamic/'
    # mean_short_pulse,x_mean_short_pulse=read_from_pickle('/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/data/short_pulse/mean_short_pulse.p')
    # check_dynamics(mean_short_pulse,x_mean_short_pulse,folder_data_check_dynamic)

    #found the nueron property



