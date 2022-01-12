import numpy as np
import matplotlib.pyplot as plt
from open_pickle import read_from_pickle
from add_figure import add_figure
import os
import pickle
from glob import glob
spine_type="mouse_spine"

initial_locs=glob('data/fit/'+spine_type+'/*')
for loc in initial_locs:
    datas=glob(loc+'/different_initial_conditions/RA0*/RA0_fit_results.p')
    for data in datas:
        if '+' in data: datas.remove(data)
    print(loc)
    print('datas', datas)
    data1,data2=datas
    #####change the locations
    for data in [data1,data2]:
        save_folder=data[:data.rfind('/')]+'/analysis'
        try:os.mkdir(save_folder)
        except FileExistsError: pass

        dict=read_from_pickle(data)
        dict1=read_from_pickle(data1)
        dict2=read_from_pickle(data2)

        RA0=np.sort([float(key.split('=')[-1]) for key in dict.keys()])
        value=[dict[key] for key in dict.keys()]
        errors=[value[i]['error'][1] for i in range(len(value))]
        RAs=[value[i]['RA'] for i in range(len(value))]
        RMs=[value[i]['RM'] for i in range(len(value))]
        CMs=[value[i]['CM'] for i in range(len(value))]


        add_figure('diffrent RA0 against error\n'+loc.split('/')[-1],'RA0','errors')
        plt.plot(RA0,errors)
        minimums_arg = np.argsort(errors)
        dict_minimums={}
        for mini in minimums_arg[:10]:
            plt.plot(RA0[mini], errors[mini], '*',
                     label='RA0=' + str(RA0[mini]) + ' RM=' + str(round(RMs[mini], 2)) + ' RA=' + str(
                         round(RAs[mini], 2)) + ' CM=' + str(
                         round(CMs[mini], 2)) + ' error=' +  str(round(errors[mini]*100, 3)))
            del value[mini]['error']
            dict_minimums['RA0=' + str(RA0[mini])]=value[mini]
        plt.legend(loc='upper left')
        plt.savefig(save_folder+'/diffrent RA0 against error.png')

        add_figure('diffrent RA against error\n'+loc.split('/')[-1],'RA','errors')
        plt.plot(RAs,errors,'.')
        for mini in minimums_arg[:10]:
            plt.plot(RAs[mini], errors[mini], '*',
                     label=' RM=' + str(round(RMs[mini], 2)) + ' RA=' + str(
                         round(RAs[mini], 2)) + ' CM=' + str(
                         round(CMs[mini], 2)) +' RA0=' + str(RA0[mini])+ ' error=' +  str(round(errors[mini]*100, 3)) )
        plt.legend(loc='upper left')
        plt.savefig(save_folder+'/diffrent RA against error.png')

        add_figure('diffrent RA0 against CM\n'+loc.split('/')[-1],'RA0','CM')
        plt.plot(RA0,CMs)
        plt.savefig(save_folder+'/diffrent RA0 against CM.png')
        add_figure('diffrent RA0 against RA after fit\n'+loc.split('/')[-1],'RA0','RA')
        plt.plot(RA0,RAs)
        plt.savefig(save_folder+'/diffrent RA0 against RA after fit.png')
        add_figure('diffrent RA0 against RM\n'+loc.split('/')[-1],'RA0','RM')
        plt.plot(RA0,RMs)
        plt.savefig(save_folder+'/diffrent RA0 against RM.png')
        pickle.dump(dict_minimums, open(save_folder + "/RA0_10_minimums.p", "wb"))
    if float(next(iter(dict1.keys())).split('=')[-1])<float(next(iter(dict2.keys())).split('=')[-1]):
        dict = dict1.copy()  # Copy the dict1 into the dict3 using copy() method
        for key, value in dict2.items():  # use for loop to iterate dict2 into the dict3 dictionary
            dict[key] = value
    else:
        dict = dict2.copy()  # Copy the dict1 into the dict3 using copy() method
        for key, value in dict1.items():  # use for loop to iterate dict2 into the dict3 dictionary
            dict[key] = value

    save_folder=loc+'/different_initial_conditions/'+data1.split('/')[-2]+'_'+data2.split('/')[-2]
    try:os.mkdir(save_folder)
    except FileExistsError: pass

    RA0=[float(key.split('=')[-1]) for key in dict.keys()]
    value=[dict[key] for key in dict.keys()]
    errors=[value[i]['error'][1] for i in range(len(value))]
    RAs=[value[i]['RA'] for i in range(len(value))]
    RMs=[value[i]['RM'] for i in range(len(value))]
    CMs=[value[i]['CM'] for i in range(len(value))]
    add_figure('diffrent RA0 against error\n'+loc.split('/')[-1],'RA0','errors')
    plt.plot(RA0,errors)
    minimums_arg=np.argsort(errors)
    # mini = np.argmin(errors)
    dict_minimums={}
    # print(loc)
    for mini in minimums_arg[:10]:
        plt.plot(RA0[mini], errors[mini], '*',label='RA0=' + str(RA0[mini]) + ' RM=' + str(round(RMs[mini], 2)) + ' RA=' + str(round(RAs[mini], 2)) + ' CM=' + str(
                     round(CMs[mini], 2)) + ' error=' + str(round(errors[mini]*100, 3)))
        dict_minimums['RA0=' + str(RA0[mini])]=dict.get('RA0=' + str(RA0[mini]), dict.get('RA0=' + str(int(RA0[mini])), None))
        if dict_minimums['RA0=' + str(RA0[mini])] is None:
            raise TypeError("dict_minimums get None")
    pickle.dump(dict_minimums, open(save_folder + "/RA0_10_minimums.p", "wb"))
    plt.legend(loc='upper left')
    plt.savefig(save_folder+'/diffrent RA0 against error.png')

    add_figure('diffrent RA against error\n'+loc.split('/')[-1],'RA','errors')
    plt.plot(RAs,errors,'.')
    for mini in minimums_arg[:10]:
        plt.plot(RAs[mini], errors[mini], '*',label= 'RM=' + str(round(RMs[mini], 2)) + ' RA=' + str(round(RAs[mini], 2)) + ' CM=' + str(
                     round(CMs[mini], 2)) + ' error=' + str(round(errors[mini]*100, 3)) )
    plt.legend(loc='upper left')
    plt.savefig(save_folder+'/diffrent RA against error.png')

    add_figure('diffrent RA0 against CM\n'+loc.split('/')[-1],'RA0','CM')
    plt.plot(RA0,CMs)
    plt.savefig(save_folder+'/diffrent RA0 against CM.png')
    add_figure('diffrent RA0 against ֵֻֻRA after fit\n'+loc.split('/')[-1],'RA0','RA')
    plt.plot(RA0,RAs)
    plt.savefig(save_folder+'/diffrent RA0 against RA after fit.png')
    add_figure('diffrent RA0 against RM\n'+loc.split('/')[-1],'RA0','RM')
    plt.plot(RA0,RMs)
    plt.savefig(save_folder+'/diffrent RA0 against RM.png')

    data=loc+'/const_param/RA/Ra_const_errors.p'
    save_folder1=data[:data.rfind('/')]+'/analysis'
    try:os.mkdir(save_folder1)
    except FileExistsError:pass
    dict3=read_from_pickle(data)
    RA0=dict3['RA']
    error=dict3['errors'][1]
    errors=dict3['errors']
    RAs=[value['RA'] for value in dict3['params']]
    RMs=[value['RM'] for value in dict3['params']]
    CMs=[value['CM'] for value in dict3['params']]
    add_figure('RA const against errors\n'+loc.split('/')[-1],'RA const','error')
    plt.plot(RA0,error)
    minimums_arg=np.argsort(error)
    dict_minimums2={}
    for mini in minimums_arg[:10]:
        plt.plot(RA0[mini], error[mini], '*',label=' RM=' + str(round(RMs[mini], 2)) + ' RA=' + str(round(RAs[mini], 2)) + ' CM=' + str(
                     round(CMs[mini], 2)) + ' error=' +  str(round(error[mini]*100, 3)))
        dict_minimums2['RA_const=' + str(RA0[mini])]={'params': {'RM': RMs[mini], 'RA': RAs[mini], 'CM': CMs[mini]},'error':[err[mini] for err in errors] }
    pickle.dump(dict_minimums2, open(save_folder1 + "/RA_const_10_minimums.p", "wb"))
    plt.legend(loc='upper left')
    plt.savefig(save_folder1+'/RA const against errors')
    plt.savefig(save_folder1+'/RA const against errors.pdf')

    end_plot=60
    add_figure('RA const against errors\n'+loc.split('/')[-1],'RA const','error')
    plt.plot(RA0[:end_plot],error[:end_plot])
    for mini in minimums_arg[:10]:
        plt.plot(RA0[mini], error[mini], '*',label=' RM=' + str(round(RMs[mini], 2)) + ' RA=' + str(round(RAs[mini], 2)) + ' CM=' + str(
                     round(CMs[mini], 2)) + ' error=' +  str(round(error[mini]*100, 3)))
    plt.legend(loc='upper left')
    plt.savefig(save_folder1+'/RA const against errors until point '+str(end_plot))

    add_figure('RA const against RMs\n'+loc.split('/')[-1],'RA const','RM')
    plt.plot(RA0,RMs)
    plt.savefig(save_folder1+'/RA const against RM')
    add_figure('RA const against RA after fit\n'+loc.split('/')[-1],'RA const','RA')
    plt.plot(RA0,RAs)
    plt.savefig(save_folder1+'/RA const against RA after fit')
    add_figure('RA const against CM\n'+loc.split('/')[-1],'RA const','CM')
    plt.plot(RA0,CMs)
    plt.savefig(save_folder1+'/RA const against CMs')
