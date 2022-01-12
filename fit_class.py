from neuron import h, gui
import os
import pickle
from simulation import FitModel

# def get_inj(T,I,V):
#     #found the begining,end and median of the injection
#     I_abs = np.abs(I)
#     inj_start = np.where(I_abs > I_abs.max() / 4.0)[0][0] - 1
#     inj_end = np.where(I_abs > I_abs.max() / 4.0)[0][-1]
#     inj = np.median(I[inj_start:inj_end])
#     return inj, T[inj_start], T[inj_end]
def fit2short_pulse(fit_simulation_data,folder="",CM=1,RM=10000,RA=100):
    try:os.mkdir(folder)
    except FileExistsError: pass
    opt_vals = h.Vector(3)
    opt_vals.x[RM_IX] =RM
    opt_vals.x[RA_IX] = RA
    opt_vals.x[CM_IX] = CM
    fit_simulation_data.change_model_pas(CM=CM, RA=RA, RM=RM, E_PAS=fit_simulation_data.E_PAS)

    imp = h.Impedance(sec=fit_simulation_data.soma)
    imp.loc(fit_simulation_data.soma(0.5))
    imp.compute(0)
    imp.input(0)

    fit_simulation_data.plot_res(CM=CM, RM=RM, RA=RA, save_name="before")
    print('the initial impadence is', imp.input(0))
    # allway run the fitting 3 time to avoid stack in local minima
    for i in range(3):
        RMSD = h.fit_praxis(fit_simulation_data.efun,opt_vals)   #@# take too much time if the fitting isn't found
        RM = opt_vals.x[RM_IX]
        RA = opt_vals.x[RA_IX]
        CM = opt_vals.x[CM_IX]

        print("RMSD", RMSD,", RM",  RM, ", RA",  RA, ", CM",  CM)
        fit_simulation_data.plot_res(CM=CM, RM=RM, RA=RA, save_name="_fit_after_" + str(i + 1))
        imp.compute(0)
        print('the impadence is',imp.input(0))

    pickle.dump({
        "RM": RM,
        "RA": RA,
        "CM": CM
    }, open(folder + '/passive_parameters.p', "wb"))
    return {"RMSD": RMSD, "RM":RM, "RA":RA, "CM":  CM}


if __name__=='__main__':
    # curr_params.CM_IX= 1
    initial_folder = "data/fit/diffrent_initial_conditions"
    try:os.mkdir(initial_folder)
    except FileExistsError: pass

    cell_file = "05_08_A_01062017_Splice_shrink_FINISHED_LABEL_Bluecell_spinec91.ASC"
    short_pulse_file="data/short_pulse/clear_short_pulse.p"
    fit_simulation_data = FitModel()

    print("free params:")
    h.attr_praxis(1e-9, 1000, 0)

    cm_folder = initial_folder + "/CM"
    try:os.mkdir(cm_folder)
    except FileExistsError: pass

    CM_IX = 0
    RM_IX = 1
    RA_IX = 2

    cm=1
    try: os.mkdir(cm_folder)
    except FileExistsError:pass
    solution=fit2short_pulse(fit_simulation_data,folder=cm_folder+'/CM=1',CM=cm)


