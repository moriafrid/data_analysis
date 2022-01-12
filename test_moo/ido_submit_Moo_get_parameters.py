from slurm_job import SlurmJobFactory
from open_pickle import read_from_pickle

dataset_jobs_folder='/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/test_moo/running'
job_factory = SlurmJobFactory(dataset_jobs_folder)
passive_val_name='RA_const' #or RA_initial
dict=read_from_pickle('/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/data/fit/different_initial_conditions/RA0_50:100:1+100:200:2/RA0_10_minimums.p')

for key in dict.keys():
    job_factory.send_job(f"job_{key}", f"python3 -u MOO_get_parameters.py --RM {dict[key]['RM']} --RA {dict[key]['RA']} --CM {dict[key]['CM']} --passive_val_name {passive_val_name}")

