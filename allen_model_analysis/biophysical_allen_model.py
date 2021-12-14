from allensdk.api.queries.biophysical_api import BiophysicalApi
from allensdk.model.biophys_sim.neuron.hoc_utils import HocUtils

class Utils(HocUtils):
    def __init__(self, description):
        super(Utils, self).__init__(description)

bp = BiophysicalApi()
bp.cache_stimulus = True # change to False to not download the large stimulus NWB file
neuronal_model_id = 515175253    # get this from the web site as above
bp.cache_data(neuronal_model_id, working_directory='neuronal_model')
description='all-active'
utils = Utils.create_utils(description)
h = utils.h