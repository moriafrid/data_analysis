from allensdk.model.biophys_sim.neuron.hoc_utils import HocUtils

class Utils(HocUtils):
    def __init__(self, description):
        super(Utils, self).__init__(description)

utils = Utils.create_utils(description)
h = utils.h