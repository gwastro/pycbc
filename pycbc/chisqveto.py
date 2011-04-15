
#class ChisqVeto(pycbc.ChisqVeto, pycbc.ChisqVetoGPU):
#    def __init__(self):
#        pass


class ChisqVeto:
    def __init__(self, number_of_bins, inverse_power_spectrum):
        pass

class ChisqVetoCPU:
    def __init__(self, number_of_bins, inverse_power_spectrum):
        """
        create the memory needed for the chisq veto
        """
        pass

    def compute_frequency_bins(self, htilde):
        pass

    def generate_chisq(self, matched_filter ):
        pass

class ChisqVetoGPU():
    def __init__(self, number_of_bins, inverse_power_spectrum):
        """
        create the memory needed for the chisq veto
        """
        pass

    def compute_frequency_bins(self, htilde):
        pass

    def generate_chisq(self, matched_filter ):
        pass
