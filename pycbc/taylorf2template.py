class TaylorF2Template():
    def __init__(self):
        pass

    def condition_data(self, data):
        """
        multiply data by f^{-7/6} on the cpu

        XXX BADRI WRITE THIS FUNCTION XXX
        """
        pass


class TaylorF2TemplateCPU():
    def __init__(self, dt, n):
        """
        Create an array of x^{-1/3} in CPU memory
        """
        pass

    def generate_filter_waveform(self, mchirp, eta):
        """
        Generate the waveform that will be used in the matched filter
        """
        pass

    def generate_strain(self, mchirp, eta):
        """
        Generate h(t) strain data for this waveform.
        """
        pass
    

class TaylorF2TemplateGPU():
    def __init__(self, dt, n):
        """
        Create an array of x^{-1/3} in GPU memory
        """
        pass

    def generate_filter_waveform(self, mchirp, eta):
        """
        Generate the waveform that will be used in the matched filter
        """
        pass

    def generate_strain(self, mchirp, eta):
        """
        Generate h(t) strain data for this waveform.
        """
        pass
