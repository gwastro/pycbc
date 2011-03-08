class MatchedFilter:
    def __init__(self, data):
    """
    Initialize the matched filter for these data
    """
        pass


class MatchedFilterCPU:
    def __init__(self, data):
    """
    Initialize the matched filter for these data: q and qtilde
    """
        pass

    def generate_snr(self, stilde, htilde):
    """
    Fill this objects internal memory with \rho(t)
    """
        pass

    def max(self):
    """
    Find the maximum of q(t)
    """
        pass


class MatchedFilterGPU:
    def __init__(self, data):
    """
    Initialize the matched filter for these data: q and qtilde
    """
        pass

    def generate_snr(self, stilde, htilde):
    """
    Fill this objects internal memory with \rho(t)
    """
    pass

    def max(self):
    """
    Find the maximum of q(t)
    """
        pass
