class SingleDetectorEvent:
    def __init__(self):
        self.__event_array = []

    def cluster_events(self):
    """
    Perform some clustering on the events
    """
        pass

    def write(self, path):
    """
    Write the event array to disk as ligolw_xml
    """
        pass


class SingleDetectorEventCPU:
    def __init__(self):
        pass

    def find_events(self, htilde, matched_filter, chisq_veto):
        pass

class SingleDetectorEventGPU:
    def __init__(self):
        pass

    def find_events(self, htilde, matched_filter, chisq_veto):
        pass
