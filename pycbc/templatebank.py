class TemplateBank():
    def __init__(self, dt, n):
        pass

    def __iter__(self):
        pass

    def read(self, url):
        """
        Read templates from a ligolw xml file
        """
        pass


class TemplateBankCPU():
    def __init__(self):
        self.create_template_frequency_series()

    def create_template_frequency_series(self):
        """
        use malloc()
        """
        pass


class TemplateBankGPU():
    def __init__(self):
        self.create_template_frequency_series()

    def create_template_frequency_series(self):
        """
        use cudaMalloc()
        """
        pass
