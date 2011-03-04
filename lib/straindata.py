class StrainData():
    def __init__(self):
        self.__time_series_data = None
        self.__frequency_series_array = []

    def __iter__(self):
    """
    provide functionality to iterate over the array of frequency elements
    """

    def __rmul__(self):
    """
    overload multiply for data objects, e.g. to allow multiplcation by a psd
    """

    def read_frames(self, channel_name, interval, cache_url):
    """
    channel_name: input gravitational wave strain channel name (e.g.  H1-LSC_STRAIN)
    interval: glue segment containing [start_time, end_time) of data to be read in
    cache_url: URL of a lal frame cache file

    This method fils self.__time_series_data with the data read in from the
    frame. It is responsible for allocating memory for the input data in a
    real_vector_t.
    """

    def to_float(self):
    """
    convert the data in self.__time_series_data to single precision
    """

    def fft_segments(self, segment_length, segment_overlap):
    """
    split the time_series data into segments and transform each segment into
    the frequency domain using fftw on the cpu.
    """

    def frequency_series_array_to_cuda_global_memory(self):
    """
    move the frequency series data from cpu memory to the gpu
    """

