from foton import Filter

def filter_foton(data, frame_filenames, filter_file, swstat_channel_name, bits, filter_name):
    '''
    A naive function to determine if the filter was on at the time
    and then filter the data.

    This module just checks the first time in the SWSTAT channel
    to see if the filter was on, it doesn't check beyond that.

    This is just for a first test on a small chunck of data.
    '''

    # loop over bits that state if the filter was on or off
#    swstat = frame.read_frame(frame_filenames, swstat_channel_name,
#                      start_time=start_time, end_time=end_time)
#    bits = bin(int(swstat[0]))[2:12]

    for i in range(10):

        print bits

        # if bit is on then filter the data
        bit = int(bits[i])
        if bit:
            print filter_name, Filter(filter_file[filter_name][i]).design.name, i, bit
            print Filter(filter_file[filter_name][i]).sections
            data = Filter(filter_file[filter_name][i]).apply(data)

    return  data


