    cur_time = datetime.datetime.now()

    #construct the argument parse and parse the arguments
    ap = argparse.ArgumentParser()
    ap.add_argument("-u", "--usertag", required=False, default=cur_time,
        help="label for given run")
    ap.add_argument("-o", "--output-dir", required=False,
        help="path to output directory")
    ap.add_argument("-n", "--normalize", required=False, default=True,
        help="normalize the energy of the output")
    ap.add_argument("-s", "--samp-freq", required=True, type=float,
        help="Sampling frequency of channel")

    args = ap.parse_args()


    #Initialize parameters
    out_dir = args.output_dir
    now = args.usertag
    os.makedirs('%s/run_%s' % (out_dir,now))  # Fail early if the dir already exists
    normalized = args.normalize # Set this as needed
    sampling = args.samp_freq #sampling frequency
    mismatch=.2
    qrange=(4,64)
    frange=(0,np.inf)

    # Read data and remove low frequency content
    fname = 'H-H1_LOSC_4_V2-1126259446-32.gwf'
    url = "https://losc.ligo.org/s/events/GW150914/" + fname
    urllib.urlretrieve(url, filename=fname)
    h1 = read_frame('H-H1_LOSC_4_V2-1126259446-32.gwf', 'H1:LOSC-STRAIN')
    h1 = TimeSeries(np.random.normal(size=64*4096), delta_t = 1. / sampling)
    h1 = highpass_fir(h1, 15, 8)

    # Calculate the noise spectrum
    psd = interpolate(welch(h1), 1.0 / 32)

    #perform Q-tiling
    Qbase, frange = qtiling(h1, qrange, frange, sampling, normalized, mismatch)

    #Choose Q-plane and plot
    qplane = Qplane(Qbase, h1, sampling, normalized, out_dir, now, frange)

    #Plot spectrogram
    plotter(qplane, out_dir, now, frange, h1, sampling)

    print 'Done!'

