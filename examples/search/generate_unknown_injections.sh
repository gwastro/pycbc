if [[ ! -f injections.hdf ]]
then
    echo -e "\\n\\n>> [`date`] Generating injections"

    rm -f ./*-SIMULATED_STRAIN-*.gwf

    ./generate_injections.py
fi

start_time=1186740069
end_time=1186743653

function simulate_strain { # detector PSD_model random_seed

    out_path="$1-SIMULATED_STRAIN-{start}-{duration}.gwf"
    if [[ ! -f $1-SIMULATED_STRAIN-1186740069-3584.gwf ]]
    then
        echo -e "\\n\\n>> [`date`] Generating simulated strain in $1"
        pycbc_condition_strain \
            --fake-strain $2 \
            --fake-strain-seed $3 \
            --output-strain-file $out_path \
            --gps-start-time $start_time \
            --gps-end-time $end_time \
            --sample-rate 16384 \
            --low-frequency-cutoff 15 \
            --channel-name $1:SIMULATED_STRAIN \
            --frame-duration 3584 \
            --injection-file injections.hdf \
            --fft-backends fftw
    fi
}

simulate_strain H1 aLIGOMidLowSensitivityP1200087 1234
simulate_strain L1 aLIGOMidLowSensitivityP1200087 2345
simulate_strain V1 AdVEarlyLowSensitivityP1200087 3456


