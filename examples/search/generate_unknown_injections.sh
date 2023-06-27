if [[ ! -f injections.hdf ]]
then
    echo -e "\\n\\n>> [`date`] Generating injections"

    rm -rf ./strain

    ./generate_injections.py
fi

start_time=1186740069
end_time=1186743653
if [[ ! -d strain ]]
then
    echo -e "\\n\\n>> [`date`] Generating simulated strain"

    function simulate_strain { # detector PSD_model random_seed
        mkdir -p strain/$1

        out_path="strain/$1/$1-SIMULATED_STRAIN-{start}-{duration}.gwf"

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
	    --fft-backends mkl


	echo "${1::1} ${1}_simulated 1186740069 3584 $pwd/strain/$1/$1-SIMULATED_STRAIN-$start_time-3584.gwf" >> strain/$1/$1-cache.lcf
    }
    simulate_strain H1 aLIGOMidLowSensitivityP1200087 1234
    simulate_strain L1 aLIGOMidLowSensitivityP1200087 2345
    simulate_strain V1 AdVEarlyLowSensitivityP1200087 3456

fi

