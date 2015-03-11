#!/bin/sh

function run_pycbc {
	DATA=$1
	SCHEME=$2
	TASKSET=$3
	PROFILE=$4
	CUSTOM_BANK=$5

	case ${DATA} in
	clean)
		GPS_START_TIME=970910772
		;;
	loud)
		GPS_START_TIME=962784401
		;;
	grumbly)
		GPS_START_TIME=969658210
		;;
	*)
		echo "Data must be one of clean, loud, or grumbly"
		exit 0
	esac


	CMD="taskset -c ${TASKSET} python "

	if [ "${PROFILE}" == 'python' ]
	then
		CMD="${CMD} -m cProfile -o profile.log "
	fi

	if [ "${PROFILE}" == 'cuda' ]
	then
		export COMPUTE_PROFILE=1
	fi

	if [ "${PROFILE}" == 'valgrind' ]
	then
		CMD="taskset -c ${TASKSET} valgrind --tool=callgrind python "
	fi

	if [ "${PROFILE}" == 'operf' ]
	then
		CMD="taskset -c ${TASKSET} operf python "
	fi

	BANK=TMPLTBANK.xml.gz

	if [ "${CUSTOM_BANK}" != "" ]
	then
		BANK=${CUSTOM_BANK}
	fi

	GPS_END_TIME=$(( $GPS_START_TIME + 2048 ))
	INSTALLED=`which pycbc_inspiral`

	CMD="${CMD} ${INSTALLED} "
	CMD="${CMD} --cluster-method window "
	CMD="${CMD} --cluster-window 1 "
	CMD="${CMD} --bank-file $BANK "
	CMD="${CMD} --approximant SPAtmplt "
	CMD="${CMD} --gps-start-time $GPS_START_TIME "
	CMD="${CMD} --gps-end-time   $GPS_END_TIME "
	CMD="${CMD} --snr-threshold 5.5 "
	CMD="${CMD} --strain-high-pass 25.0 "
	CMD="${CMD} --chisq-bins 16 "
	CMD="${CMD} --psd-inverse-length 16 "
	CMD="${CMD} --psd-segment-stride 128 "
	CMD="${CMD} --psd-segment-length 256 "
	CMD="${CMD} --psd-estimation median "
	CMD="${CMD} --segment-length 256 "
	CMD="${CMD} --segment-start-pad 112 "
	CMD="${CMD} --segment-end-pad 16 "
	CMD="${CMD} --low-frequency-cutoff 30.0 "
	CMD="${CMD} --pad-data 8 "
	CMD="${CMD} --sample-rate 4096 "
	CMD="${CMD} --order 7 "
	CMD="${CMD} --frame-files ${DATA}/*.gwf "
	CMD="${CMD} --channel-name H1:LDAS-STRAIN "
	CMD="${CMD} --output /usr1/${USER}/profiling_results/${DATA}/test_${TASKSET}.hdf5  "
	CMD="${CMD} --processing-scheme ${SCHEME} "

	# Might also need arguments like:
	#--fft-backends fftw \
	#--fftw-measure-level 0 \
	#--fftw-threads-backend openmp \
	#--fftw-input-float-wisdom-file ${WISDOM}

	$CMD &
}


function verify_args {
	for key in cluster-method cluster-window bank-file approximant gps-start-time gps-end-time snr-threshold strain-high-pass chisq-bins psd-inverse-length psd-segment-stride psd-segment-length psd-estimation segment-length egment-start-pad segment-end-pad low-frequency-cutoff pad-data sample-rate order frame-files channel-name output-format processing-scheme
	do
		if [ "${args["$key"]}" == "" ]
		then
			echo "Configuration is missing required argument $key"
			exit 0
		fi
	done
}


function parse_args {
	in_f=$1

	for line in `cat $in_f | grep -v '^#' | sed 's+ +QZX+g'`
	do
		l2=`echo $line | sed 's+QZX+ +g'`
		name=`echo $l2 | cut -f1 -d= | sed -e 's/^ *//' -e 's/ *$//'`
		value=`echo $l2 | cut -f2- -d= | sed -e 's/^ *//' -e 's/ *$//'`

		args["${name}"]="${value}"
	done
	
	for line in $*
	do
		echo ${line} | grep -q =
		if [ $? == 0 ]
		then
			l2=`echo $line | sed 's+QZX+ +g'`
			name=`echo $l2 | cut -f1 -d= | sed -e 's/^ *//' -e 's/ *$//'`
			value=`echo $l2 | cut -f2- -d= | sed -e 's/^ *//' -e 's/ *$//'`
			args["${name}"]="${value}"
		fi

	done
}


function get_schemes {
	num_cuda=$( [ "$1" == "" ] && echo 0 || echo $1 )
	mkl_fill=$( [ "$2" == "" ] && echo 0 || echo $2 )


	schemes=''

	cores_per_socket=`lscpu | grep Core\(s\) | cut -f2 -d:  | sed 's/ //g'`
	sockets=`lscpu | grep Socket\(s\) | cut -f2 -d:  | sed 's/ //g'`
	cuda_every=-1
	step_count=-1

	if [ $num_cuda != 0 ]
	then
		cuda_every=$(( $cores_per_socket * $sockets / $num_cuda ))
		step_count=0
	fi

	socket=0
	core=0
	cpu=0

	while [ $socket -lt $sockets ]
	do
		core=0

		while [ $core -lt $cores_per_socket ]
		do
			if [ $step_count == 0 ]
			then
				schemes="$schemes cuda|${cpu}"
			elif [ $mkl_fill != 0 ]
			then
				schemes="$schemes mkl:1|${cpu}"
			fi

			if [ $num_cuda != 0 ]
			then
				step_count=$(( ($step_count + 1) % $cuda_every ))
			fi

			core=$(( $core + 1 ))
			cpu=$(( $cpu + 1 ))
		done

		socket=$(( $socket + 1 ))
	done

	echo $schemes
}


function run_tests {
	data=$1
	schemes=$2
	outfile=$3
	profile=$4

	# Regen function caches
	rm -rf /usr1/${USER}/*
	mkdir -p /usr1/${USER}/profiling_results/${data}
	run_pycbc ${data} mkl:1 1 no TMPLTBANK_SMALL.xml.gz 
	wait
	run_pycbc ${data} cuda  1 no TMPLTBANK_SMALL.xml.gz 
	wait
	rm -f /usr1/${USER}/profiling_results/${data}/*

	start=`date +%s`

	for scheme_and_taskset in `echo $schemes`
	do
		scheme=`echo $scheme_and_taskset | cut -f1 -d\|`
		taskset=`echo $scheme_and_taskset | cut -f2 -d\|`

		run_pycbc $data $scheme $taskset $profile
		profile=''
	done

	wait

	for file in /usr1/${USER}/profiling_results/${data}/*hdf5
	do
		end=`stat -c%Z ${file}`
		echo $(( $end - $start )) >> $outfile
	done
}

