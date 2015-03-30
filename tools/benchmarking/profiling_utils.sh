#!/bin/sh

function run_pycbc {
	DATA=$1
	SCHEME=$2
	TASKSET=$3
	PROFILE=$4
	CUSTOM_BANK=$5

	case ${DATA} in
	clean)
		GPS_START_TIME=970910372
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

	DEVICE=""
	if  `echo ${SCHEME} | grep -q cuda:`
	then
		DEVICE=`echo ${SCHEME} | cut -f2 -d:`
		SCHEME=cuda
	fi

	CMD="taskset -c ${TASKSET} python "

	if [ "${PROFILE}" == 'python' ]
	then
		CMD="${CMD} -m cProfile -o profile_${args["tag"]}_${DATA}.log "
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
		echo_if_dryrun 'mkdir -p `pwd`/${args["tag"]}_${DATA}'
		CMD="taskset -c ${TASKSET} operf -d `pwd`/${args["tag"]}_${DATA}  python "
	fi

	BANK=${args["template-bank"]}

	if [ "${CUSTOM_BANK}" != "" ]
	then
		BANK=${CUSTOM_BANK}
	fi

	# Calculate end time
	segment_length=${args["segment-length"]}
	segment_start_pad=${args["segment-start-pad"]}
	segment_end_pad=${args["segment-end-pad"]}
	num_segments=${args["num-segments"]}

	GPS_END_TIME=$(( ${GPS_START_TIME}+(${segment_length}-${segment_start_pad}-${segment_end_pad})*${num_segments}+${segment_start_pad}+${segment_end_pad} ))
 
	INSTALLED=`which pycbc_inspiral`
	profile_root=`python -c 'import tempfile ; print tempfile.gettempdir()'`
	mkdir -p ${profile_root}/profiling_results/${DATA}

	CMD="${CMD} ${INSTALLED} "
	CMD="${CMD} --cluster-method ${args["cluster-method"]}"
	CMD="${CMD} --cluster-window ${args["cluster-window"]} "
	CMD="${CMD} --bank-file $BANK "
	CMD="${CMD} --approximant ${args["approximant"]} "
	CMD="${CMD} --gps-start-time $GPS_START_TIME "
	CMD="${CMD} --gps-end-time   $GPS_END_TIME "
	CMD="${CMD} --snr-threshold ${args["snr-threshold"]} "
	CMD="${CMD} --strain-high-pass ${args["strain-high-pass"]} "
	CMD="${CMD} --chisq-bins ${args["chisq-bins"]} "
	CMD="${CMD} --psd-inverse-length ${args["psd-inverse-length"]} "
	CMD="${CMD} --psd-segment-stride ${args["psd-segment-stride"]} "
	CMD="${CMD} --psd-segment-length ${args["psd-segment-length"]} "
	CMD="${CMD} --psd-estimation ${args["psd-estimation"]} "
	CMD="${CMD} --segment-length ${args["segment-length"]} "
	CMD="${CMD} --segment-start-pad ${args["segment-start-pad"]} "
	CMD="${CMD} --segment-end-pad ${args["segment-end-pad"]} "
	CMD="${CMD} --low-frequency-cutoff ${args["low-frequency-cutoff"]} "
	CMD="${CMD} --pad-data ${args["pad-data"]} "
	CMD="${CMD} --sample-rate ${args["sample-rate"]} "
	CMD="${CMD} --order ${args["order"]} "
	CMD="${CMD} --frame-files ${DATA}/*.gwf "
	CMD="${CMD} --channel-name ${args["channel-name"]} "
	CMD="${CMD} --output ${profile_root}/profiling_results/${DATA}/test_${TASKSET}.hdf5  "
	CMD="${CMD} --processing-scheme ${SCHEME} "
	CMD="${CMD} ${args["extra_args"]} "

	if [ "${DEVICE}" != "" ]
	then
		CMD="${CMD} --processing-device-id ${DEVICE} "
	fi

	if [ "${args["wisdom"]}" != "" ]
	then
		CMD="${CMD} --fft-backends fftw "
		CMD="${CMD} --fftw-measure-level 0 "
		CMD="${CMD} --fftw-threads-backend openmp "
		CMD="${CMD} --fftw-input-float-wisdom-file ${args["wisdom"]}"
	fi

	if [ "${args["dryrun"]}" != "" ]
	then
		echo "$CMD &"
		touch ${profile_root}/profiling_results/${DATA}/test_${TASKSET}.hdf5
		sleep 1 &
	else
		$CMD &
	fi
}


function verify_args {
	for key in cluster-method cluster-window approximant snr-threshold strain-high-pass chisq-bins psd-inverse-length psd-segment-stride psd-segment-length psd-estimation segment-length segment-start-pad segment-end-pad low-frequency-cutoff pad-data sample-rate order channel-name processing-schemes tag num-segments template-bank
	do
		if [ "${args["$key"]}" == "" ]
		then
			echo "Configuration is missing required argument $key"
			exit 0
		fi
	done
}

function parse_arg {
	line=$1

	if  `echo ${line} | grep -q =`
	then
		l2=`echo $line | sed 's+QZX+ +g'`
		name=`echo $l2 | cut -f1 -d= | sed -e 's/^ *//' -e 's/ *$//'`
		value=`echo $l2 | cut -f2- -d= | sed -e 's/^ *//' -e 's/ *$//'`

		if `echo ${name} | grep -q pycbc:`
		then
			name=`echo $name | sed 's+pycbc:++'`
			args['extra_args']="${args['extra_args']} --${name} ${value}"
		else
			args["${name}"]="${value}"
		fi
	fi
}
	
function parse_args {
	in_f=$1

	for line in `cat $in_f | grep -v '^#' | sed 's+ +QZX+g'`
	do
		parse_arg "${line}"
	done
	
	for line in $*
	do
		parse_arg "${line}"
	done
}

function one_per_core {
	type=$1
	schemes=''

	cores_per_socket=`lscpu | grep Core\(s\) | cut -f2 -d:  | sed 's/ //g'`
	sockets=`lscpu | grep Socket\(s\) | cut -f2 -d:  | sed 's/ //g'`

	socket=0
	core=0
	cpu=0

	while [ $socket -lt $sockets ]
	do
		core=0

		while [ $core -lt $cores_per_socket ]
		do
			schemes="$schemes ${type}:1|${cpu}"
			core=$(( $core + 1 ))
			cpu=$(( $cpu + 1 ))
		done

		socket=$(( $socket + 1 ))
	done

	echo $schemes
}

function one_per_socket {
	type=$1
	schemes=''

	cores_per_socket=`lscpu | grep Core\(s\) | cut -f2 -d:  | sed 's/ //g'`
	sockets=`lscpu | grep Socket\(s\) | cut -f2 -d:  | sed 's/ //g'`

	socket=0
	core=0
	cpu=0

	while [ $socket -lt $sockets ]
	do
		schemes="$schemes ${type}:${cores_per_socket}|${cpu}-$(( ${cpu} + ${cores_per_socket} - 1)) "
		cpu=$(( ${cpu} + ${cores_per_socket} ))
		socket=$(( $socket + 1 ))
	done

	echo $schemes
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


function echo_if_dryrun {
	if [ "${args["dryrun"]}" != "" ]
	then
		echo $1
	fi
}

function run_tests {
	data=$1
	schemes=$2
	outfile=$3
	profile=$4

	# Regen function caches, clean up previous runs
	cachedir=`python -c 'import os ; import sys ; import tempfile ; _python_name =  "python%d%d_compiled" % tuple(sys.version_info[:2]) ; _tmp_dir = tempfile.gettempdir() ; _cache_dir_name = repr(os.getuid()) + "_" + _python_name ; _cache_dir_path = os.path.join(_tmp_dir, _cache_dir_name); print _cache_dir_path' `

	profile_root=`python -c 'import tempfile ; print tempfile.gettempdir()'`

	rm -rf ${profile_root}/profiling_results
	mkdir -p ${profile_root}/profiling_results/${data}


	if [ "${args["regenerated_caches"]}" == "" ]
	then
		rm -rf ${cachedir}
		echo_if_dryrun '# Starting small runs to regenerate function caches'

		run_pycbc ${data} mkl:1 1 no TMPLTBANK_SMALL.xml.gz 
		wait
		echo_if_dryrun 'wait'
		run_pycbc ${data} cuda  1 no TMPLTBANK_SMALL.xml.gz 
		wait
		echo_if_dryrun 'wait'
		run_pycbc ${data} cpu 1 no TMPLTBANK_SMALL.xml.gz 
		wait
		echo_if_dryrun 'wait'

		rm -f ${profile_root}/profiling_results/${data}/*

		args["regenerated_caches"]="true"

	fi

	start=`date +%s`

	echo_if_dryrun ""
	echo_if_dryrun "# Starting test runs for $data"
	echo_if_dryrun ""

	if  `echo ${schemes} | grep -q one_per_core`
	then
		type=`echo ${schemes} | awk '{print $2}'`
		schemes=`one_per_core ${type}`
	fi

	if  `echo ${schemes} | grep -q one_per_socket`
	then
		type=`echo ${schemes} | awk '{print $2}'`
		schemes=`one_per_socket ${type}`
	fi

	if [ "${doing_profile}" != "" ]
	then
		profile=${doing_profile}
	fi

	doing_profile=""
	if [ "${profile}" != "" ]
	then
		echo_if_dryrun 'start=`date +%s`'
		doing_profile=${profile}
	fi

	for scheme_and_taskset in `echo $schemes`
	do
		scheme=`echo $scheme_and_taskset | cut -f1 -d\|`
		taskset=`echo $scheme_and_taskset | cut -f2 -d\|`

		run_pycbc $data $scheme $taskset $profile

		if [ "${profile}" != "" ]
		then
			output="${profile_root}/profiling_results/${data}/test_${taskset}.hdf5  "
			profile=''
		fi
	done

	echo_if_dryrun 'wait'
	wait

	if [ "${doing_profile}" != "" ]
	then
		echo_if_dryrun "end=\`stat -c%Z ${output}\`"
		echo_if_dryrun "echo \$(( \$end - \$start ))  >> ${args["tag"]}_runtimes.dat"
	fi


	rm -f ${outfile}

	for file in ${profile_root}/profiling_results/${data}/*hdf5
	do
		end=`stat -c%Z ${file}`
		echo $(( $end - $start )) >> $outfile
	done
}

