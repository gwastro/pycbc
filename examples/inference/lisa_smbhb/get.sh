for channel in A E T
do
    strain_file=${channel}_TDI_v2.gwf
    test -f ${strain_file} && continue
    curl -O --silent https://zenodo.org/record/7078835/files/${strain_file}

    psd_file=${channel}_psd.txt
    test -f ${psd_file} && continue
    curl -O --silent https://zenodo.org/record/7078835/files/${psd_file}
done

params_file=MBHB_params_v2.pkl 
test -f ${params_file} && continue
curl -O --silent https://zenodo.org/record/7078835/files/${params_file}
