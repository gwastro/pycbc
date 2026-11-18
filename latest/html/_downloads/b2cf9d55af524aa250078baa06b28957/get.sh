set -e

download_if_absent() {
    local URL="$1"
    local FILENAME=$(basename "$URL")
    if [ ! -f "$FILENAME" ]; then
        echo "Downloading $FILENAME"
        curl -O -L --show-error --silent "$URL"
    else
        echo "File $FILENAME already exists, download skipped"
    fi
}

for channel in A E T
do
    strain_file=${channel}_TDI_v2.gwf
    download_if_absent https://zenodo.org/record/7497853/files/${strain_file}

    psd_file=${channel}_psd.txt
    download_if_absent https://zenodo.org/record/7497853/files/${psd_file}
done

params_file=MBHB_params_v2_LISA_frame.pkl 
download_if_absent https://zenodo.org/record/7497853/files/${params_file}
