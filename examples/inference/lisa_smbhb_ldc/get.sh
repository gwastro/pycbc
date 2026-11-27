set -e

download_if_absent() {
    local FILENAME=$1
    if [ ! -f "$FILENAME" ]; then
        # The original data lives under https://zenodo.org/record/7497853/files,
        # but that URL occasionally produces 403 errors on GitHub Actions, so
        # we use our own mirror on GitHub.
        local URL="https://github.com/gwastro/pycbc_data/raw/refs/heads/master/lisa/${FILENAME}"
        echo "Downloading $FILENAME"
        curl -O -L --show-error --silent "$URL"
    else
        echo "File $FILENAME already exists, download skipped"
    fi
}

for channel in A E T
do
    strain_file=${channel}_TDI_v2.gwf
    download_if_absent ${strain_file}

    psd_file=${channel}_psd.txt
    download_if_absent ${psd_file}
done

params_file=MBHB_params_v2_LISA_frame.pkl
download_if_absent ${params_file}
