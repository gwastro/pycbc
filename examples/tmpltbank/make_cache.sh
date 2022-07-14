mkdir -p cache
gw_data_find \
    --observatory H \
    --type H1_NINJA2_G1000176_EARLY_RECOLORED \
    --gps-start-time 900000024 \
    --gps-end-time 900010677 \
    --url-type file \
    --lal-cache > cache/H-H1_NINJA2_G1000176_EARLY_RECOLORED_CACHE-900000024-10653.lcf
