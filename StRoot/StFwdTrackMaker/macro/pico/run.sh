
cmd="${1:-NONE}"
label=${2:-muon_minus}
data_prefix="${3:-single_}"
tt="${4:-0}" # default to 0, can be set to 1 or 2 for different track types (used for QA)
run_eff() {
    echo root -l -b -q 'qa.C( "/Users/brandenburg.89/star/ssw/data/'"${label}"'/'"${data_prefix}"'", 0, "'"${label}"'" )'  
    root -l -b -q 'qa.C( "/Users/brandenburg.89/star/ssw/data/'"${label}"'/'"${data_prefix}"'", 0, "'"${label}"'" )'  
    root -l -b -q 'qa.C( "/Users/brandenburg.89/star/ssw/data/'"${label}"'/'"${data_prefix}"'", 1, "'"${label}"'" )'  
    root -l -b -q 'qa.C( "/Users/brandenburg.89/star/ssw/data/'"${label}"'/'"${data_prefix}"'", 2, "'"${label}"'" )'  
}

run_compare_tt() {
    root -l -b -q 'compare_tt.C( "'"${label}"'" )' 
}

if [ "$cmd" = "qa" ]; then
    root -l -b -q 'qa.C( "/Users/brandenburg.89/star/ssw/data/'"${label}"'/'"${data_prefix}"'", '${tt}', "'"${label}"'" )'
elif [ "$cmd" = "eff" ]; then
    run_eff "$label"
elif [ "$cmd" = "compare_tt" ]; then
    run_compare_tt "$label"
elif [ "$cmd" = "all" ]; then
    run_eff "$label"
    run_compare_tt "$label"
else
    cat <<EOF
Usage: $0 [cmd] [label] [data_prefix] [track_type]

  cmd          = qa | eff | compare_tt | all   (default: all)
  label        = dataset label (e.g. muon_minus, muon_plus) (default: muon_minus)
  data_prefix  = prefix for data directory (default: single_)
  track_type   = 0 (all), 1 (primaries), 2 (secondaries) — used only with cmd=qa (default: 0)

Examples:
  $0 eff muon_plus
  $0 qa muon_minus single_ 2
  $0 all muon_plus altprefix_

EOF
fi

