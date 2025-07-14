
cmd="${1:-all}"
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
    echo "Unknown command: $cmd"
    echo "Usage: $0 [muon_minus|muon_plus] [eff|compare_tt|all]"
fi

