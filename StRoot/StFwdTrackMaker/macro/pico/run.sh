
label=${1:-muon_minus}
cmd="${2:-all}"
run_eff() {
    root -l -b -q 'qa.C( "/Users/brandenburg.89/star/ssw/data/'"$1"'/single_", 0, "'"$1"'" )'  
    root -l -b -q 'qa.C( "/Users/brandenburg.89/star/ssw/data/'"$1"'/single_", 1, "'"$1"'" )'  
    root -l -b -q 'qa.C( "/Users/brandenburg.89/star/ssw/data/'"$1"'/single_", 2, "'"$1"'" )'  
}

run_compare_tt() {
    root -l -b -q 'compare_tt.C( "'"$1"'" )' 
}

if [ "$cmd" = "qa" ]; then
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

