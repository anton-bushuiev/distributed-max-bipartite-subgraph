# ============================Arguments=========================================
# 1 - path for input directory with results
# 2 - path for output file
# ==============================================================================

function exactly_one_new_line() 
{
    [[ $(wc -l "$1") == "1 $1" ]]
}

# Write csv header
echo "Submit datetime, Solver type, Nodes type, #CPUs, Queue type, Input, #Threads, Sort, Real time" \
> "$2"

# Write lines (files)
for f in results/final/*; do
  if exactly_one_new_line "$f"; then
    cat "$f" >> "$2"
  else
    cat "$f" | head -1 | cut -d, -f 1-8 | tr -d '\n' >> "$2"
    echo ", inf" >> "$2"
  fi
done

