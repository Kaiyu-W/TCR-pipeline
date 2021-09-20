#!/bin/bash

# run the command and capture the stdout
# $@ is a java -jar code, so both stdout and stderr will output into terminal
out=$(command $@)

echo "$out"

if [[ $(echo "$out" | grep -c reduce) -gt 0 ]]; then
    echo ""
    echo "Original Rscript failed, try reduce = 0 ..."
    echo "Modifying vj_pairing_plot.r ..."
    sed -i 's/mat, annotationTrack/mat, reduce = 0, annotationTrack/' vj_pairing_plot.r
    cmd=$(echo "$out" | grep "\\[RUtil\\] Executing" | sed 's/\[RUtil\] Executing //')
    echo "Run again: $cmd"
    
    # $cmd is a Rscript code, so you need to set it up manually that stdout and stderr should output into terminal, by 2>&1
    out2=$(command $cmd 2>&1)
    echo "$out2"

    if [[ $(echo "$out2" | grep Error | wc -l) -gt 0 ]]; then
        echo ""
        echo "Original Rscript failed, try gap.degree = 0.5 ..."
        echo "Modifying vj_pairing_plot.r ..."
        sed -i "/circos.par/s/rep(1/rep(0.5/g" vj_pairing_plot.r

        echo "Run again: $cmd"
        command $cmd
    fi
else
    echo ""
    echo "Not a 'reduce' failure for original script, nothing to do."
fi
