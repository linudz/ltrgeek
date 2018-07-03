#!/bin/bash -ue
# 1.1: Launch suffixerator. Creates an index that is needed for LTRharvest run

cp /Users/linudz/WorkFlows/LTRGeek/ltrgeek/peach_reduced.fa .
gt suffixerator -db /Users/linudz/WorkFlows/LTRGeek/ltrgeek/peach_reduced.fa -indexname /Users/linudz/WorkFlows/LTRGeek/ltrgeek/peach_reduced.fa -tis -suf -lcp -des -ssp -sds -dna

# 1.2: Launch LTRharvest with the parameters defined in the config file

echo "gt ltrharvest     -index /Users/linudz/WorkFlows/LTRGeek/ltrgeek/peach_reduced.fa     -seed 30      -minlenltr 100     -maxlenltr 1000     -mindistltr 1000     -maxdistltr 15000     -similar 85     -mintsd 4     -maxtsd 4     -vic 60     -overlaps best     -xdrop 5     -mat 2     -mis -2     -ins -3     -del -3 "
