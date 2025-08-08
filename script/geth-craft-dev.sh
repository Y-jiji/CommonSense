pushd data/geth
head -n 10000 state-a.txt > state-a-dev.txt
head -n 10000 state-b.txt > state-b-dev.txt
omni diff-sorted-string --input-a state-a-dev.txt --input-b state-b-dev.txt --output-a-minus-b state-a-minus-b-dev.txt --output-b-minus-a state-b-minus-a-dev.txt --output-intersect state-intersect-dev.txt
for x in $(ls *-dev.txt)
do
    wc -l $x
done