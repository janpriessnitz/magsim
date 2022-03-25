
#!/bin/bash

mkdir -p data
rm commands.txt

for b in `seq 1 20`; do
    for d in `seq 1 20`; do
        echo "./run_config.sh ${d}0 $b" >> commands.txt
    done
done

cat commands.txt | parallel

rm commands.txt
