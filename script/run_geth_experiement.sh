echo "PATH TO DOWNLOADED ETHEREUM DATA: $P"
RESULT=$(pwd)/result

# 2th & 4th state, iblt and cdma

for seed in $(seq 44 48); do
    echo '----'
    echo '{
        "root path": "#PWD",
        "seed": #sd,
        "k": 7,
        "d": 233000000,
        "counting": true,
        "A path": "#P/state-22020359-sha256-sorted.txt",
        "B path": "#P/state-22399992-sha256-sorted.txt",
        "Intersect path": "#P/state-22020359-intersect-22399992.txt",
        "resolving round": 0,
        "result filename": "#RESULT/geth_full_2nd_with_4th_sd#sd.json"
    }' \
    | sed "s/#sd/${seed}/g"\
    | sed "s/#P/${P}/g"\
    | sed "s/#RESULT/${RESULT}/g"\
    > $(pwd)/config/geth_cdma_2nd_with_4th_sha256_sd${seed}.json
    build/bin/two_party $(pwd)/config/geth_cdma_2nd_with_4th_sha256_sd${seed}.json
done

for seed in $(seq 44 48); do
    echo '----'
    echo '{
        "root path": "#PWD",
        "seed": #sd,
        "A path": "#P/state-22020359-sha256-sorted.txt",
        "B path": "#P/state-22399992-sha256-sorted.txt",
        "Intersect path": "#P/state-22020359-intersect-22399992.txt",
        "result filename": "#RESULT/geth_iblt_2nd_with_4th_sd#sd.json"
    }' \
    | sed "s/#sd/${seed}/g"\
    | sed "s/#P/${P}/g"\
    | sed "s/#RESULT/${RESULT}/g"\
    > $(pwd)/config/geth_iblt_2nd_with_4th_sha256_sd${seed}.json
    build/bin/two_party_iblt $(pwd)/config/geth_iblt_2nd_with_4th_sha256_sd${seed}.json
done

# 3rd & 4th state, iblt and cdma

for seed in $(seq 44 48); do
    echo '----'
    echo '{
        "root path": "#PWD",
        "seed": #sd,
        "k": 7,
        "d": 10000000,
        "counting": true,
        "A path": "#P/state-22392874-sha256-sorted.txt",
        "B path": "#P/state-22399992-sha256-sorted.txt",
        "Intersect path": "#P/state-22392874-intersect-22399992.txt",
        "resolving round": 0,
        "result filename": "#RESULT/geth_full_3rd_with_4th_sd#sd.json"
    }' \
    | sed "s/#sd/${seed}/g"\
    | sed "s/#P/${P}/g"\
    | sed "s/#RESULT/${RESULT}/g"\
    > $(pwd)/config/geth_cdma_3rd_with_4th_sha256_sd${seed}.json
    build/bin/two_party $(pwd)/config/geth_cdma_3rd_with_4th_sha256_sd${seed}.json
done

for seed in $(seq 44 48); do
    echo '----'
    echo '{
        "root path": "#PWD",
        "seed": #sd,
        "A path": "#P/state-22392874-sha256-sorted.txt",
        "B path": "#P/state-22399992-sha256-sorted.txt",
        "Intersect path": "#P/state-22392874-intersect-22399992.txt",
        "result filename": "#RESULT/geth_iblt_3rd_with_4th_sd#sd.json"
    }' \
    | sed "s/#sd/${seed}/g"\
    | sed "s/#P/${P}/g"\
    | sed "s/#RESULT/${RESULT}/g"\
    > $(pwd)/config/geth_iblt_3rd_with_4th_sha256_sd${seed}.json
    build/bin/two_party_iblt $(pwd)/config/geth_iblt_3rd_with_4th_sha256_sd${seed}.json
done
