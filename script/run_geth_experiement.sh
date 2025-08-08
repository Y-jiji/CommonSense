# 2th & 4th state, iblt and cdma

for seed in $(seq 44 48); do
    echo '----'
    echo '{
        "root path": "/media/gtnetuser/Dell/CDMA",
        "seed": #sd,
        "k": 7,
        "d": 233000000,
        "counting": true,
        "A path": "/media/gtnetuser/HDD_16TB_ATLAS/geth-2nd-state/state-22020359-sha256-sorted.txt",
        "B path": "/media/gtnetuser/HDD_16TB_ATLAS/geth-4th-state/state-22399992-sha256-sorted.txt",
        "Intersect path": "/media/gtnetuser/HDD_16TB_ATLAS/state-22020359-intersect-22399992.txt",
        "resolving round": 0,
        "result filename": "/media/mydrive/CDMA/result/geth_full_2nd_with_4th_sd#sd.json"
    }' \
    | sed "s/#sd/${seed}/g"\
    > /media/mydrive/CDMA/config/geth_cdma_2nd_with_4th_sha256_sd${seed}.json
    build/bin/two_party /media/mydrive/CDMA/config/geth_cdma_2nd_with_4th_sha256_sd${seed}.json
done

for seed in $(seq 44 48); do
    echo '----'
    echo '{
        "root path": "/media/gtnetuser/Dell/CDMA",
        "seed": #sd,
        "A path": "/media/gtnetuser/HDD_16TB_ATLAS/geth-2nd-state/state-22020359-sha256-sorted.txt",
        "B path": "/media/gtnetuser/HDD_16TB_ATLAS/geth-4th-state/state-22399992-sha256-sorted.txt",
        "Intersect path": "/media/gtnetuser/HDD_16TB_ATLAS/state-22020359-intersect-22399992.txt",
        "result filename": "/media/mydrive/CDMA/result/geth_iblt_2nd_with_4th_sd#sd.json"
    }' \
    | sed "s/#sd/${seed}/g"\
    > /media/mydrive/CDMA/config/geth_iblt_2nd_with_4th_sha256_sd${seed}.json
    build/bin/two_party_iblt /media/mydrive/CDMA/config/geth_iblt_2nd_with_4th_sha256_sd${seed}.json
done

# 3rd & 4th state, iblt and cdma

for seed in $(seq 44 48); do
    echo '----'
    echo '{
        "root path": "/media/gtnetuser/Dell/CDMA",
        "seed": #sd,
        "k": 7,
        "d": 10000000,
        "counting": true,
        "A path": "/media/gtnetuser/HDD_16TB_ATLAS/geth-3rd-state/state-22392874-sha256-sorted.txt",
        "B path": "/media/gtnetuser/HDD_16TB_ATLAS/geth-4th-state/state-22399992-sha256-sorted.txt",
        "Intersect path": "/media/gtnetuser/HDD_16TB_ATLAS/state-22392874-intersect-22399992.txt",
        "resolving round": 0,
        "result filename": "/media/mydrive/CDMA/result/geth_full_3rd_with_4th_sd#sd.json"
    }' \
    | sed "s/#sd/${seed}/g"\
    > /media/mydrive/CDMA/config/geth_cdma_3rd_with_4th_sha256_sd${seed}.json
    build/bin/two_party /media/mydrive/CDMA/config/geth_cdma_3rd_with_4th_sha256_sd${seed}.json
done

for seed in $(seq 44 48); do
    echo '----'
    echo '{
        "root path": "/media/gtnetuser/Dell/CDMA",
        "seed": #sd,
        "A path": "/media/gtnetuser/HDD_16TB_ATLAS/geth-3rd-state/state-22392874-sha256-sorted.txt",
        "B path": "/media/gtnetuser/HDD_16TB_ATLAS/geth-4th-state/state-22399992-sha256-sorted.txt",
        "Intersect path": "/media/gtnetuser/HDD_16TB_ATLAS/state-22392874-intersect-22399992.txt",
        "result filename": "/media/mydrive/CDMA/result/geth_iblt_3rd_with_4th_sd#sd.json"
    }' \
    | sed "s/#sd/${seed}/g"\
    > /media/mydrive/CDMA/config/geth_iblt_3rd_with_4th_sha256_sd${seed}.json
    build/bin/two_party_iblt /media/mydrive/CDMA/config/geth_iblt_3rd_with_4th_sha256_sd${seed}.json
done
