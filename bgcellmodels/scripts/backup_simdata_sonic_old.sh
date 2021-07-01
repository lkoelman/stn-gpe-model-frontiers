#!/bin/bash

# File synchronization using Google Drive command line utility "gdrive":
# https://github.com/gdrive-org/gdrive

# Setting up a sync is explained in this answer:
# https://github.com/gdrive-org/gdrive/issues/205#issuecomment-423759052
# $ mkdir test
# $ gdrive upload -r test
# $ gdrive list -m 1 # -> get ID
# $ gdrive sync upload --keep-local test <ID>

if [[ $(hostname) == login ]]; then
    sonic_datadir=$HOME/scratch
    gdrive_syncid="TODO"
elif [[ $(hostname) == sonic ]]; then
    sonic_datadir=$HOME/storage
    gdrive_syncid=1pn3tqxZeOIE4JVuvd2EqmUTazmwhVh1d
else
    echo "Not on SONIC cluster. Aborting."
    exit 0
fi


gdrive sync upload --keep-local $sonic_datadir $gdrive_syncid