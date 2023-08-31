#!/bin/bash
source ~/.bashrc
micromamba config prepend channels conda-forge
micromamba config set channel_priority strict
micromamba self-update
micromamba create -f /home/siconosenv.yml
