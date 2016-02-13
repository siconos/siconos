#!/bin/sh

if test "$TRAVIS" = true; then 
    pip install --user PyYAML
    sudo apt-get -o Acquire::CompressionTypes::Order::=bz2 update
    sudo apt-get -y install software-properties-common
#    sudo add-apt-repository ppa:smspillaz/cmake-2.8.12 -y
    sudo apt-get update -qq
#    sudo apt-get purge cmake -qq
    sudo apt-get -y install make
    sudo free -m
fi
