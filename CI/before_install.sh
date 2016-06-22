#!/bin/sh

if test "$TRAVIS" = true; then 
    pip install --user PyYAML
    sudo apt-get -o Acquire::CompressionTypes::Order::=bz2 update
    sudo apt-get -y install software-properties-common
    sudo apt-get update -qq
    sudo apt-get -y install make
    sudo free -m
fi
