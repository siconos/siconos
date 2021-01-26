FROM ubuntu:20.04
ENV TZ=Europe/Paris
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone && apt update  && apt install -y -qq \
        ntp \
        cmake \
        git-core \
        make \
        libboost-dev \
        libboost-filesystem-dev \
        libboost-timer-dev \
        libboost-chrono-dev \
        libgmp-dev \
        gcc \
        gfortran \
        libgfortran-10-dev \
        g++ \
        libopenblas-dev \
        liblapacke-dev \
        libcppunit-dev \
        vim && apt autoclean -y && apt autoremove -y && rm -rf /var/lib/apt/lists/*
