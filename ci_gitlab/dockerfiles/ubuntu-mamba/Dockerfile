FROM ubuntu:22.04
ENV TZ=Europe/Paris
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone && apt update  && apt upgrade -y && apt install -y -qq \
        git-core \
        make \
	g++ \
        gfortran \
	# lp-solve \
        # liblpsolve55-dev \
	valgrind \
	vim \
	bzip2 \
	curl \
	wget
RUN apt autoclean -y && apt autoremove -y&& rm -rf /var/lib/apt/lists/*
WORKDIR /home
ENV CI_PROJECT_DIR /home
# -- Install micromamba --
COPY ci_gitlab/dockerfiles/install_micromambaenv.sh  /home/
RUN (curl -L  micro.mamba.pm/install.sh) | bash
COPY ci_gitlab/dockerfiles/siconosenv.yml  /home/
RUN  grep -v '[ -z "\$PS1" ] && return' /root/.bashrc  > /root/newbashrc  && mv  /root/newbashrc /root/.bashrc

# -- Install siconosenv env --
RUN chmod +x /home/install_micromambaenv.sh  && /home/install_micromambaenv.sh 
