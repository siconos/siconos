FROM intel/oneapi-hpckit:latest
RUN apt update  && apt upgrade -y &&  apt install -y -qq \
libgl1 
ENV CXX icpx
ENV CC icx
ENV FC ifx
# ENV PATH /opt/intel/oneapi/intelpython/latest/bin/:$PATH
ENV CI_PROJECT_DIR /home
RUN mkdir -p $CI_PROJECT_DIR/build
RUN mkdir -p $CI_PROJECT_DIR/deps

# -- Install micromamba --
COPY ci_gitlab/dockerfiles/install_intel_micromambaenv.sh  /home/
RUN (curl -L  micro.mamba.pm/install.sh) | bash
COPY ci_gitlab/dockerfiles/siconosintelenv.yml  /home/
RUN  grep -v '[ -z "\$PS1" ] && return' /root/.bashrc  > /root/newbashrc  && mv  /root/newbashrc /root/.bashrc

# -- Install siconosenv env --
RUN chmod +x /home/install_intel_micromambaenv.sh  && /home/install_intel_micromambaenv.sh 
