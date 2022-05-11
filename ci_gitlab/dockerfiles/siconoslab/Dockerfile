FROM IMAGENAME
USER root
RUN cd /builds/nonsmooth/siconos/build; make install
#COPY $CI_PROJECT_DIR/build/requirements.txt /home
RUN python3 -m pip install -U pyhull lxml matplotlib vtk
USER $NB_USER
ENV SICONOS_INSTALL_DIR=/home/install-siconos PATH=/home/install-siconos/bin/:$PATH
ENV JUPYTER_ENABLE_LAB yes
RUN cd $HOME; git clone https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos-tutorials.git


