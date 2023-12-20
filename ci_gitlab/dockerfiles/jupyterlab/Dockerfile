##  - Prepare an image including (and starting) Jupyter lab
ARG REGISTRY=gricad-registry.univ-grenoble-alpes.fr
ARG PROJECT=nonsmooth/siconos
ARG DEFAULT_IMAGE=sources/ubuntu22.04
ARG IMAGENAME=$REGISTRY/$PROJECT/$DEFAULT_IMAGE
# - Prepare a jupyter notebook, based on a Docker image with a complete Siconos install. Add (clone) the siconos-tutorials repository
FROM $IMAGENAME AS jupylab
RUN pip install --no-cache-dir -U jupyterlab ipykernel
ENV WORK=/home/siconos-user
WORKDIR $WORK
SHELL ["/bin/bash", "-o", "pipefail", "-c"]
#ENV SHELL /bin/bash
ENV JUPYTER_PORT=8888
ENV JUPYTER_ENABLE_LAB yes
EXPOSE $JUPYTER_PORT
ENTRYPOINT ["siconosenv/bin/jupyter", "lab","--ip=0.0.0.0","--allow-root", "--no-browser"]


