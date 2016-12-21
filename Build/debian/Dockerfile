FROM debian:sid

USER root

RUN apt-get update -y && apt-get install -y packaging-dev

# Modify to taste
RUN apt-get install -y emacs

# See README to generate this file
ADD gpg-secrets /root/gpg-secrets

WORKDIR /root
CMD /bin/bash
