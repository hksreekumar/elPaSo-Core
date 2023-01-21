# This Dockerfile creates a Docker image with tubslatex installed.
# The base image is Ubuntu 20.04 and TeX Live will be installed using
# the packages provided by Ubuntu.


# Installation instructions:
# 1. Save this file in an empty directory
# 2. Download the latest tubslatex_intaller from https://www.tu-braunschweig.de/latex
#    and save it in the same directory
# 3. Execute in this directory: `docker build -t tubslatex .`

# Usage instructions:
# Linux:
#   Execute `docker run --rm -v $(pwd):/latex tubslatex`
#   in the directory where your latex file is.
#   This corresponds to `latexmk -pdf`.
#   If you want to build a specific file you have to specify the exact command.
#   For example: `docker run --rm -v $(pwd):/latex tubslatex latexmk -pdf main.tex`
#   or `docker run --rm -v $(pwd):/latex tubslatex pdflatex main.tex`.
#
# Windows:
#   In general the same but you might get problems with `$(pwd)` depending on the Shell you are using.


# For problems, suggestions or anything else visit https://git.rz.tu-bs.de/tubslatex/tubslatex/-/issues
# or contact us: tubs-latex@tu-braunschweig.de


FROM ubuntu:20.04
LABEL "description"="tubslatex on Ubuntu 20.04 with TeX Live" \
      "copyright"="TU Braunschweig" \
      "contact"="tubs-latex@tu-braunschweig.de" \
      "website"="https://www.tu-braunschweig.de/latex" \
      "issues"="https://git.rz.tu-bs.de/tubslatex/tubslatex/-/issues" \
      "git"="https://git.rz.tu-bs.de/tubslatex/tubslatex" \
      "packages"="apt-utils texlive-base texlive-latex-extra \
                  texlive-fonts-extra texlive-latex-recommended \
                  texlive-fonts-recommended texlive-lang-german \
                  texlive-extra-utils cm-super texlive-bibtex-extra \
                  latexmk biber build-essential wget git make \
                  zip unzip curl"

# set to noninteractive to avoid warnings
ENV DEBIAN_FRONTEND noninteractive

# install packages
RUN apt-get update -q \
    && apt-get install -qy apt-utils \
    # texlive
    && apt-get install -qy \
        texlive-base texlive-latex-extra texlive-fonts-extra \
        texlive-latex-recommended texlive-fonts-recommended \
        texlive-lang-german texlive-extra-utils cm-super \
        texlive-bibtex-extra latexmk biber \
    # build
    && apt-get install -qy \
        build-essential wget git make zip unzip curl \
    && rm -rf /var/lib/apt/lists/*
  
# copy tubslatex installer
COPY tubslatex_installer_*.sh /tmp/

# install tubslatex
RUN sh /tmp/tubslatex_installer_*.sh --noninteractive --del-local-font --getnonfreefonts --install-required \
    && rm /tmp/tubslatex_installer_*.sh

# reset to default
ENV DEBIAN_FRONTEND newt
# set workdir to /latex
WORKDIR /latex
CMD ["latexmk", "-pdf"]