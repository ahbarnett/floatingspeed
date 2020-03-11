FROM ubuntu:18.04

# Avoid warnings by switching to noninteractive
ENV DEBIAN_FRONTEND=noninteractive

# This Dockerfile adds a non-root user with sudo access. Use the "remoteUser"
# property in devcontainer.json to use it. On Linux, the container user's GID/UIDs
# will be updated to match your local UID/GID (when using the dockerFile property).
# See https://aka.ms/vscode-remote/containers/non-root-user for details.
ARG USERNAME=vscode
ARG USER_UID=1000
ARG USER_GID=$USER_UID

# Uncomment the following COPY line and the corresponding lines in the `RUN` command if you wish to
# include your requirements in the image itself. It is suggested that you only do this if your
# requirements rarely (if ever) change.
# COPY requirements.txt /tmp/pip-tmp/

# Configure apt and install packages
RUN apt-get update \
    && apt-get -y install --no-install-recommends apt-utils dialog 2>&1 \
    #
    # Verify git, process tools, lsb-release (common in install instructions for CLIs) installed
    && apt-get -y install git openssh-client less iproute2 procps lsb-release \
    #
    # Update Python environment based on requirements.txt
    # && pip --disable-pip-version-check --no-cache-dir install -r /tmp/pip-tmp/requirements.txt \
    # && rm -rf /tmp/pip-tmp \
    #
    # Create a non-root user to use if preferred - see https://aka.ms/vscode-remote/containers/non-root-user.
    && groupadd --gid $USER_GID $USERNAME \
    && useradd -s /bin/bash --uid $USER_UID --gid $USER_GID -m $USERNAME \
    # [Optional] Add sudo support for the non-root user
    && apt-get install -y sudo \
    && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME\
    && chmod 0440 /etc/sudoers.d/$USERNAME \
    #
    # Clean up
    && apt-get autoremove -y \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/*

# Let the installed vscode extensions persist between rebuilds
# See: https://code.visualstudio.com/docs/remote/containers-advanced#_avoiding-extension-reinstalls-on-container-rebuild
RUN mkdir -p /home/$USERNAME/.vscode-server/extensions && chown -R $USERNAME /home/$USERNAME/.vscode-server

RUN apt-get update && apt-get install -y nano

RUN apt-get install -y software-properties-common

#########################################
# Install gfortran-9
RUN add-apt-repository -y ppa:ubuntu-toolchain-r/test
RUN apt-get update && apt-get install -y gfortran-9 gfortran-7

#########################################
# Install g++
RUN apt-get update && apt-get install -y g++

#########################################
# Install Julia
RUN apt-get install -y wget
RUN wget -qO- https://julialang-s3.julialang.org/bin/linux/x64/1.2/julia-1.2.0-linux-x86_64.tar.gz | tar -xzv
ENV PATH="/julia-1.2.0/bin:${PATH}"
USER vscode
RUN julia -e 'import Pkg; Pkg.add("BenchmarkTools")'
USER root

#########################################
### Install Python                                                               
RUN apt-get update && apt-get -y install git wget build-essential
RUN apt-get install -y python3 python3-pip
### Make sure we have python3 and a working locale
RUN (rm /usr/bin/python || true) && ln -s python3 /usr/bin/python && (rm /usr/bin/pip || true) && ln -s pip3 /usr/bin/pip
ENV LANG='en_US.UTF-8' LANGUAGE='en_US:en' LC_ALL='en_US.UTF-8'
RUN apt-get install -y locales && locale-gen en_US.UTF-8
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y python3-tk
# numpy and numba
RUN pip install numpy numba

#########################################
# Install octave
RUN apt-get update && apt-get -y install octave

ENV DEBIAN_FRONTEND=dialog
