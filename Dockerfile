FROM ubuntu

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update && apt-get install -y \
    libeigen3-dev \
    cmake \
    build-essential \
    git \
    cmake \
    libx11-dev \
    mesa-common-dev libgl1-mesa-dev libglu1-mesa-dev \
    libxrandr-dev \
    libxi-dev \
    libxmu-dev \
    libblas-dev \
    libxinerama-dev \
    libxcursor-dev \
    python2


## DOWNLOAD VTK SOURCE
RUN mkdir -p /projects
WORKDIR /projects
RUN git clone --progress --verbose http://github.com/Kitware/VTK.git
WORKDIR VTK

## UPDATE THE CODE
# RUN git fetch origin
# RUN git rebase origin/master

## CONFIGURE (similar to the use of ccmake) VTK DEBUG BUILD
ADD setup.py /projects/VTK/
RUN chmod +x setup.py
RUN mkdir -p /projects/build
WORKDIR /projects/build
RUN echo "Following advice from http://www.vtk.org/Wiki/VTK/Configure_and_Build#Configure_VTK_with_CMake"
RUN python2 /projects/VTK/setup.py --type=debug --timings -DBUILD_SHARED_LIBS:BOOL=ON -DBUILD_EXAMPLES:BOOL=ON -DBUILD_TESTING:BOOL=ON /projects/build

# ## COMPILE VTK
WORKDIR /projects/build
RUN make -j4
RUN make install

COPY code /code/
WORKDIR /code
RUN sh build.sh