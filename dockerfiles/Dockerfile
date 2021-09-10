FROM gcc:7.5

ARG LIBDIR=/home/libraries
ARG SRCDIR=/home/src

WORKDIR $SRCDIR


# 3rd party dependencies
#########################################################################################################
# Install gfortran, blas, lapack
RUN apt-get update \
    &&  apt-get install -y libblas-dev liblapack-dev gfortran-7

# Install cfitsio
RUN wget http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-3.47.tar.gz \
    &&	tar -xvf cfitsio-3.47.tar.gz \
    &&	mv cfitsio-3.47 cfitsio \
    &&	cd cfitsio \
    &&	./configure --prefix=$LIBDIR/cfitsio \
    &&	make \
    &&	make install

# Install CCfits
RUN wget https://heasarc.gsfc.nasa.gov/fitsio/CCfits-2.5/CCfits-2.5.tar.gz \
    &&	tar -xvf CCfits-2.5.tar.gz \
    &&	cd CCfits \
    &&	./configure --prefix=$LIBDIR/CCfits --with-cfitsio=$LIBDIR/cfitsio \
    &&	make \
    &&	make install

# Install gmp
RUN wget https://gmplib.org/download/gmp/gmp-6.1.2.tar.xz \
    &&	tar -xvf gmp-6.1.2.tar.xz \
    &&	cd gmp-6.1.2 \
    &&	./configure --prefix=$LIBDIR/gmp \
    &&	make \
    &&	make install

# Install mpfr
RUN wget http://www.mpfr.org/mpfr-4.0.2/mpfr-4.0.2.tar.gz \
    &&	tar -xvf mpfr-4.0.2.tar.gz \
    &&	cd mpfr-4.0.2 \
    &&	./configure --with-gmp=$LIBDIR/gmp --prefix=$LIBDIR/mpfr \
    &&	make \
    &&	make install

# Install boost
RUN wget https://boostorg.jfrog.io/artifactory/main/release/1.66.0/source/boost_1_66_0.tar.gz \
    &&	tar -xvf boost_1_66_0.tar.gz \
    &&	cd boost_1_66_0 \
    &&	./bootstrap.sh --with-libraries=thread \
    &&	./b2 --prefix=$LIBDIR/boost install

# Install cmake
RUN wget https://github.com/Kitware/CMake/releases/download/v3.15.5/cmake-3.15.5-Linux-x86_64.sh \
    -q -O /tmp/cmake-install.sh \
    &&	chmod u+x /tmp/cmake-install.sh \
    &&	mkdir /usr/bin/cmake \
    &&	/tmp/cmake-install.sh --skip-license --prefix=/usr/bin/cmake \
    &&	rm /tmp/cmake-install.sh
ENV PATH=$PATH:/usr/bin/cmake/bin

# Install CGAL
RUN wget https://github.com/CGAL/cgal/archive/releases/CGAL-4.11.2.tar.gz \
    &&	tar -xvf CGAL-4.11.2.tar.gz \
    &&	cd cgal-releases-CGAL-4.11.2 \
    &&	mkdir build_CGAL \
    &&	cd build_CGAL \
    &&	cmake \
    -DGMP_INCLUDE_DIR=$LIBDIR/gmp/include \
    -DGMP_LIBRARIES=$LIBDIR/gmp/lib \
    -DMPFR_INCLUDE_DIR=$LIBDIR/mpfr/include \
    -DMPFR_LIBRARIES=$LIBDIR/mpfr/lib \
    -DBOOST_INCLUDEDIR=$LIBDIR/boost/include \
    -DCMAKE_INSTALL_PREFIX=$LIBDIR/CGAL \
    .. \
    &&	make \
    &&	make install

# Install GSL
RUN wget ftp://ftp.gnu.org/gnu/gsl/gsl-2.7.tar.gz \
    &&  tar -xvf gsl-2.7.tar.gz \
    &&  cd gsl-2.7 \
    &&  ./configure --prefix=$LIBDIR/gsl \
    &&  make \
    && make install

# Install jsoncpp
RUN git clone https://github.com/open-source-parsers/jsoncpp.git \
    &&	cd jsoncpp \
    &&	mkdir -p build \
    &&	cd build \
    &&	cmake \
    -DCMAKE_BUILD_TYPE=release \
    -DBUILD_SHARED_LIBS=ON \
    -DARCHIVE_INSTALL_DIR=. \
    -DCMAKE_INSTALL_PREFIX=$LIBDIR/jsoncpp \
    -G "Unix Makefiles" \
    .. \
    &&	make \
    &&	make install

# Install openmpi
RUN wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.1.tar.gz \
    &&  tar -xvf openmpi-4.1.1.tar.gz \
    &&  cd openmpi-4.1.1 \
    &&  ./configure --prefix=$LIBDIR/openmpi \
    &&  make all install

# Install Eigen 3
RUN wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.bz2 \
    &&  tar -xvf eigen-3.4.0.tar.bz2 \
    &&  cd eigen-3.4.0 \
    &&  mkdir build_dir \
    &&  cd build_dir \
    &&  cmake -DCMAKE_INSTALL_PREFIX=$LIBDIR/eigen .. \
    &&  make install

# Install Multinest
RUN git clone https://github.com/JohannesBuchner/MultiNest.git \
    &&  cd MultiNest/build \
    &&  cmake -DCMAKE_PREFIX_PATH=$LIBDIR/openmpi -DCMAKE_INSTALL_PREFIX=$LIBDIR/multinest .. \
    &&  make \
    &&  make install


### Uncomment the line below to use cache up to here and update the libraries and code with their latest versions
ARG DUMMY_DATE=2011-11-04

ENV PATH $PATH:$LIBDIR/openmpi/bin:$LIBDIR/MultiNest/bin
ENV CPATH $CPATH:$LIBDIR/cfitsio/include:$LIBDIR/CCfits/include:$LIBDIR/gmp/include:$LIBDIR/mpfr/include:$LIBDIR/boost/include:$LIBDIR/CGAL/include:$LIBDIR/jsoncpp/include:$LIBDIR/gsl/include:$LIBDIR/eigen/include/eigen3:$LIBDIR/openmpi/include:$LIBDIR/multinest/include
ENV LIBRARY_PATH $LIBRARY_PATH:$LIBDIR/cfitsio/lib:$LIBDIR/CCfits/lib:$LIBDIR/gmp/lib:$LIBDIR/mpfr/lib:$LIBDIR/boost/lib:$LIBDIR/CGAL/lib:$LIBDIR/jsoncpp/lib:$LIBDIR/gsl/lib:$LIBDIR/openmpi/lib:$LIBDIR/multinest/lib
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:$LIBDIR/cfitsio/lib:$LIBDIR/CCfits/lib:$LIBDIR/gmp/lib:$LIBDIR/mpfr/lib:$LIBDIR/boost/lib:$LIBDIR/CGAL/lib:$LIBDIR/jsoncpp/lib:$LIBDIR/gsl/lib:$LIBDIR/openmpi/lib:$LIBDIR/multinest/lib
#########################################################################################################





# Install 'gerlumphpp' and 'vkl_lib' libraries
#########################################################################################################
# Install vkl_lib
RUN cd /home/ \
    &&	git clone https://github.com/gvernard/verykool.git \
    &&	cd verykool \ 
    &&	make verykool
#########################################################################################################

RUN rm -r -f $SRCDIR
WORKDIR /home/verykool
