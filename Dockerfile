FROM ubuntu:18.04
RUN apt-get update
RUN apt-get install -y cmake \         
	&& apt-get install    -y    flex bison  \
	&& apt-get install    -y    gfortran    \            
	&& apt-get install    -y    git         \           
	&& apt-get install    -y    g++                   \
	&& apt-get install    -y    libboost-dev \        
	&& apt-get install    -y     libboost-python-dev \
	&& apt-get install    -y     libboost-filesystem-dev \
	&& apt-get install    -y     libboost-program-options-dev \
	&& apt-get install    -y     libboost-signals-dev   \   
	&& apt-get install    -y     libboost-system-dev    \  
	&& apt-get install    -y    libboost-thread-dev     \ 
	&& apt-get install    -y    libcfitsio-dev     \     
	&& apt-get install    -y    libffi-dev       \      
	&& apt-get install    -y    libfftw3-dev    \      
	&& apt-get install    -y    libgsl-dev     \        
	&& apt-get install    -y  liblog4cxx-dev   \         
	&& apt-get install    -y  libopenblas-dev   \        
	&& apt-get install    -y  libopenmpi-dev     \       
	&& apt-get install    -y  libpython-dev    \ 
	&& apt-get install    -y  make              \ 
	&& apt-get install    -y  patch             \   
	&& apt-get install    -y  python-pip       \      
	&& apt-get install    -y  subversion       \        
	&& apt-get install    -y  wget             \        
        && apt-get install    -y  docker           \
	&& apt-get install    -y  wcslib-dev                
RUN cd /usr/local 
RUN wget ftp://ftp.astron.nl/outgoing/Measures/WSRT_Measures.ztar
RUN mv WSRT_Measures.ztar WSRT_Measures.tar.gz
RUN gunzip WSRT_Measures.tar.gz
RUN tar -xvf WSRT_Measures.tar
RUN rm WSRT_Measures.tar
RUN mkdir /var/lib/jenkins
RUN mkdir /var/lib/jenkins/workspace



