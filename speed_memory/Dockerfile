FROM r-base

# Install python
RUN apt-get update && apt-get install -y python3-pip python3-dev python3-tk
RUN cd /usr/local/bin \
	&& ln -s /usr/bin/python3 python \
	&& pip3 install --upgrade pip

RUN apt-get install -y nano vim emacs git 
RUN apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev libcairo2-dev libgeos-dev libudunits2-dev libgdal-dev procps default-jre

# Install ggplot2
RUN Rscript -e "install.packages(c('scales','ggplot2','devtools','Rcpp','RcppArmadillo','quadprog','glmnet'),dependencies=TRUE, repos='http://cran.rstudio.com/')"

# Install python packages
RUN pip3 install ipython pandas==0.23 h5py matplotlib seaborn scipy scikit-learn numpy==1.16.4
RUN pip3 install scanpy dca magic-impute

# Install DeepImpute
RUN git clone "https://github.com/lanagarmire/deepimpute.git" && cd deepimpute && pip3 install .

# scImpute
RUN Rscript -e "devtools::install_github('Vivianstats/scImpute@v0.0.8')"

# SAVER
RUN Rscript -e "devtools::install_github('mohuangx/SAVER')"

# DrImpute
RUN Rscript -e "install.packages('DrImpute')"

# VIPER
RUN Rscript -e "devtools::install_github('ChenMengjie/VIPER')"

# Install nextflow
RUN wget -qO- https://get.nextflow.io | bash
RUN mv nextflow /usr/local/bin

WORKDIR /workspace
RUN mkdir -p results/speed_memory paper_data/speed_memory speed_memory

COPY . speed_memory/

CMD /bin/bash
