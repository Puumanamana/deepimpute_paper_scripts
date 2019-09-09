FROM r-base

RUN apt-get update

# Install python
RUN apt-get install -y python3-pip python3-dev python3-tk
RUN cd /usr/local/bin \
	&& ln -s /usr/bin/python3 python \
	&& pip3 install --upgrade pip

RUN apt-get install -y nano vim emacs git 
RUN apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev libcairo2-dev libgeos-dev libudunits2-dev libgdal-dev

# Install ggplot2
RUN Rscript -e "install.packages(c('scales','ggplot2'),dependencies=TRUE, repos='http://cran.rstudio.com/')"

# Install python packages
RUN pip3 install ipython pandas h5py matplotlib seaborn scipy scikit-learn numpy==1.16.4
RUN pip3 install scanpy leidenalg

# Install DeepImpute
RUN git clone "https://github.com/lanagarmire/deepimpute.git" && cd deepimpute && pip3 install .

WORKDIR /workspace

RUN mkdir results && cd results && mkdir accuracy downstream fish dropout_effect speed_memory training_w_subsets && cd ..

COPY . .

CMD /bin/bash
