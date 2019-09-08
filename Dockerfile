FROM gw000/keras:2.1.4-py3-tf-cpu

RUN apt-get update
RUN apt-get install --no-install-recommends -y python-matplotlib

RUN apt-get install -y dirmngr --install-recommends
RUN apt-get install -y apt-transport-https software-properties-common gnupg2

RUN apt-key adv --keyserver keys.gnupg.net --recv-key 'E19F5F87128899B192B1A2C2AD5F960A256A04AF'
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/debian stretch-cran35/'
RUN apt-get update && apt-get install -y r-base

# Install python packages
RUN pip3 install pandas h5py seaborn scipy scikit-learn scanpy

# Install ggplot2
RUN apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev libcairo2-dev libgeos-dev libudunits2-dev libgdal-dev
RUN Rscript -e "install.packages(c('scales','ggplot2'),dependencies=TRUE, repos='http://cran.rstudio.com/')"

# Install DeepImpute
RUN git clone "https://github.com/lanagarmire/deepimpute.git" && cd deepimpute && pip3 install .

RUN pip3 install ipython
RUN apt-get install nano

RUN mkdir -p results/accuracy results/downstream results/fish results/dropout_effect results/speed_memory results/training_w_subsets

COPY . .

CMD /bin/bash
