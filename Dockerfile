FROM bioconductor/bioconductor_docker:latest
MAINTAINER Nuno Agostinho <nunodanielagostinho@gmail.com>

RUN apt-get update && apt-get -y upgrade && apt-get -y autoremove

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('remotes')"

# Copy description
WORKDIR psichomics
ADD . .

# Install R package from source
RUN Rscript -e "remotes::install_local()"
RUN rm -rf *

# # To start an R session with this version installed:
# docker run -ti [docker image] R
# library(R)
