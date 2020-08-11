FROM bioconductor/bioconductor_docker:latest
MAINTAINER Nuno Agostinho <nunodanielagostinho@gmail.com>
    
RUN apt-get update && apt-get -y upgrade && apt-get -y autoremove

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile

# Install required packages to make it quicker to test new changes
RUN Rscript -e "install.packages('remotes')"
RUN Rscript -e "install.packages('shiny')"
RUN Rscript -e "install.packages('shinyBS')"
RUN Rscript -e "install.packages('DT')"
RUN Rscript -e "install.packages('data.table')"
RUN Rscript -e "install.packages('reshape2')"
RUN Rscript -e "install.packages('dplyr')"
RUN Rscript -e "install.packages('XML')"
RUN Rscript -e "install.packages('highcharter')"
RUN Rscript -e "install.packages('Rcpp')"

RUN Rscript -e "BiocManager::install('limma')"
RUN Rscript -e "BiocManager::install('AnnotationHub')"
RUN Rscript -e "BiocManager::install('recount')"
RUN Rscript -e "BiocManager::install('edgeR')"

RUN Rscript -e "install.packages('miniUI')"
RUN Rscript -e "install.packages('R.oo')"
RUN Rscript -e "install.packages('fastICA')"
RUN Rscript -e "install.packages('fastmatch')"
RUN Rscript -e "install.packages('ggrepel')"
RUN Rscript -e "install.packages('shinyjs')"
RUN Rscript -e "install.packages('R.utils')"
RUN Rscript -e "install.packages('xfun')"
RUN Rscript -e "install.packages('ps')"
RUN Rscript -e "install.packages('glue')"
RUN Rscript -e "install.packages('backports')"
RUN Rscript -e "install.packages('pairsD3')"

RUN Rscript -e "BiocManager::install('org.Hs.eg.db')"
RUN Rscript -e "install.packages('colourpicker')"
RUN Rscript -e "install.packages('ellipsis')"

# Copy description
WORKDIR psichomics
ADD DESCRIPTION .
ADD CONDUCT.md .
ADD Dockerfile .
ADD LICENSE .
ADD man man
ADD NAMESPACE .
ADD NEWS.md .
ADD R R
ADD README.md .
ADD src src
ADD tests tests
ADD vignettes vignettes

# Install R package from source
RUN Rscript -e "remotes::install_local()"

# # To start an R session with this version installed:
# docker run -ti [docker image] R
# library(R)
