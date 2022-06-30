FROM rocker/r-bspm:20.04 as base

RUN Rscript -e "install.packages('Plumber', Ncpus=2)"

COPY ./R /app/R
COPY ./DESCRIPTION /app
COPY ./NAMESPACE /app
WORKDIR /app

RUN Rscript -e "remotes::install_local('.', dependencies=T, Ncpus=2)"

WORKDIR /
COPY ./Api /app/api 
WORKDIR /app/api

EXPOSE 80

ENTRYPOINT ["Rscript", "app.R"]