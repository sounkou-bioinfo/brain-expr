# Auto-hébergement de l'app :  docker build -t brain-expr . && docker run -p 3838:3838 brain-expr
# puis http://localhost:3838
FROM rocker/shiny:4.4.1
RUN install2.r --error --skipinstalled data.table ggplot2 DT bslib && rm -rf /tmp/downloaded_packages
RUN rm -rf /srv/shiny-server/*
COPY app/ /srv/shiny-server/
RUN chown -R shiny:shiny /srv/shiny-server
EXPOSE 3838
CMD ["/usr/bin/shiny-server"]
