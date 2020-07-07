FROM php:5.6-fpm
MAINTAINER yhy <yuanhy@bcc.ac.cn>

######
# You can install php extensions using docker-php-ext-install
######
RUN echo "deb http://mirrors.163.com/debian/ stretch main non-free contrib"
RUN echo "deb http://mirrors.163.com/debian/ jessie-updates main non-free contrib"
RUN echo "deb http://mirrors.163.com/debian/ jessie-backports main non-free contrib"
RUN echo "deb-src http://mirrors.163.com/debian/ jessie main non-free contrib"
RUN echo "deb-src http://mirrors.163.com/debian/ jessie-updates main non-free contrib"
RUN echo "deb-src http://mirrors.163.com/debian/ jessie-backports main non-free contrib"
RUN echo "deb http://mirrors.163.com/debian-security/ jessie/updates main non-free contrib"
RUN echo "deb-src http://mirrors.163.com/debian-security/ jessie/updates main non-free contrib"
RUN apt-get update && apt-get install -y \
        libfreetype6-dev \
        libjpeg62-turbo-dev \
        libmcrypt-dev \
        libpng-dev \
    && docker-php-ext-install pdo_mysql \
    && docker-php-ext-install -j$(nproc) iconv mcrypt \
    && docker-php-ext-configure gd --with-freetype-dir=/usr/include/ --with-jpeg-dir=/usr/include/ \
    && docker-php-ext-install -j$(nproc) gd

RUN docker-php-ext-install pdo_mysql \
 && apt-get install -y r-base 
RUN Rscript -e "install.packages('ggplot2', repos='https://mirrors.ustc.edu.cn/CRAN/')" 
RUN Rscript -e "source('http://bioconductor.org/biocLite.R');options(BioC_mirror='https://mirrors.tuna.tsinghua.edu.cn/bioconductor');biocLite();biocLite('ggtree')" 



