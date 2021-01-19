#!/bin/bash

apt-get update
apt-get install npm
apt-get install make
apt-get install g++
apt-get install build-essential
apt-get install sqlite3
apt-get install libsqlite3-dev
apt-get install zlib1g-dev
apt-get install unzip

git clone https://github.com/mapbox/tippecanoe.git
cd tippecanoe && make -j && make install

npm install -g mapshaper

cd ..
conda create --file conda_env.yml
conda activate espana_en_cifras
