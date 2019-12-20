#!/bin/bash

echo "This script needs wget"
echo "Downloading databases from http://klif.uu.nl/download/metagenomics_db/"
echo "Please be patient. This will take some time (but only once)"
echo
echo "Downloading Kraken2 db"
wget http://klif.uu.nl/download/metagenomics_db/Kraken2_Nanodb.tar.gz
echo "Downloading host filter db (Mash db)"
wget http://klif.uu.nl/download/metagenomics_db/mash_db.tar.gz
echo "Downloading KMA db"
wget http://klif.uu.nl/download/metagenomics_db/KMA_ResFinder.tar.gz

echo "Untarring files"
tar xzvf Kraken2_Nanodb.tar.gz
tar xzvf mash_db.tar.gz
tar xzvf KMA_ResFinder.tar.gz

echo "Done"
echo ""
echo "please do 'rm *.tar.gz' to remove temporary download files"
exit 1