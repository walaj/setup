sudo apt update && sudo apt upgrade -y
sudo apt-get install -y emacs git tree bcftools bzip2 tmux \
     samtools libcurl4-openssl-dev bwa cmake autoconf automake cloud-guest-utils \
     libtool libboost-all-dev htop iotop dstat g++ libbz2-dev \
     libxml2-dev libssl-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
     libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
     libhts-dev libopenslide0 libopenslide-dev 
export BOOST_ROOT=/usr/include/boost

# Reference genome
mkdir -p ~/ref/
cd ~/ref
gsutil cp gs://osteosarc-genomics/ref_genome/* .
cd
## FYI these came from below, but are now in the above
#gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta .
#gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai .
#gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict .
## NB: these are identical to the GATK.GRCh38 human genome from sarek:
## aws s3 cp s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta     . --no-sign-request


## get customizations from github
mkdir git
cd git
git clone https://github.com/walaj/setup
cat ~/git/setup/bash/.bashrc >> ~/.bashrc
git config --global user.name "Jeremiah Wala"
git config --global user.email "jeremiah.wala@gmail.com"
cd

## setup ssh key
ssh-keygen -t rsa -b 4096 -C "jeremiah.wala@gmail.com"
## add key to ssh agent
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_rsa
## copy to clipboard
cat ~/.ssh/id_rsa.pub
##Add the key to GitHub
#Go to: https://github.com/settings/keys
#Click New SSH key
#Give it a title like "Linux VM" or "Ubuntu dev box"
#Paste the key you just copied

## setup emacs and tmux customizations
mkdir ~/.emacs.d
cp ~/git/setup/emacs/init.el ~/.emacs.d/
cp ~/git/setup/emacs/lite.el ~/.emacs.d/
cp ~/git/setup/tmux/.tmux.conf ~/

############
# GCloud
############

## login to google cloud
gcloud auth login

## attach disk
gcloud compute instances attach-disk VM_NAME --disk=DISK_NAME --zone=ZONE

## format if needed (if first time used)
#sudo mkfs.ext4 /dev/sdb

## create mount point and make persistent on rebood
sudo mkdir -p /mnt/wgs
sudo mount /dev/sdb /mnt/wgs
## have to do this manually
sudo blkid /dev/sdb ## make note of the UUID
## with sudo emacs -nw /etc/fstab add this: UUID=1234-ABCD /mnt/wgs ext4 defaults,nofail 0 2
## test: sudo mount -a -- if no errors then OK  or test with df -h

## give user permsissions on the disk
sudo chown -R $USER:$USER /mnt/wgs

############
# Bioinformatics tools
############

## install SRA toolkit
cd /tmp
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -xvzf sratoolkit.current-ubuntu64.tar.gz
sudo mv sratoolkit.*-ubuntu64 /opt/sratoolkit
echo 'export PATH=/opt/sratoolkit/bin:$PATH' >> ~/.bashrc
cd

## install VEP ensembl
sudo apt-get install -y \
    git unzip curl build-essential \
    zlib1g-dev libperl-dev \
    libdbi-perl libdbd-mysql-perl \
    libjson-perl liblzma-dev \
    cpanminus tabix

## install MaSurca
#wget tarball from github
# in tar dir, BOOST_ROOT=/usr ./install.sh

###########
# PYTHON
###########

## install conda and mamba
curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-*.sh    ## select "yes"
source ~/.bashrc
conda install -n base -c conda-forge mamba

############
# NF CORE
############

## install nextflow
sudo apt install openjdk-17-jdk
curl -s https://get.nextflow.io | bash
mkdir ~/software/
mv nextflow ~/software/
echo 'export PATH=$HOME/software:$PATH' >> ~/.bashrc
source ~/.bashrc

## install Docker
sudo apt install apt-transport-https ca-certificates curl software-properties-common  ## dependencies
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -   ## add key
sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/debian $(lsb_release -cs) stable"
sudo apt install docker-ce
docker --version ## verify

## start docker damen
sudo systemctl start docker
sudo systemctl enable docker ## auto start on reboot
sudo systemctl status docker ## verify
sudo usermod -aG docker $USER ## run without sudo
newgrp docker ## apply changes

## pull nf-core module
nextflow pull nf-core/sarek


##########
## SHINY and R
##########

# Install R and tools
sudo apt-get install -y r-base gdebi-core


# Install Shiny package in R
sudo su - -c "R -e \"install.packages('shiny', repos='https://cloud.r-project.org/')\""
sudo su - -c "R -e \"install.packages('data.table', repos='https://cloud.r-project.org/')\""
sudo su - -c "R -e \"install.packages('tidyverse', repos='https://cloud.r-project.org/')\""
sudo su - -c "R -e \"install.packages('plotly', repos='https://cloud.r-project.org/')\""
sudo su - -c "R -e \"install.packages('devtools', repos='https://cloud.r-project.org/')\""
sudo su - -c "R -e \"install.packages('BiocManager', repos='https://cloud.r-project.org/')\""
sudo su - -c "R -e \"BiocManager::install('GenomicRanges')\""
sudo su - -c "R -e \"BiocManager::install('rtracklayer')\""

# Install Shiny Server
wget https://download3.rstudio.org/ubuntu-20.04/x86_64/shiny-server-1.5.23.1030-amd64.deb
sudo gdebi -n shiny-server-1.5.23.1030-amd64.deb

# Start Shiny Server
sudo systemctl enable --now shiny-server

# Download Sarcoma Shiny app and data
mkdir ~/data
cd ~/data && gsutil cp gs://osteosarc-genomics/web/shiny_data/* .
cd ~/git && git clone https://github.com/walaj/sarcoma-shiny

# Setup shiny app
sudo mkdir -p /srv/shiny-server/sarcoma-app
