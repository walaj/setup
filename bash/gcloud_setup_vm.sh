## basic installations
sudo apt-get update
sudo apt-get install emacs 
sudo apt-get install git
sudo apt-get install tree
sudo apt-get install tmux
sudo apt-get install samtools

## get customizations from github
mkdir git
cd git
git clone https://github.com/walaj/setup
cat ~/git/setup/bash/.bashrc >> ~/.bashrc

## setup emacs and tmux customizations
mkdir ~/.emacs.d
cp ~/git/setup/emacs/init.el ~/.emacs.d/
cp ~/git/setup/emacs/lite.el ~/.emacs.d/
cp ~/git/setup/tmux/.tmux.conf ~/

## install conda and mamba
curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-*.sh    ## select "yes"
source ~/.bashrc
conda install -n base -c conda-forge mamba

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
blkid /dev/sdb ## make note of the UUID
## with sudo emacs -nw /etc/fstab add this: UUID=1234-ABCD  /mnt/wgs  ext4  defaults  0 2
## test: sudo mount -a -- if no errors then OK  or test with df -h

## give user permsissions on the disk
sudo chown -R $USER:$USER /mnt/wgs
