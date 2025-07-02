sudo apt-get install -y \                                                                                                                                        
     git unzip curl build-essential \                                                                                                                             
     zlib1g-dev libperl-dev \                                                                                                                                     
     libdbi-perl libdbd-mysql-perl \                                                                                                                              
     libjson-perl liblzma-dev \                                                                                                                                   
     cpanminus tabix

## install
conda create -n vep_env -c bioconda -c conda-forge ensembl-vep
conda activate vep_env
cpan install Compress::Raw::Zlib

## download VEP cache
vep_install -a cf -s homo_sapiens -y GRCh38 -c ~/.vep
## or
#mkdir -p ~/.vep
#cd ~/.vep
#wget https://ftp.ensembl.org/pub/release-113/variation/indexed_vep_cache/homo_sapiens_vep_113_GRCh38.tar.gz


## Filter Mutect2 for somatic
bcftools view -f PASS -i 'INFO/TLOD>6 && INFO/NLOD>2' input.vcf > somatic.vcf

## Run VEP
vep -i somatic.vcf -o somatic.annotated.vcf \
  --vcf \
  --cache \
  --dir_cache ~/.vep \
  --offline \
  --assembly GRCh38 \
  --everything \
  --fork 4


