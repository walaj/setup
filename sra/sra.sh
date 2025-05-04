## install SRA
cd /tmp
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -xvzf sratoolkit.current-ubuntu64.tar.gz
sudo mv sratoolkit.*-ubuntu64 /opt/sratoolkit
echo 'export PATH=/opt/sratoolkit/bin:$PATH' >> ~/.bashrc
cd

## Create a list of SRR accession numbers
# - go to the SRA Run Selector for PRJN*** and download accession list to get all SRR numbers

## fetch sra files and extract to fastq
while read srr; do
  echo "Processing $srr"
  prefetch $srr --output-directory ./sra_files
  fasterq-dump ./sra_files/$srr --outdir ./fastq_files --split-files
done < SRR_Acc_List
