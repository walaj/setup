## on HMS server
module load omero
## or otherwise have to install: https://omero.readthedocs.io/en/stable/users/cli/installation.html

# on an interactive node, login to Omero - using Harvard Key password
omero login -u $USER -s omero-app.hms.harvard.edu -t 86400

#
omero sessions group 'CycIF Sharing'

omero download Image:1614258 .

## Then run that command for each image ID you need. If an image comprises multiple files, they'll all get downloaded. For cycif ome-tiffs it will always be one file, but some H&E scanner formats like .vsi will come as a set of multiple files. In that case you might want to download into a specific directory by changing that final . argument to some other path. you might want to run the downloads as slurm jobs so you aren't tied to the ssh session while they run

sbatch --mem 4g -t 0-8 -p short --wrap 'omero download Image:1614258 .'

# (be sure to put the command to run in quotes when you pass it as the --wrap argument like this)
