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



### f you need to download from other groups, change the omero sessions group command. omero sessions are tied to a single group at a time. You just give it the string name of the group (case sensitive)

### systematically get all image ids

I can give you 2 options. One hack is to shift-click select all the images in the left-hand menu, then click the "link" button in the upper right. the URL it gives you will include all the IDs but you need to parse the string a bit.


the other way is omero cli

omero hql -q --style csv 'select img.id from Dataset ds join ds.imageLinks dil join dil.child img where ds.id=17812'



### all together with sbatch
omero hql -q --style csv \
  'select img.id from Dataset ds join ds.imageLinks dil join dil.child img where ds.id=17812' \
| tail -n +2 | cut -d, -f2 \
| while read -r id; do
    sbatch --job-name omero_${id} --mem 4g -t 0-8 -p short \
      --wrap "omero download Image:${id} ."
  done


### all together without sbatcy
omero hql -q --style csv \
  'select img.id from Dataset ds join ds.imageLinks dil join dil.child img where ds.id=17812' \
| tail -n +2 | cut -d, -f2 \
| while read -r id; do
    omero download Image:${id} .
  done


### all together, view only
omero hql -q --style csv \
  'select img.id from Dataset ds join ds.imageLinks dil join dil.child img where ds.id=17812' \
| tail -n +2 \
| cut -d, -f2 \
| while read -r id; do
    echo "omero download Image:${id} ."
  done
