#!/bin/bash

#seq_name='V033mat'
#sim_dir='/export/home/WG-shahlo/projects/HIV/simulations/published_antibodies/V033mat/'
#cdrh3_sequence='ARVDGDDYGYFDTVPGDSKKYYFKH' # V033

#seq_name='CH04'
#sim_dir='/export/home/WG-shahlo/projects/HIV/simulations/published_antibodies/CH04/'
#cdrh3_sequence='ARGTDYTIDDQGIRYQGSGTFWYFDV' # CH04

#seq_name='seq22299'
#sim_dir='/export/home/WG-shahlo/projects/HIV/simulations/antibodies_from_OAS/rhesus_like_axe/search_3/seq22299/'
#cdrh3_sequence='GKDRTGDDYGFWGGQTKQPYYFDY' # seq22299

#seq_name='seq10941'
#sim_dir='/export/home/WG-shahlo/projects/HIV/simulations/antibodies_from_OAS/ch01_like_axe/search_7/seq10941/'
#cdrh3_sequence='ARVRTPSRRYYYGSGSYLGYYFDY' # seq10941 

#seq_name='seq10821'
#sim_dir='/export/home/WG-shahlo/projects/HIV/simulations/antibodies_from_OAS/ch01_like_axe/search_7/seq10821/'
#cdrh3_sequence='ARSSATHGRYYYGSGSYYKRGFDP' # seq10821

#seq_name='seq215'
#sim_dir='/export/home/WG-shahlo/projects/HIV/simulations/antibodies_from_OAS/rhesus_like_axe/search_3/seq215/'
#cdrh3_sequence='ARDGDEDDYGGTSDIDS' # seq215

#seq_name='seq17216'
#sim_dir='/export/home/WG-shahlo/projects/HIV/simulations/antibodies_from_OAS/rhesus_like_axe/search_3/seq17216/'
#cdrh3_sequence='ARGDDYGDYDRLGGEDIPKNWFDP' # seq17216

#seq_name='seq15061'
#sim_dir='/export/home/WG-shahlo/projects/HIV/simulations/antibodies_from_OAS/ch01_like_axe/search_7/seq15061/'
#cdrh3_sequence='ARVLNKRGLIVSRNRTYYYGSGSYSMQFDP' # seq15061

#seq_name='seq20621'
#sim_dir='/export/home/WG-shahlo/projects/HIV/simulations/antibodies_from_OAS/ch01_like_axe/search_7/seq20621/'
#cdrh3_sequence='ARAPLPEYCSAGTCYQGSGYFDL' # seq20521

seq_name='seq12552'
sim_dir='/export/home/WG-shahlo/projects/HIV/simulations/antibodies_from_OAS/rhesus_like_axe/search_3/seq12552/'
cdrh3_sequence='AKSGRGDDYGYVWGSYRYNRGDFDY' # seq12552

#seq_name='Ab41328'
#sim_dir='/export/home/WG-shahlo/projects/HIV/simulations/published_antibodies/Ab41328/'
#cdrh3_sequence='ARGSIYYEDDDGYYYSEATYLHLHLW' # Ab41328

for i in {1..10}; do
	which_run='run'$i;
	echo $which_run;
	sbatch slurm_submission.sh $seq_name $sim_dir $cdrh3_sequence $which_run;
done;

