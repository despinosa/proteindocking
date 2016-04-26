gmx grompp -f em.mdp -o em.tpr -c conf.gro 
gmx trjconv -f conf.gro -o conf.pdb -s em.tpr 
