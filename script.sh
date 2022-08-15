for i in {1..4}
	do 
		cd system$i
		pwd
		cp ../input_script.vmd .
		sed -i "s/pcca1/pcca$i/g" input_script.vmd
		vmd ../pcca1.pdb -dispdev text -e input_script.vmd
		rm input_script.vmd
	       	cd .. 	
	done
wait
