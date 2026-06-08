sif_img="/home/data1/ryanyip/project/RFdiffusion/rosetta_relax.sif"

input_dir=$1  
pdbfiles="${input_dir}/pdbfiles.txt"

find ${input_dir} -iname "*.pdb" -maxdepth 1 > ${pdbfiles}


#  reference : https://blog.csdn.net/m0_55097528/article/details/132350336
echo "rosetta InterfaceAnalyzer"
singularity exec --nv ${sif_img} InterfaceAnalyzer -in:file:l ${pdbfiles}  \
	-use_input_sc \
	--compute_interface_sc true \
	-compute_packstat true \
	-tracer_data_print true \
	-out:file:score_only ${input_dir}/pack_input_score.sc \
	-pack_input true \
	-pack_separated true \
	-add_regular_scores_to_scorefile true \
	-overwrite true \
	-atomic_burial_cutoff 0.01 \
	-sasa_calculator_probe_radius 1.4  \
	-pose_metrics::interface_cutoff 8.0 \
        -fixedchains A B	


cat ${input_dir}/pack_input_score.sc | awk '{s=""; for(i=2;i<=NF;i++){s=s $i","} print substr(s,1,length(s)-1)}' > ${input_dir}/InterfaceAnalyzer_output.csv


