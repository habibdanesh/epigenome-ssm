configfile: "config.json"

SRC = workflow.basedir


rule all:
    input:
        f"{config['results_dir']}/annotate_done.txt"


# Prepare training data
rule prepare_training_data:
    output:
        marker = f"{config['training']['out_dir']}/prepare_done.txt"
    
    threads: config["n_cores"]
    
    shell:
        """
        # Call the child Snakefile for training
        snakemake -s {SRC}/workflow/Snakefile_prepare \
            --cores {threads} \
            --config \
                root_path={SRC} \
                in_files_locator={config[in_files_locator]} \
                regions_file={config[training][regions_file]} \
                chromosomes="{config[training][chromosomes]}" \
                bin_size={config[bin_size]} \
                out_dir={config[training][out_dir]}

        # Touch a marker file to signal that this step is done
        touch {output.marker}
        """


# Prepare annotation data
rule prepare_annotation_data:
    output:
        marker = f"{config['annotation']['out_dir']}/prepare_done.txt"
    
    threads: config["n_cores"]
    
    shell:
        """
        # Call the child Snakefile for annotation
        snakemake -s {SRC}/workflow/Snakefile_prepare \
            --cores {threads} \
            --config \
                root_path={SRC} \
                in_files_locator={config[in_files_locator]} \
                regions_file={config[annotation][regions_file]} \
                chromosomes="{config[annotation][chromosomes]}" \
                bin_size={config[bin_size]} \
                out_dir={config[annotation][out_dir]}
                
        # Touch a marker file to signal that this step is done
        touch {output.marker}
        """


# Run training
rule run_training:
    input:
        marker_tr_prep = f"{config['training']['out_dir']}/prepare_done.txt"
    output:
        marker = f"{config['results_dir']}/train_done.txt"
    threads: config["n_cores"]
    shell:
        """
        snakemake -s {SRC}/workflow/Snakefile_train \
            --cores {threads} \
            --config \
                root_path={SRC} \
                in_dir={config[training][out_dir]} \
                model_type={config[model_type]} \
                n_features={config[n_features]} \
                max_iter={config[max_iter]} \
                out_dir={config[results_dir]} \
                model=None \
                x_array=None \
                debug={config[debug]}
                
        touch {output.marker}
        """


# Run annotation
rule run_annotation:
    input:
        marker_an_prep = f"{config['annotation']['out_dir']}/prepare_done.txt",
        marker_train = f"{config['results_dir']}/train_done.txt"
    output:
        marker = f"{config['results_dir']}/annotate_done.txt"
    threads: config["n_cores"]
    shell:
        """
        snakemake -s {SRC}/workflow/Snakefile_annotate \
            --cores {threads} \
            --config \
                root_path={SRC} \
                in_dir={config[annotation][out_dir]} \
                out_dir={config[results_dir]} \
                model=None \
                n_chunks={config[n_chunks]} \
                ref_genome={config[ref_genome]} \
                debug={config[debug]}
                
        touch {output.marker}
        """