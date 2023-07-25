

rule run_bcl2fastq:
    input:
        sample_sheet="/vagrant/230721_KCCR-Run1_Test.csv"
    output:
        out_dir=directory("/vagrant/KCCR_Run1/230721_MN02032_0002_A000H5KF7K/fastq_dir")
    params:
        run_dir="/vagrant/KCCR_Run1/230721_MN02032_0002_A000H5KF7K",
    conda:
        "envs/bcl2fastq.yaml"
    shell:
        "bcl2fastq --sample-sheet {input.sample_sheet} --runfolder-dir {params.run_dir} " 
        "--no-lane-splitting --output-dir {output.out_dir}"