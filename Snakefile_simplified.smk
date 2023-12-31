from snakemake.io import glob_wildcards

raw_dir = "/vagrant/ghana_test_files/"
result_folder = "/vagrant/ghana_result/"

(SAMPLES,) = glob_wildcards(raw_dir + "{sample}_R1_001.fastq.gz")
print(SAMPLES)

# link to download a reference fasta file for diamond
reference_data_link = "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz"
# link to download the linking file for megan, linking protein-ids to taxonomic ids
megan_db_link ="https://software-ab.cs.uni-tuebingen.de/download/megan6/megan-map-Feb2022.db.zip"
megan_mapping_fn ='megan-map-Feb2022.db'
reference_data = 'protein_ref.faa'
database_name = result_folder + 'protein_ref_db'
diamond_db_fn =  f"{database_name}.dmnd"

rule all:
    input:
        expand(result_folder + "diamond/{sample}.daa", sample=SAMPLES)

rule fastp:
    input:
        R1=raw_dir + "{sample}_R1_001.fastq.gz",
        R2=raw_dir + "{sample}_R2_001.fastq.gz"
    output:
        R1=result_folder + "fastp/{sample}_R1.fastp.fastq.gz",
        R2=result_folder + "fastp/{sample}_R2.fastp.fastq.gz",
        unpaired=temp(result_folder + 'fastp/{sample}_unpaired.fastp.fastq.gz'),
        merge=temp(result_folder + 'fastp/{sample}_merged.fastp.fastq.gz'),
        singletons=result_folder + 'fastp/{sample}_singletons.fastp.fastq.gz',
        all=result_folder + 'fastp/{sample}_all.fastp.fastq.gz',
        json='fastp' + '/{sample}.fastp.json',
        html='fastp' + '/{sample}.html'
    conda:
        "envs/preprocess.yaml"
    log:
        stdout="logs/fastp/{sample}.stdout.log",
        stderr="logs/fastp/{sample}.stderr.log"
    threads:
        8
    shell:
         "fastp --in1 {input.R1} --in2 {input.R2} --out1 {output.R1} --out2 {output.R2} "
         "--unpaired1 {output.unpaired} --unpaired2 {output.unpaired} --merge "
         "--merged_out {output.merge} --dedup --length_required 30 -p "
         "--json {output.json} --html {output.html} --thread {threads} > {log.stdout} 2> {log.stderr};"
         "touch {output.unpaired} {output.merge};"
         "cat {output.unpaired} {output.merge} > {output.singletons};"
         "cat {output.singletons} {output.R1} {output.R2} > {output.all}"

rule download_protein_ref:
    output:
        fasta=reference_data
    params:
        dest=f"{reference_data}.gz",
        link=reference_data_link
    conda:
        "envs/download_ref.yaml"
    shell:
        'wget {params.link} -O {params.dest}; gunzip {params.dest}'

rule download_megan_mapping_file:
    output:
        fasta=megan_mapping_fn
    params:
        dest=f"{megan_mapping_fn}.zip",
        link=megan_db_link
    conda:
        "envs/download_ref.yaml"
    shell:
        'wget {params.link}; unzip {params.dest}'

rule create_diamond_db:
    input:
        db=reference_data,
    output:
        daa=diamond_db_fn,
    conda:
        "envs/diamond_env.yaml"
    log:
        stdout="logs/diamond/createdb.stdout.log",
        stderr="logs/diamond/createdb.stderr.log"
    params:
        db_name=database_name
    threads:
        32
    shell:
        "diamond makedb --in {input.db} -d {database_name} > {log.stdout} 2> {log.stderr}"

rule diamond:
    input:
        db=diamond_db_fn,
        megan_map=megan_mapping_fn,
        all=result_folder + 'fastp/{sample}_all.fastp.fastq.gz',
    output:
        daa=result_folder + "diamond/{sample}.daa",
    conda:
        "envs/diamond_env.yaml"
    log:
        stdout="logs/diamond/run_diamond_{sample}.stdout.log",
        stderr="logs/diamond/run_diamond_{sample}.stderr.log"
    threads:
        32
    shell:
        "diamond blastx --db {input.db} --query {input.all} --min-orf 1 --threads {threads} "
        "--more-sensitive  -e 0.001 -p {threads} -a {output.daa}> {log.stdout} 2> {log.stderr}; "
        "daa-meganizer -i {output.daa} -mdb {input.megan_map}"


