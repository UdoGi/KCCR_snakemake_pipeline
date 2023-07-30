

rule all:
    input:
        ['output_data/1-karamazov_counts.txt',
         'output_data/2-trial_counts.txt',
         'output_data/3-ulysses_counts.txt']

rule count_lines:
    input:
        text_file="data/{sample}.txt"
    output:
        txt_out="output_data/{sample}_counts.txt"
    shell:
        "wc -l {input.text_file} > {output.txt_out}"

