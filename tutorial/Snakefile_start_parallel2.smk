
(SAMPLES,) = glob_wildcards("data/{sample}.txt")
print(SAMPLES)
print(expand("output_data/{sample}_counts.txt", sample=SAMPLES))

rule all:
    input:
        expand("output_data/{sample}_counts.txt", sample=SAMPLES)

rule count_lines:
    input:
        text_file="data/{sample}.txt"
    output:
        txt_out="output_data/{sample}_counts.txt"
    shell:
        "wc -l {input.text_file} > {output.txt_out}"

