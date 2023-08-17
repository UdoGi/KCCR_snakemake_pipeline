
wildcards = glob_wildcards("data/{sample}.txt")

print(wildcards.sample)
print(expand("output_data2/{sample}_counts.txt", sample=wildcards.sample))

rule all:
    input:
        expand("output_data2/{sample}_counts.txt", sample=wildcards.sample)

rule count_lines:
    input:
        text_file="data/{sample}.txt"
    output:
        txt_out="output_data2/{sample}_counts.txt"
    shell:
        "wc -l {input.text_file} > {output.txt_out}"

