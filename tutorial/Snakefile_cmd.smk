
(SAMPLES,) = glob_wildcards("data/{sample}.txt")
print(SAMPLES)

rule all:
    input:
        "output_data/overall_line_count.txt"

rule count_lines:
    input:
        text_file="data/{prefix}.txt"
    output:
        txt_out="output_data/{prefix}_counts.txt"
    shell:
        "wc -l {input.text_file} > {output.txt_out}"


rule summarize_line_count:
    input:
        txt_in=expand("output_data/{sample}_counts.txt",
                      sample=SAMPLES)
    output:
        text_out="output_data/overall_line_count.txt"
    shell:
        "cat {input} > {output.text_out}"