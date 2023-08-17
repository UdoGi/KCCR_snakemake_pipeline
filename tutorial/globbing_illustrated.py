from snakemake.io import glob_wildcards

wildcards = glob_wildcards("data/{number}-{text}.txt")
for number in wildcards.number:
    print(number)
for text in wildcards.text:
    print(text)
