with open(snakemake.input["city"], "r") as f:
    city_name = f.readline().strip()

with open(snakemake.output["city"], "w") as f:
    print(f"Hello {city_name}!", file=f)