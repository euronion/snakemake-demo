import pandas as pd

locations = list()

for city_file in snakemake.input["tickets"]:
    with open (city_file, "r") as f:
        city_name = f.readline().strip()

        locations.append(city_name)

df = pd.DataFrame(locations, columns=["city"])
df.to_csv(snakemake.output["tickets"])