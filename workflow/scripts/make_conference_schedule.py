import pandas as pd

tickets = pd.read_csv(snakemake.input["tickets"])
conferences = pd.read_csv(snakemake.input["conferences"])

schedule = conferences[conferences["City"].isin(tickets["city"])]

schedule.to_csv(snakemake.output["schedule"])