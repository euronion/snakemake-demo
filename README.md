# Demo for Snakemake (DPG spring meeting 2022, Erlangen, Hacky Hour)

## Installation

*Note:*
> * A Python installation with `snakemake` is required for the demo.
> * Clone or download this repository.
> * Commands shown should be assumed to be executed via termin from within the downloaded repository as working directory.

* It is recommended to install Python via `miniconda`. [Link to instructions](https://docs.conda.io/en/latest/miniconda.html)
* After installing and setting up `miniconda`, create a dedicated environment and install `snakemake` via
    ```
    conda create -f envs/environment.yaml
    ```

    Alternatively you can just install `snakemake` in an existing environment using
    ```
    conda install -c bioconda snakemake
    ```

## Examples

1. Simple rule with output to STDOUT, no input/output files
    ```
    snakemake -call say_hello
    ```

2. Rule to produce a welcoming file containing "Welcome Gießen!" in `results/greetings_giessen.txt` based on the city name as written in `data/giessen.txt`
    ```
    snakemake -call greet_giessen
    ```

3. Rule to produce a welcoming file containing "Welcome Waschmaschine!" in `results/greetings_berlin.txt` based on the contents of `data/berlin.txt`
    ```
    snakemake -call greet_berlin
    ```

4. Run the same command from 3. again:
    ```
    snakemake -call greet_berlin
    ```
    Nothing happens, as the output file already exists and is up-to-date.
    Now change the contents of `data/berlin.txt`.
    Running the same command again will re-run:
    ```
    snakemake -call --reason greet_berlin
    ```
    the workflow and create a new `results/greetings_berlin.txt` file representing the changes made to `data/berlin.txt`.
    The output from `Snakemake` now also shows the `reason`, why the rule was run again and the output file was recreated: The input file was changed and the output file thus no longer up-to-date.

5. Generalisation of `greet_giessen` and `greet_berlin` using `city` as a wildcard. For this to work as intended, change `if True` in `Snakefile` in line 6 to `if False`. The rule to be run is automatically determined based on the output file specified. The output file path also is responsible for setting the `{city}` wildcard to the value `giessen` or `berlin`:
    ```
    snakemake -call results/greetings_giessen.txt results/greetings_berlin.txt
    ```
6. More complex chain of rules reading the files in `data/tickets/`, combining them into `resources/tickets.csv` and then matching the entries there against `data/conference_schedule.csv` to create `results/attendence_schedule.csv`
    ```
    snakemake -call results/attendence_schedule.csv
    ```
7. Plot the DAG (direct acyclic graph) `figures/dag.png` for 6.
    ```
    snakemake --dag results/attendence_schedule.csv | dot -Tpng > figures/dag.png
    ```