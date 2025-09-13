# igm-cgg-mci
Tool for aggregating case data for the Molecular Characterization Initiative into a single summary file.

Reads JSONs in a target folder or folders, identifies valid MCI clinical reports and COG case data files, and outputs a combined Excel sheet or JSON file.

## Usage

Run via `python Parse_MCI-JSON.py`

Commandline:
| Parameter	| Description	| Required/Optional	| Default |
| - | - | - | - |
| --input-json-dirs	| Comma-separated list of directories from which to gather JSON files. | Required	| N/A |
| --output-prefix	| Prefix for writing output file(s), including output directory.	| Required	| N/A |
| --data-dict-reference	 | Tab-separated data dictionary template file.	| Optional |	Internal file |
| --blank-field-indicator	 |Character to indicate fields which were included in data but left blank.	| Optional	| "." (pass "" to leave empty) |
| --output-type	 | Type of output file to produce (Options: Excel, JSON, both). |	Optional |	Excel |

## Docker

Available Dockerhub at https://hub.docker.com/r/nationwidechildrens/mci-data-aggregator.
Update Docker repo via:

`docker buildx build --tag nationwidechildrens/mci-data-aggregator:staging --platform=linux/amd64 --sbom=true --provenance=true --push .`


## CGC Seven Bridges
Available on Seven Bridges platform as "MCI Data Aggregator Tool" app.
