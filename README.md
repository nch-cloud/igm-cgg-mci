# igm-cgg-mci
Tool for aggregating case data for the Molecular Characterization Initiative into a single summary file.

Reads JSONs in a target folder or folders, identifies valid MCI clinical reports and COG case data files, and outputs a combined Excel sheet or JSON file.

## Usage

Commandline:
Parameter	Description	Required/Optional	Default
--input-json-dirs	Comma-separated list of directories from which to gather JSON files.	Required	-
--output-prefix	Prefix for writing output file(s), including output directory.	Required	-
--data-dict-reference	Tab-separated data dictionary template file.	Optional	Internal file
--blank-field-indicator	Character to indicate fields which were included in data but left blank.	Optional	"." (pass "" to leave empty)
--output-type	Type of output file to produce (Options: Excel, JSON, both).	Optional	Excel

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

The things you need before installing the software.

* You need this
* And you need this
* Oh, and don't forget this

### Installation

A step by step guide that will tell you how to get the development environment up and running.

```
$ First step
$ Another step
$ Final step
```

Available Dockerhub at https://hub.docker.com/r/nationwidechildrens/mci-data-aggregator, 
and on Seven Bridges platform as "MCI Data Aggregator Tool" app.

## Usage

A few examples of useful commands and/or tasks.

```
$ First example
$ Second example
$ And keep this in mind
```

## Deployment

Additional notes on how to deploy this on a live or release system. Explaining the most important branches, what pipelines they trigger and how to update the database (if anything special).

### Server

* Live:
* Release:
* Development:

### Branches

* Master:
* Feature:
* Bugfix:
* etc...

## Additional Documentation and Acknowledgments

* Project folder on server:
* Confluence link:
* Asana board:
* etc...

