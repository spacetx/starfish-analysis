#!/usr/local/env bash

# extract the starfish version
STARFISH_COMMIT=$(cat $2 | jq -r .commit)

# store current directory
CWD=$(pwd)

# checkout and install the correct starfish version
CD $1 && git checkout ${STARFISH_COMMIT}
source .venv/bin/activate
pip install -e .

# install any additional requirements
ADDITIONAL_REQUIREMENTS=$(cat $2 | jq -r .requirements)

# run the script
SCRIPT=$(cat $2 | jq -r .script_name)
python3 ${SCRIPT}
