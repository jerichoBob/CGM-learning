#!/bin/bash

# Define directories to be created
DIRS=("src/utils" "src/models" "data" "tests" "assets" "scripts")

# Check if command line argument is provided
if [ -z "$1" ]
then
    # If not provided, prompt the user for a project name
    echo "Please enter a project name:"
    read ROOT
else
    # If provided, use the argument as the project name
    ROOT="$1"
fi

# Create root directory
mkdir $ROOT

# Create directories
for dir in "${DIRS[@]}"; do
    mkdir -p "$ROOT/$dir"
done

# Create .streamlit directory and config.toml
mkdir "$ROOT/.streamlit"
touch "$ROOT/.streamlit/config.toml"

# Create other files in the root directory
touch "$ROOT/src/main.py"
touch "$ROOT/src/utils.py"
touch "$ROOT/requirements.txt"
touch "$ROOT/.gitignore"
touch "$ROOT/README.md"

# Add standard Streamlit .gitignore content
echo -e "*.py[cod]\n*.swap\n*.bak\n__pycache__/\n*.log\n.history/\n.venv/\n.env\n.DS_Store\n" > "$ROOT/.gitignore"
