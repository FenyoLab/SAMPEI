# Agnostic Search

Clone the repository
```shell
git clone https://github.com/FenyoLab/AgnosticSearch.git
```

Navigate to the cloned repository
```shell
cd AgnosticSearch
```

If you have not installed Command Line Tools run this command
```shell
xcode-select --install
```

Create the virtual environment (requires make which comes with Command Line Tools)
```shell
make venv
```

Activate the virtual environment (Do this every time you open a new terminal)
```shell
source venv/bin/activate
```

Run agnostic search
```shell
python -m src.agnostic_search.cli <query_mgf_path> <target_mgf_path> <id_path>
```

Run this for help on the command line arguments for setting the output dir, errors etc
```shell
python -m src.agnostic_search.cli --help
```