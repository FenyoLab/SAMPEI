venv: venv/bin/activate

venv/bin/activate: requirements.txt
	test -d venv || python3 -m venv venv
	venv/bin/pip install -Ur requirements.txt
	touch venv/bin/activate

pc_data.csv: data/raw/masterfile_supreme_2.csv src/data/loader.py
	venv/bin/python src/data/loader.py "data/raw/masterfile_supreme_2.csv" "data/interim/pc_data.csv" 

interactive: batch/jupyter_forward.job
	sbatch batch/jupyter_forward.job

token:
	find logs/ -type f -printf '%T@ %p\0' | sort -rz | sed -Ezn '1s/[^ ]* //p' | xargs --null grep 'ssh -N -L\|\]  or http:'

clean:
	rm logs/*.log
