.PHONY: test reports all

test:
	python -m pytest -q

reports:
	cd scripts && python run_time_share_study.py && python run_calibration.py && python run_validation.py

all: test reports
