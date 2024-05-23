export PYTHONPATH := .:$$PYTHONPATH
model:
	python3 main.py
profile:
	./profile/run.sh
debug:
	python3 debug.py
tikz:
	cd tikz && make all
.PHONY: model profile tikz
