PYTHON ?= /usr/bin/python3
PYTHONPATH := src
MPLCONFIGDIR ?= /tmp/mpl
PUBLIC_DATA := data/public/atract_analysis_public.csv
DICT_PATH := metadata/public_data_dictionary.csv
RAW_PATH ?= $(ATRACT_RAW_XLSX)

.PHONY: public-data analysis all test check

public-data:
	@if [ -z "$(RAW_PATH)" ]; then echo "ATRACT_RAW_XLSX must point to the local raw workbook"; exit 1; fi
	PYTHONPATH=$(PYTHONPATH) MPLCONFIGDIR=$(MPLCONFIGDIR) $(PYTHON) -m atract_analysis public-data --raw "$(RAW_PATH)" --output "$(PUBLIC_DATA)" --dictionary "$(DICT_PATH)"

analysis:
	PYTHONPATH=$(PYTHONPATH) MPLCONFIGDIR=$(MPLCONFIGDIR) $(PYTHON) -m atract_analysis analysis --input "$(PUBLIC_DATA)"

all:
	@if [ -n "$(RAW_PATH)" ]; then \
		PYTHONPATH=$(PYTHONPATH) MPLCONFIGDIR=$(MPLCONFIGDIR) $(PYTHON) -m atract_analysis all --raw "$(RAW_PATH)" --input "$(PUBLIC_DATA)"; \
	else \
		PYTHONPATH=$(PYTHONPATH) MPLCONFIGDIR=$(MPLCONFIGDIR) $(PYTHON) -m atract_analysis analysis --input "$(PUBLIC_DATA)"; \
	fi

check:
	PYTHONPATH=$(PYTHONPATH) MPLCONFIGDIR=$(MPLCONFIGDIR) $(PYTHON) -m atract_analysis check --input "$(PUBLIC_DATA)"

test:
	PYTHONPATH=$(PYTHONPATH) MPLCONFIGDIR=$(MPLCONFIGDIR) $(PYTHON) -m unittest discover -s tests -v
