build:
	cargo build --bin seqalign --release

test:
	cargo test

matrices:
	cd ./src/matrices/data && \
		wget -r ftp://ftp.ncbi.nlm.nih.gov/blast/matrices/ && \
		mv ftp.ncbi.nlm.nih.gov/blast/matrices/* . && \
		rm -rf ftp.ncbi.nlm.nih.gov
	python3 src/matrices/codegen.py
