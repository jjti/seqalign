edit:
	git ls-files -z | xargs -0 nvim +only -o

test:
	cargo test

matricies:
	cd matricies && \
		wget -r ftp://ftp.ncbi.nlm.nih.gov/blast/matrices/ && \
		mv ftp.ncbi.nlm.nih.gov/blast/blast/matricies/* . && \
		rm -rf ftp.ncbi.nlm.nih.gov

